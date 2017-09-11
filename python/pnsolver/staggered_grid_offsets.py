
import numpy as np
import util
import rte_terms
import meh
import itertools
import pnsolver
import os
import stencil




# This function is a copy of stencil.eval_term with some modifications.
# The purpose of this function here is to analyse the depencies of
# coefficients l,m on derivatives of other coefficients l+-a, m+-b
# This is what constitues the staggered grid structure and the placements
# of unknowns. This is why we dont do any interpolation when querying
# unknowns, but rather assume and return the unknown location which was
# queried.
def eval_term( expr, info, level=0 ):
	istr = meh.indent_string(level)

	if expr.__class__ == meh.Number:
		if info.debug == True:
			print("{}eval_term::Number".format(istr))
		return expr
	elif expr.__class__ == meh.Negate:
		if info.debug == True:
			print("{}eval_term::Negate".format(istr))
		return meh.neg(eval_term(expr.getExpr(), info, level+1))
	elif isinstance(expr, meh.Variable):
		if info.debug == True:
			print("{}eval_term::Variable {}".format(istr, expr.getSymbol()))

		result = None
		if expr.__class__ == meh.ImaginaryUnit:
			result = expr
		elif expr.getSymbol() in info.vars:
			result = info.vars[expr.getSymbol()]
		else:
			# we simply will pass along unknown variables to the cpp code and
			# expect the user to deal with them
			result = expr
			#raise ValueError("unable to resolve variable {} of type {}".format(expr.getSymbol(), expr.__class__.__name__))

		if info.debug == True:
			print("{}result={}".format(istr, str(result)))
		return result
	elif isinstance(expr, meh.Function):
		if info.debug == True:
			print("{}eval_term::Function".format(istr))
		numArgs = expr.numArguments()
		# evaluate all arguments
		args = []
		for i in range(numArgs):
			arg = eval_term(expr.getArgument(i), info, level+1)
			args.append( arg )
		if expr.getSymbol() == info.unknown_symbol:
			# currently, we assume that:
			# args[0] -> l;args[1] -> m;args[2] -> x
			l_expr = args[0]
			m_expr = args[1]

			if not l_expr.__class__ == meh.Number and not l_expr.__class__ == meh.Negate:
				raise ValueError("expected l to be of expression type meh.Number or meh.Negate")

			if not m_expr.__class__ == meh.Number and not m_expr.__class__ == meh.Negate:
				raise ValueError("expected m to be of expression type meh.Number or meh.Negate")

			# TODO: do we need sure to have ints here?
			l = l_expr.evaluate()
			m = m_expr.evaluate()

			if l < 0 or l > info.gen.order:
				info.term_vanishes = True
				return 0.0

			coeff_index = info.gen.sh_index(l, m)

			# check bounds
			if coeff_index == None or coeff_index >=info.gen.num_coeffs():
				info.term_vanishes = True
				return 0.0


			numDimensions = 2

			# check if the location, at which to evaluate the unknown,
			# matches the actual grid location of the unknown
			# this is true for first order equation with non-anisotropic media
			location_offset = info.location.getOffset()
			#print("\tTEST: {} {}".format(location_offset[0], location_offset[1]))
			u = stencil.Unknown(l, m, info.location.getVoxel())
			u.location_offset = location_offset
			result = stencil.UnknownSet([u])
		else:
			# we have a function to evaluate which is not the unknown
			# does the function has a body?
			if not expr.getBody() is None:
				return meh.num(expr.getBody()(*[arg.evaluate() for arg in args]))

			# kind-of-hack:
			# here we have to check weather the l,m arguments to rte functions, such as q or f_p are valid
			# otherwise we will filter out the term as l,m are invalid
			# this is basially checking the l,m boundaries for functions which depend on them
			if expr.getSymbol() == "q":
				l_expr = args[0]
				m_expr = args[1]

				if not l_expr.__class__ == meh.Number and not l_expr.__class__ == meh.Negate:
					raise ValueError("expected l to be of expression type meh.Number or meh.Negate")

				if not m_expr.__class__ == meh.Number and not m_expr.__class__ == meh.Negate:
					raise ValueError("expected m to be of expression type meh.Number or meh.Negate")

				# TODO: do we need sure to have ints here?
				l = l_expr.evaluate()
				m = m_expr.evaluate()

				if l < 0 or l > info.gen.order or m < -l or m > l:
					info.term_vanishes = True
					return 0.0

			# else keep and return the function expression
			return meh.fun( expr.getSymbol(), *args)

		if info.debug == True:
			print("{}result={}".format(meh.indent_string(level), str(result)))

		#print("test2 {}".format(result.__class__.__name__))
		return result
	elif expr.__class__ == meh.Derivation:
		if info.debug == True:
			print("{}eval_term::Derivation".format(istr))

		if expr.getVariable().getSymbol() == "x":
			# NB: grid has x-dimension along columns (therefore dim=1)
			dimension = 0
		elif expr.getVariable().getSymbol() == "y":
			dimension = 1
		elif expr.getVariable().getSymbol() == "z":
			dimension = 2
		else:
			raise ValueError("unable to identify derivation variable")

		#print("??? dim={}".format(dimension))
		max_dimension = 2
		if dimension >= max_dimension:
			# we have a derivative in z although we only work in 2d domain
			info.term_vanishes = True
			return 0.0	

		# stepsize determines the stepsize of the stencil in number of half-voxels
		stepsize = info.gen.stencil_half_steps
		step = np.zeros(max_dimension, dtype=int)
		step[dimension] = stepsize
		# 
		location = info.location

		#voxelsize = np.array([7.0/70, 7.0/70])
		#central_difference_weight = 1.0/(stepsize*voxelsize[dimension])
		#central_difference_weight = meh.num(1.0/(stepsize*voxelsize[dimension]))
		central_difference_weight = meh.var("h_inv[{}]".format(dimension))

		nested_expr = expr.getExpr()

		#print("??? step={} {}".format(step[0], step[1]))
		#print(location)
		info.location = location.getShiftedLocation(-step)
		#print(info.location)
		#print("TEST: {} {}".format(info.location.getOffset()[0], info.location.getOffset()[1]))
		info.vars["\\vec{x}"] = info.location
		a = eval_term( meh.mul(meh.neg(central_difference_weight), nested_expr.deep_copy()), info, level+1)

		info.location = location.getShiftedLocation(step)
		info.vars["\\vec{x}"] = info.location
		b = eval_term( meh.mul(central_difference_weight, nested_expr.deep_copy()), info, level+1)

		info.location = location
		info.vars["\\vec{x}"] = info.location

		#if a.__class__ == UnknownSet and b.__class__ == UnknownSet:
		#	return a+b
		#elif a.__class__ == UnknownSet or b.__class__ == UnknownSet:
		#	raise ValueError("expecting either a and b or none of both to be an UnknownSet")
		#else:
		result = meh.add(a, b)
		#result = a
		if info.debug == True:
			print("{}result={}".format(meh.indent_string(level), str(result)))
		return result
	elif expr.__class__ == meh.Multiplication:
		if info.debug == True:
			print("{}eval_term::Multiplication".format(meh.indent_string(level)))
		numOperands = expr.numOperands()
		result = []
		for i in range(numOperands):
			result.append(eval_term(expr.getOperand(i), info, level+1))
		if info.debug == True:
			print("{}result={}".format(meh.indent_string(level), str(result)))
		if len(result) <= 1:
			raise ValueError("expected at least two operands")
		return meh.Multiplication(result)
	elif expr.__class__ == meh.Addition:
		if info.debug == True:
			print("{}eval_term::Addition".format(meh.indent_string(level)))
		numOperands = expr.numOperands()
		result = []
		for i in range(numOperands):
			result.append(eval_term(expr.getOperand(i), info, level+1))
		if info.debug == True:
			print("{}result={}".format(meh.indent_string(level), str(result)))
		if len(result) <= 1:
			raise ValueError("expected at least two operands")
		return meh.Addition(result)
	elif expr.__class__ == meh.Power:
		base_expr = eval_term(expr.getBase(), info, level+1)
		exponent_expr = eval_term(expr.getExponent(), info, level+1)
		return meh.Power( base_expr, exponent_expr )
	elif expr.__class__ == meh.Quotient:
		num_expr = eval_term(expr.getNumerator(), info, level+1)
		denom_expr = eval_term(expr.getDenominator(), info, level+1)
		return meh.Quotient( num_expr, denom_expr )
	else:
		raise ValueError("unable to handle expression of type {}".format(expr.__class__.__name__))



def analyse( order, terms ):

	pni = stencil.PNInfo2D(order)
	numCoeffs = pni.num_coeffs()

	todo = [0]

	while todo:
		coeff_index = todo[0]
		(l,m) = pni.lm_index(coeff_index)
		print("coeff={} (l={} m={})==========".format(coeff_index, l, m))
		todo.remove(coeff_index)

		# get the grid position which is associated with the current component/row
		coeff_offset = pni.unknown_offset(coeff_index)
	
		# for each term
		for term_index in range(len(terms)):
		#for term_index in range(1, 2):
			print("\tterm={} ---".format(term_index))
			term = terms[term_index]
			#meh.print_expr(term)

			info = stencil.EvalInfo()
			info.unknown_symbol = "L"
			info.gen = pni
			info.location = stencil.GridLocation2D(np.array([0, 0]), coeff_offset)
			info.vars = {}
			info.vars["\\vec{x}"] = info.location

			# the following snippet will replace l and m with its concrete numbers
			# then it will apply constant folding
			# this is not done for the location variable, as it will change when there
			# are derivatives in the expression
			current_term = term.deep_copy()

			current_term = meh.apply_recursive(current_term, meh.Substitute(meh.var("l'"), meh.num(l)))
			current_term = meh.apply_recursive(current_term, meh.Substitute(meh.var("m'"), meh.num(m)))
			current_term = meh.apply_recursive(current_term, meh.FoldConstants())

			# eval_term will parse the term and replace all occurences of L^{lm} with a special
			# unknown token, holding the information about spatial offset and l,m offset
			# it will also take care of replacing general variables and applying derivation stencils
			# the result is an expression tree where leaves can be UnknownSets
			result = eval_term(current_term, info)

			# terms vanishes in 2d if (l+m)%2 == 1
			# or if l,m indices are out of bounds
			if info.term_vanishes == True:
				#print("term vanishes coeff_index={}".format(coeff_index))
				continue

			# This visitor is applied to the expression tree. It will collapse the unknownsets in Multiplications
			# and Additions. The result is a final unknownset where each unknown weight is an expression tree,
			# or an expression itsself, if no unknowns are present (RHS terms)
			result = meh.apply_recursive(result, stencil.FactorizeUnknowns())
			#print(result)

			if result.__class__ == stencil.UnknownSet:

				print("\tdependency on:".format(coeff_index))

				# the result of term evaluation is a list of unknowns (including voxel and l,m)
				# for each unknown coefficient in A, we have an expression which we can more easily translate
				# into cpp
				u_offset = None
				u_sh_index = None
				for u in result.unknowns:
					
					if u.weight_expr is None:
						expr = meh.num(1.0)
					else:
						expr = u.weight_expr.deep_copy()

					# the weight of each unknown is represented by an expression
					# here we make sure that we collapse numbers 
					expr = meh.apply_recursive(expr, meh.CleanupSigns())
					expr = meh.apply_recursive(expr, meh.FoldConstants())

					# Here we check if the coefficient is reduced to zero. This happens alot.
					if expr.canEvaluate() and np.abs(expr.evaluate()) < 1.0e-8:
						if not u_offset is None:
							raise ValueError("weight expressions need to be vanish for all unknowns")
						# coefficient is vanishes
						continue

					# here we turn negative m into positive m
					# this is done because the step of creating a real valued matrix from
					# the complex-valued coefficients will cause all negative m coefficients
					# to vanish and be folded into coefficients for l, abs(m)
					#if u.m < 0:
					#	u.m = -u.m

					if not u_offset is None:
						if pni.sh_index(u.l, u.m) != u_sh_index:
							raise ValueError("all unknowns expected to have same sh index")
					else:
						u_offset = u.location_offset
						u_sh_index = pni.sh_index(u.l, u.m)

				(u_l,u_m) = pni.lm_index(u_sh_index)
				print("\t\tcoefficient {} (l={} m={})  offset={} {}".format(u_sh_index, u_l, u_m, u_offset[0], u_offset[1]))





if __name__ == "__main__":
	order = 1
	staggered = True
	#filename = "c:/projects/epfl/epfl17/cpp/pnsolver/src/stencil.cpp"








	terms = []
	terms.append(rte_terms.fopn.transport_term())
	terms = rte_terms.splitAddition( terms )

	analyse(order, terms)

