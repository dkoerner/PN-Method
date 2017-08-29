# This module represents the application which uses pnbuilder to solve some problem.

import numpy as np
import problems
import pnbuilder
import util
import lspn
import meh
import itertools
import pnbuilder_cpp
import scipy.io



class Unknown(object):
	def __init__(self, l, m, voxel, weight_expr = None):
		self.l = l
		self.m = m
		self.voxel = voxel
		self.weight_expr = weight_expr
	def __str__(self):
		if self.weight_expr is None:
			weight_expr_str = "None"
		elif isinstance(self.weight_expr, meh.Expression):
			weight_expr_str = self.weight_expr.toLatex()
		else:
			weight_expr_str = weight_expr.__class__.__name__
		return "l={} m={} voxel={} {} weight_expr={}\n".format(self.l, self.m, self.voxel[0], self.voxel[1], weight_expr_str)



class UnknownSet(object):
	def __init__(self, unknowns):
		self.unknowns = unknowns
	def __str__(self):
		result = "unknowns:\n"
		for u in self.unknowns:
			result += "\t {}\n".format(str(u))
		return result
	def __mul__(self, other):
		if other.__class__ == UnknownSet:
			raise ValueError("unable to multiply two unknownsets")
		for u in self.unknowns:
			if u.weight_expr is None:
				u.weight_expr = other
			else:
				u.weight_expr = meh.Multiplication([u.weight_expr, other])
		return self
	def __rmul__(self, lhs):
		return self * lhs
	def __add__(self, other):
		return UnknownSet( self.unknowns + other.unknowns )
	def __radd__(self, lhs):
		return self + lhs
	'''
	def __sub__(self, other):
		return self+ (-1.0)*other
	'''
	def __neg__(self):
		return meh.num(-1.0)*self


class EvalInfo(object):
	def __init__(self):
		self.term_vanishes = False
		self.debug = False

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
		# TODO: if we come across the unknown, then return a stencil point
		#if expr.__class__ == meh.ImaginaryUnit:
		#	result = complex(0.0, 1.0)
		if expr.__class__ == meh.ImaginaryUnit:
			result = expr
		elif expr.getSymbol() in info.vars:
			result = info.vars[expr.getSymbol()]
		else:
			raise ValueError("unable to resolve variable {} of type {}".format(expr.getSymbol(), expr.__class__.__name__))

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
			l = l_expr.getValue()
			m = m_expr.getValue()

			if l < 0 or l > info.pnb.N:
				info.term_vanishes = True
				return 0.0

			coeff_index = info.pnb.sh_index(l, m)

			# check bounds
			if coeff_index == None or coeff_index >=info.pnb.num_coeffs():
				info.term_vanishes = True
				return 0.0

			numDimensions = 2

			# check if the location, at which to evaluate the unknown,
			# matches the actual grid location of the unknown
			# this is true for first order equation with non-anisotropic media
			location_offset = info.location.getOffset()
			unknown_offset = info.pnb.unknown_info[coeff_index]['offset']
			if (location_offset == unknown_offset).all():
				#print("check0")
				# unknown location and eval location are the same spot
				# no interpolation needed
				u = Unknown(l, m, info.location.getVoxel())
				result = UnknownSet([u])
			elif location_offset[0] == unknown_offset[0]:
				#print("check1")
				# in the previous if-clause, we checked for location and unknown to be exactly equal
				# now if their offset matches only in one dimension, then we can simply interpolate
				# between the two neighbouring datapoints in that dimension

				# TODO: generalize to 3D

				u0 = Unknown(l, m, info.location.getShiftedLocation(np.array([0, 1])).getVoxel(), meh.num(0.5))
				u1 = Unknown(l, m, info.location.getShiftedLocation(np.array([0, -1])).getVoxel(), meh.num(0.5))
				return UnknownSet([u0,u1])
			elif location_offset[1] == unknown_offset[1]:
				#print("check2")
				u0 = Unknown(l, m, info.location.getShiftedLocation(np.array([1, 0])).getVoxel(), meh.num(0.5))
				u1 = Unknown(l, m, info.location.getShiftedLocation(np.array([-1, 0])).getVoxel(), meh.num(0.5))
				result = UnknownSet([u0,u1])
			else:
				#print("check3")
				# the offsets of evaluation and unknown do not match in any dimension, therefore 
				# we can conclude that the location of the unknown is on the diagonals to the eval location

				# the unknown is located diagonal to the position at which we want to evaluate it
				# therefore we will do an equal weighted sum of all for surrounding unknowns
				offset_combinations = itertools.product(*[[-1, 1] for d in range(numDimensions)])
				num_offset_combinations = 2**numDimensions
				weight =  1.0/num_offset_combinations
				unknowns = []
				for o in offset_combinations:
					u = Unknown(l, m, info.location.getShiftedLocation(np.array(o)).getVoxel(), meh.num(weight))
					unknowns.append(u)
				result = UnknownSet(unknowns)
		else:
			# we have a function to evaluate which is not the unknown
			# does the function has a body?
			if not expr.getBody() is None:
				return meh.num(expr.getBody()(*[arg.getValue() for arg in args]))
			# else keep and return the function expression
			return meh.fun( expr.getSymbol(), *args)
		#elif expr.getSymbol() in info.functions:
		#	# evaluate function
		#	result = info.functions[expr.getSymbol()]( *args )
		#elif not expr.body2 is None:
		#	result = expr.body2(*args)
		#else:
		#	raise ValueError("function {} not defined for evaluation".format(expr.getSymbol()))
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

		if dimension >= info.pnb.domain.voxelSize().shape[0]:
			# we have a derivative in z although we only work in 2d domain
			info.term_vanishes = True
			return 0.0	

		# stepsize determines the stepsize of the stencil in number of half-voxels
		stepsize = info.pnb.stencil_half_steps
		step = np.zeros(info.pnb.domain.voxelSize().shape[0], dtype=int)
		step[dimension] = stepsize

		# 
		location = info.location

		central_difference_weight = 1.0/(stepsize*info.pnb.domain.voxelSize()[dimension])

		nested_expr = expr.getExpr()

		info.location = location.getShiftedLocation(-step)
		info.vars["\\vec{x}"] = info.location
		a = eval_term( meh.mul(meh.num(-central_difference_weight), nested_expr.deep_copy()), info, level+1)

		info.location = location.getShiftedLocation(step)
		info.vars["\\vec{x}"] = info.location
		b = eval_term( meh.mul(meh.num(central_difference_weight), nested_expr.deep_copy()), info, level+1)

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


def to_cpp( expr, level=0 ):
	istr = meh.indent_string(level)
	#print("{}{}".format(istr, expr.__class__.__name__))
	if expr.__class__ == meh.Number:
		expr_value = expr.getValue()
		if expr_value.__class__ == np.complex128:
			return "std::complex<double>({}, {})".format(np.real(expr_value), np.imag(expr_value))
		return str(expr_value)
	if expr.__class__ == meh.Negate:
		return str( "-{}".format(to_cpp(expr.getExpr(), level+1) ))
	elif isinstance(expr, meh.Function):
		numArgs = expr.numArguments()
		# turn all arguments into cpp
		arg_str = ""
		for i in range(numArgs):
			arg_str += to_cpp(expr.getArgument(i), level+1)
			if i < numArgs-1:
				arg_str += ", "

		field_id = expr.getSymbol()
		if field_id == "\\sigma_t":
			field_id = "sigma_t"
		elif field_id == "\\sigma_a":
			field_id = "sigma_a"
		elif field_id == "\\sigma_s":
			field_id = "sigma_s"
		return "fields.{}->eval({})".format(field_id, arg_str)
	#elif expr.__class__ == PositionToken:
	#	return "loc.getPWS() + V2d({}, {})".format(pnb.domain.voxelSize()[0]*expr.offset[0]*0.5, pnb.domain.voxelSize()[1]*expr.offset[1]*0.5);
	elif expr.__class__ == pnbuilder.GridLocation2D:
		return "sys.voxelToWorld(vd+V2d({}, {}))".format(expr.getVoxel()[0]+expr.getOffset()[0]*0.5, expr.getVoxel()[1]+expr.getOffset()[1]*0.5 )
	elif expr.__class__ == meh.Derivation:
		raise ValueError("didnt expect any object of type Derivation")
	elif expr.__class__ == meh.Addition:
		numOperands = expr.numOperands()
		cpp_ops = []
		max_str_le = 0
		for i in range(numOperands):
			cpp_op = to_cpp(expr.getOperand(i), level+1)
			cpp_ops.append(cpp_op)
			max_str_le = max(max_str_le, len(cpp_op))

		result = "("
		for i in range(numOperands):
			result +=  cpp_ops[i]
			if i<numOperands-1:
				result+= "+"
				if max_str_le > 20:
					result += "\n\t\t\t";
		result += ")"
		return result
	elif expr.__class__ == meh.Multiplication:
		numOperands = expr.numOperands()
		result = "("
		for i in range(numOperands):
			op = expr.getOperand(i)
			op_cpp = to_cpp(op, level+1)
			result += op_cpp
			if i<numOperands-1:
				result+= "*"
		result += ")"
		return result
	elif expr.__class__ == meh.Power:
		base_cpp = to_cpp(expr.getBase())
		exponent_cpp = to_cpp(expr.getExponent())
		result = "std::pow({}, {})".format(base_cpp, exponent_cpp)
		#result = "std::pow({}, 1.0)".format(base_cpp, exponent_cpp)
		return result
	else:
		#pass
		raise ValueError("unable to handle expression of type {}".format(expr.__class__.__name__))
		#print("test {}".format(expr.__class__.__name__))



#class PositionToken(object):
#	def __init__(self, location):
#		self.location = location
#	def toLatex(self):
#		return self.__class__.__name__

class FactorizeUnknowns(object):
	def visit_Negate(self, expr):
		if expr.getExpr().__class__ == UnknownSet:
			return -expr.getExpr()
		return expr
	def visit_Quotient(self, expr):
		num = expr.getNumerator()
		denom = expr.getDenominator()
		return meh.num(num.getValue()/denom.getValue())
	def visit_Multiplication(self, expr):
		# look if there are unknowns among the operands
		numOps = expr.numOperands()
		non_unknowns = []
		unknowns = []
		for i in range(numOps):
			op = expr.getOperand(i)
			if op.__class__ == UnknownSet:
				unknowns.append(op)
			else:
				non_unknowns.append(op)
		if len(unknowns) == 0:
			return expr
		elif len(unknowns) == 1:
			if len(non_unknowns) == 0:
				raise ValueError("multiplication with only one factor")
			elif len(non_unknowns) == 1:
				return non_unknowns[0]*unknowns[0]
			else:
				return meh.Multiplication(non_unknowns)*unknowns[0]
		else:
			raise ValueError("unable to multiply unknownsets")
	def visit_Addition(self, expr):
		# look if there are unknowns among the operands
		numOps = expr.numOperands()
		non_unknowns = []
		unknowns = []
		for i in range(numOps):
			op = expr.getOperand(i)
			if op.__class__ == UnknownSet:
				unknowns.append(op)
			else:
				non_unknowns.append(op)
		if len(unknowns) == 0:
			return expr
		elif len(non_unknowns)==0:
			if len(unknowns) < 2:
				raise ValueError("addition with less than two terms")
			result = unknowns[0]
			for i in range(1, len(unknowns)):
				result = result+unknowns[i]
			return result
		else:
			raise ValueError("unable to add unknowns and non_unknowns")






def test_pnsystem():
	problem = problems.checkerboard()
	sys = pnbuilder_cpp.PNSystem(problem["domain_cpp"], problem["order"])

	sys.setField( "sigma_t", problem["sigma_t"] )
	sys.setField( "sigma_a", problem["sigma_a"] )
	sys.setField( "sigma_s", problem["sigma_s"] )

	sys.build()
	#sys.solve()

	A_sys = sys.get_A()


	debugVoxel = np.array([35,35])
	index_i = sys.getGlobalIndex(debugVoxel, 0)



	#print("PNSystem result:")
	#print(A_sys[index_i:index_i+sys.getNumCoefficients(),:])


	# load 
	filename = "C:/projects/epfl/epfl17/python/sopn/system2_{}.mat".format("checkerboard")
	data = scipy.io.loadmat(filename)
	A_complex = data["A"]
	'''
	print("PNBuilder result:")
	for i in range(sys.getNumCoefficients()):
		for j in range(sys.getNumVoxels()*sys.getNumCoefficients()):
			v = A_complex[index_i+i, j]
			if np.abs(v) > 0.0:
				print("  ({}, {}) {}".format(i, j, v))
	'''

	abs_diff = np.abs(A_sys[index_i:index_i+sys.getNumCoefficients(),:]-A_complex[index_i:index_i+sys.getNumCoefficients(),:])
	diff_flatten = abs_diff.flatten()
	mm = diff_flatten.max(axis=1)
	print(mm[0,0])



if __name__ == "__main__":
	test_pnsystem()
	exit(1)



	problem = problems.checkerboard()
	pnb = pnbuilder.PNBuilder(problem["order"], problem["domain"])

	#'''
	# staggered grid (and 3 point stencil)
	if problem["staggered"] == True:
		pnb.place_unknown( 0, (1,1) )
		pnb.place_unknown( 1, (1,0) )
		pnb.place_unknown( 2, (0,1) )
		pnb.set_stencil_half_steps(1)
	#'''

	'''
	pnb.add_terms( lspn.term0_projected_expr() )
	#pnb.add_terms( lspn.term1_projected_expr() )
	pnb.add_terms( lspn.term2_projected_expr() )
	pnb.add_terms( lspn.term3_projected_expr() )
	pnb.add_terms( lspn.term4_projected_expr() )
	pnb.add_terms( lspn.term5_projected_expr() )
	pnb.add_terms( lspn.term6_projected_expr() )

	(A,b) = pnb.build_global( problem )

	util.write_pn_system( pnb, problem, A, b )
	'''

	# terms from lsrte ------------
	#term_gg = lspn.term0_projected_expr().getOperand(1)
	term_gg = lspn.term0_projected_expr()
	#term_gg = meh.add( lspn.term0_projected_expr().getOperand(1), lspn.term0_projected_expr().getOperand(2) )
	#meh.print_expr(term_gg)
	#term = term_gg
	#term = term_gg.getOperand(1)
	#term = lspn.term2_projected_expr()
	terms = []
	#for i in range(0,50):
	for i in range(term_gg.numOperands()):
		terms.append(term_gg.getOperand(i))


	# terms for testing ------------------
	'''
	x = meh.tensor("\\vec{x}", rank=1, dimension=3)
	x.setComponent(0, meh.var("x"))
	x.setComponent(1, meh.var("y"))
	x.setComponent(2, meh.var("z"))
	x.collapsed = True

	sigma_t = meh.fun("\\sigma_t", x)

	#term = meh.SHCoefficient( "L", meh.var("l'"), meh.var("m'"), x )
	#term = meh.mul(meh.num(2.0), meh.SHCoefficient( "L", meh.var("l'"), meh.var("m'"), x ))
	#term = meh.mul(sigma_t, meh.SHCoefficient( "L", meh.var("l'"), meh.var("m'"), x ))
	#term = meh.deriv(sigma_t, x.getComponent(0), is_partial=True)
	#term = meh.deriv(meh.SHCoefficient( "L", meh.var("l'"), meh.var("m'"), x ), x.getComponent(0), is_partial=True)
	#term = meh.deriv(meh.SHCoefficient( "L", meh.var("l'"), meh.var("m'"), x ), x.getComponent(1), is_partial=True)
	#term = meh.mul(meh.deriv(sigma_t, x.getComponent(0), is_partial=True), meh.SHCoefficient( "L", meh.var("l'"), meh.var("m'"), x ))
	#term = meh.add(meh.SHCoefficient( "L", meh.var("l'"), meh.var("m'"), x ), meh.SHCoefficient( "L", meh.var("l'"), meh.sub(meh.var("m'"), meh.num(1)), x ))
	# derivative of different coefficient index
	#term = meh.deriv(meh.SHCoefficient( "L", meh.add(meh.var("l'"), meh.num(1)), meh.var("m'"), x ), x.getComponent(1), is_partial=True)
	#term = meh.deriv(meh.SHCoefficient( "L", meh.add(meh.var("l'"), meh.num(1)), meh.add(meh.var("m'"), meh.num(1)), x ), x.getComponent(0), is_partial=True)
	#term = meh.deriv(meh.SHCoefficient( "L", meh.add(meh.var("l'"), meh.num(1)), meh.add(meh.var("m'"), meh.num(1)), x ), x.getComponent(1), is_partial=True)
	#TODO: multiple derivatives
	#dx_Llm = meh.deriv(meh.SHCoefficient( "L", meh.add(meh.var("l'"), meh.num(1)), meh.add(meh.var("m'"), meh.num(1)), x ), x.getComponent(0), is_partial=True)
	#dx_Llm = meh.deriv(meh.SHCoefficient( "L", meh.var("l'"), meh.sub(meh.var("m'"), meh.num(2)), x ), x.getComponent(0), is_partial=True)
	#term = dx_Llm
	#term = meh.deriv(dx_Llm, x.getComponent(0), is_partial=True)
	#term = meh.deriv(dx_Llm, x.getComponent(1), is_partial=True)
	'''


	'''
	if term.__class__ == meh.Addition:
		numOperands = term.numOperands()
		#print(numOperands)
		#for i in range(numOperands):
		#for i in range(1,2):
		#	terms.append(term.getOperand(i))
		terms.append(term.getOperand(1))
	else:
		terms.append(term)
	'''


	arg_sys = "PNSystem::Voxel& sys"
	arg_fields = "PNSystem::Fields& fields"
	prototype = "void set_system_row({},\n\t\t\t\t\t{})\n".format(arg_sys, arg_fields)

	# source file
	path = "c:/projects/epfl/epfl17/cpp/pnb_cpp/src"
	file = open(path+"/stencil.cpp", "w")

	file.write("#include <PNSystem.h>\n\n")
	file.write( prototype )
	file.write( "{\n" )
	file.write( "\tV2i vi = sys.getVoxel();\n" )
	file.write( "\tV2d vd = sys.getVoxel().cast<double>();\n\n" )


	unknowns = []
	rhs = []
	# for each term
	for term in terms:
		# for each coefficient->l,m
		for coeff_index in range(pnb.num_coeffs()):
		#for coeff_index in range(1):
		#for coeff_index in range(2, 3):
			(l,m) = pnb.lm_index(coeff_index)

			info = EvalInfo()
			#info.debug = True
			info.unknown_symbol = "L"
			info.pnb = pnb
			info.location = pnbuilder.GridLocation2D(pnb.domain, np.array([0, 0]), pnb.unknown_offset(coeff_index))

			info.vars = {}
			info.vars["\\vec{x}"] = info.location

			# the following snippet will replace l and m with its concrete numbers
			# then it will apply constant folding
			# this is not done for the location variable, as it will change when there
			# are derivatives in the expression
			current_term = term.deep_copy()

			if info.debug == True:
				meh.print_expr(current_term)
			current_term = meh.apply_recursive(current_term, meh.Substitute(meh.var("l'"), meh.num(l)))
			current_term = meh.apply_recursive(current_term, meh.Substitute(meh.var("m'"), meh.num(m)))
			current_term = meh.apply_recursive(current_term, meh.FoldConstants())
			if info.debug == True:
				meh.print_expr(current_term)

			# eval_term will parse the term and replace all occurences of L^{lm} with a special
			# unknown token, holding the information about spatial offset and l,m offset
			# it will also take care of replacing general variables and applying derivation stencils
			# the result is an expression tree where leaves can be UnknownSets
			#print(meh.tree_str(current_term))
			result = eval_term(current_term, info)

			#print(meh.tree_str(result))

			# terms vanishes in 2d if (l+m)%2 == 1
			# or if l,m indices are out of bounds
			if info.term_vanishes == True:
				#print("term vanishes coeff_index={}".format(coeff_index))
				continue

			# the next stepp will collapse the unknownsets in Multiplications and Additions
			# the result could be a final unknownset or an expression
			#print(meh.tree_str(result))
			result = meh.apply_recursive(result, FactorizeUnknowns())
			#print("after")
			#print(meh.tree_str(result))
			#'''
			if result.__class__ == UnknownSet:
				# the result of term evaluation is a list of unknowns (including voxel and l,m)
				# for each unknown coefficient in A, we have an expression which we can more easily translate
				# into cpp
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
					if hasattr(expr, "getValue") and np.abs(expr.getValue()) < 1.0e-8:
						# coefficient is vanishes
						continue


					# the expression consists of function calls or numerical values
					# these are very easily translated into cpp code
					# the t_cpp function parses the expression tree and replaces all meh expression nodes
					# with a string of cpp code
					#print(meh.tree_str(expr))
					
					expr_cpp = to_cpp(expr)

					#print(u)

					#print("voxel={} {} unknown_index={}-> {}".format(u.voxel[0], u.voxel[1], pnb.sh_index(u.l,u.m), pnb.global_index(35+u.voxel[0], 35+u.voxel[1], pnb.sh_index(u.l,u.m))))

					file.write( "\tsys.A( {}, vi + V2i({},{}), {} ) += {};".format(coeff_index, u.voxel[0], u.voxel[1], pnb.sh_index(u.l,u.m), expr_cpp) )
					file.write("\n")

					#print(u)

			else:
				# the result is a simple expression which does not contain any unknown
				# therefore it goes straight to the RHS
				# NB: we do not take care of the sign here. That should be considered
				# by the user
				#rhs.append((coeff_index, result))
				result_cpp = to_cpp(result)
				file.write( "\tsys.b({}) += {};".format(coeff_index, result_cpp) )
				file.write("\n")
			#'''




	file.write( "}\n" )

	file.close()



