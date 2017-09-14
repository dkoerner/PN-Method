# This module is in charge of generating cpp code for applying the stencil to a particular voxel.
# The generated code is used by the pnsolver (cpp) to generate Ax=b.
# The generator uses terms from the factorized version of the expanded RTE and generates a cpp stencil
# which respects spatial derivatives across voxels and coupling of different SH coefficients within a voxel.


import numpy as np
import util
import rte_terms
import meh
import itertools
import pnsolver
import os





class GridLocation2D(object):
	def __init__(self, voxel, offset):
		voxel_i = voxel[0]
		voxel_j = voxel[1]

		(voxel_offset_i, offset_i) = divmod( offset[0], 2 )
		(voxel_offset_j, offset_j) = divmod( offset[1], 2 )

		self.voxel_i = voxel_i + voxel_offset_i
		self.voxel_j = voxel_j + voxel_offset_j
		self.voxel = np.array([voxel_i + voxel_offset_i,voxel_j + voxel_offset_j])
		self.offset = np.array([offset_i,offset_j])
	def getOffset(self):
		return self.offset
	def getVoxel(self):
		return self.voxel
	def pVS(self):
		return self.voxel + self.offset*0.5
	def getShiftedLocation(self, offset):
		return GridLocation2D(self.voxel, self.offset+offset)
	def __str__(self):
		return "voxel={} {} offset={} {}".format(self.voxel[0], self.voxel[1], self.offset[0], self.offset[1])
	def toLatex(self):
		return self.__class__.__name__
	def deep_copy(self):
		return GridLocation2D(self.voxel, self.offset)




class Unknown(object):
	def __init__(self, coeff_index, voxel, weight_expr = None):
		self.coeff_index = coeff_index
		self.voxel = voxel
		self.weight_expr = weight_expr
	def __str__(self):
		if self.weight_expr is None:
			weight_expr_str = "None"
		elif isinstance(self.weight_expr, meh.Expression):
			weight_expr_str = self.weight_expr.toLatex()
		else:
			weight_expr_str = weight_expr.__class__.__name__
		return "coeff_index={} voxel={} {} weight_expr={}\n".format(self.coeff_index, self.voxel[0], self.voxel[1], weight_expr_str)



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


# eval_term will parse the term and replace all occurences of L^{lm} with a special
# unknown token, holding the information about spatial offset and l,m offset
# it will also take care of replacing general variables and applying derivation stencils
# the result is an expression tree where leaves can be UnknownSets
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


			coeff_index = info.getUnknownIndex(args)

			if coeff_index is None:
				info.term_vanishes = True
				return 0.0

			numDimensions = 2

			# check if the location, at which to evaluate the unknown,
			# matches the actual grid location of the unknown
			# this is true for first order equation with non-anisotropic media
			location_offset = info.location.getOffset()
			unknown_offset = info.pni.getOffset(coeff_index)

			if (location_offset == unknown_offset).all():
				# unknown location and eval location are the same spot
				# no interpolation needed
				u = Unknown(coeff_index, info.location.getVoxel())
				u.interpolated = False
				result = UnknownSet([u])
			elif location_offset[0] == unknown_offset[0]:
				# in the previous if-clause, we checked for location and unknown to be exactly equal
				# now if their offset matches only in one dimension, then we can simply interpolate
				# between the two neighbouring datapoints in that dimension

				# TODO: generalize to 3D

				u0 = Unknown(coeff_index, info.location.getShiftedLocation(np.array([0, 1])).getVoxel(), meh.num(0.5))
				u0.interpolated = True
				u1 = Unknown(coeff_index, info.location.getShiftedLocation(np.array([0, -1])).getVoxel(), meh.num(0.5))
				u1.interpolated = True
				return UnknownSet([u0,u1])
			elif location_offset[1] == unknown_offset[1]:
				u0 = Unknown(coeff_index, info.location.getShiftedLocation(np.array([1, 0])).getVoxel(), meh.num(0.5))
				u0.interpolated = True
				u1 = Unknown(coeff_index, info.location.getShiftedLocation(np.array([-1, 0])).getVoxel(), meh.num(0.5))
				u1.interpolated = True
				result = UnknownSet([u0,u1])
			else:
				# the offsets of evaluation and unknown do not match in any dimension, therefore 
				# we can conclude that the location of the unknown is on the diagonals to the eval location

				# the unknown is located diagonal to the position at which we want to evaluate it
				# therefore we will do an equal weighted sum of all for surrounding unknowns
				offset_combinations = itertools.product(*[[-1, 1] for d in range(numDimensions)])
				num_offset_combinations = 2**numDimensions
				weight =  1.0/num_offset_combinations
				unknowns = []
				for o in offset_combinations:
					u = Unknown(coeff_index, info.location.getShiftedLocation(np.array(o)).getVoxel(), meh.num(weight))
					u.interpolated = True
					unknowns.append(u)
				result = UnknownSet(unknowns)
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

				if l < 0 or l > info.pni.order or m < -l or m > l:
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

		max_dimension = 2
		if dimension >= max_dimension:
			# we have a derivative in z although we only work in 2d domain
			info.term_vanishes = True
			return 0.0	

		# stepsize determines the stepsize of the stencil in number of half-voxels
		stepsize = info.pni.stencil_half_steps
		step = np.zeros(max_dimension, dtype=int)
		step[dimension] = stepsize

		# 
		location = info.location

		#voxelsize = np.array([7.0/70, 7.0/70])
		#central_difference_weight = 1.0/(stepsize*voxelsize[dimension])
		#central_difference_weight = meh.num(1.0/(stepsize*voxelsize[dimension]))
		central_difference_weight = meh.var("h_inv[{}]".format(dimension))

		nested_expr = expr.getExpr()

		info.location = location.getShiftedLocation(-step)
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


def to_cpp( expr, level=0 ):
	istr = meh.indent_string(level)
	#print("{}{}".format(istr, expr.__class__.__name__))
	if expr.__class__ == meh.Number:
		expr_value = expr.evaluate()
		#print(expr_value.__class__.__name__)
		#if expr_value.__class__ == np.complex128 or :
		if np.iscomplex(expr_value):
			return "std::complex<double>({}, {})".format(np.real(expr_value), np.imag(expr_value))
		return str(expr_value)
	elif expr.__class__ == meh.Negate:
		return str( "-{}".format(to_cpp(expr.getExpr(), level+1) ))
	elif expr.__class__ == meh.Variable:
		# this variable could not be resolved at stencil generation time
		# therefore we directly pass it along to the cpp code and
		# expect the user/client code to deal with it
		return expr.getSymbol()
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
	elif expr.__class__ == GridLocation2D:
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


# This visitor is applied to expression trees. It will collapse the unknownsets in Multiplications
# and Additions. The result is a final unknownset where each unknown weight is an expression tree,
# or an expression itsself if no unknowns are present (RHS terms)
class FactorizeUnknowns(object):
	def visit_Negate(self, expr):
		if expr.getExpr().__class__ == UnknownSet:
			return -expr.getExpr()
		return expr
	def visit_Quotient(self, expr):
		num = expr.getNumerator()
		denom = expr.getDenominator()
		return meh.num(num.evaluate()/denom.evaluate())
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






class PNInfo2D(object):
	def __init__(self, order, staggered=True):
		self.order = order

		# the following dictionaries help mapping lm to shindex and vice versa ---
		# NB: sh_index refers to the complex valued SH coefficient
		# find out how much coefficients we have in 2d
		self.index_to_lm = [] # this will map equation index to l,m indices
		self.lm_to_index = {} # this dict will map l,m indices to the equation index
		# NB: we dont use l(l+1)+m because in 2d, we skipp odd (l+m) entries, which changes the sequence.
		# iterate all sh coefficients for given truncation order
		for l in range(0, self.order+1):
			for m in range(-l, l+1):
				# in 2d, we only need to solve for moments where l+m is even
				if (l+m) % 2 == 0:
					self.index_to_lm.append( (l, m) )
					self.lm_to_index[(l,m)] = len(self.index_to_lm)-1
		self.numCoeffs = len(self.index_to_lm)

		# the following dictionaries help mapping lm to uindex and vice versa ---
		# NB: u_index refers to the real-valued coefficient of the solution vector u which
		# represents real and imaginary part (both are real numbers) of sh coefficients lm with m>=0
		self.u_index_to_lm = []
		self.lm_to_u_index = {}
		local_u_index = 0 # this counter represents the current row in the solution vector for a single voxel
		for l in range(0, self.order+1):
			# NB: the starmap code builds the solution vector from reverse order, we do it the same,
			# as we want our conversion matrix S (and therefore the structure of Mx, My) to be indentical
			for m in range(l, -1, -1):

				# In 2d, we only need to solve for moments where l+m is even.
				# This is why all other coefficients are not even part of the global system Ax=b.
				if (l+m) % 2 != 0:
					continue

				# handle real valued part of the current sh coefficient
				key = (l,m,0)
				self.u_index_to_lm.append(key)
				self.lm_to_u_index[key] = local_u_index
				local_u_index += 1

				# handle imaginary part of the current sh coefficient
				# NB: if m==0, then the imaginary part will always be zero, this is why we skip it
				if m>0:
					key = (l,m,1)
					self.u_index_to_lm.append(key)
					self.lm_to_u_index[key] = local_u_index
					local_u_index += 1


		# we keep track of where the unknowns are placed
		self.unknown_info = [ {} for i in range(self.numCoeffs)]

		# by default we place all unknowns at the cell centers
		for i in range(self.numCoeffs):
			self.place_unknown(i, (1, 1))
		# and we use full voxel central differences
		self.stencil_half_steps = 2

		if staggered == True:
			# TODO: generalize this to arbitrary order
			# see starmap python implementation on how to do it
			# but how to do it without M matrix?
			self.place_unknown( 0, (1,1) )
			self.place_unknown( 1, (0,1) )
			self.place_unknown( 2, (1,0) )
			self.stencil_half_steps = 1

	def num_coeffs(self):
		return self.numCoeffs

	def u_index( self, l, m, part ):
		# the third argument identifies the real(0) or imaginary(1) part of the complex number
		key = (l,m,part)
		if key in self.lm_to_u_index:
			return self.lm_to_u_index[key]
		return None

	def lmp_index(self, coeff_index):
		return self.u_index_to_lm[coeff_index]


	'''
	def sh_index(self, l, m):
		key = (l,m)
		#NB: we dont use l(l+1)+m because in 2d, we skipp odd (l+m) entries, which changes the sequence.
		if key in self.lm_to_index:
			return self.lm_to_index[key]
		return None

	def lm_index(self, coeff_index):
		return self.index_to_lm[coeff_index]
	'''


	def place_unknown( self, coeff_index, grid_id ):
		self.unknown_info[coeff_index]['grid_id'] = grid_id
		self.unknown_info[coeff_index]['offset'] = np.array( [grid_id[0], grid_id[1]] , dtype=int)

	def unknown_offset(self, coeff_index):
		return self.unknown_info[coeff_index]['offset']

	def getOffset(self, coeff_index):
		return self.unknown_info[coeff_index]['offset']

	def getLocation(self, voxel, coeff_index):
		return GridLocation2D(voxel, self.unknown_offset(coeff_index)) 



class SubstituteL(object):
	def __init__(self):
		self.success = False

	def visit_Expression(self, expr):
		if expr.__class__ == meh.SHCoefficient and expr.getSymbol() == "L":
			self.success = True
			# we came across an SH coefficient Llm

			# we replace Llm by Rlm and Ilm (see bottom of page 5 in the starmap paper)
			# the real valued coefficients are identified by the function u. It has the same l,m arguments
			# as the Llm it replaces. However, the third argument is an integer which identifies if that
			# particular cofficient relates to the real or imaginary part.
			# therefore the rte terms will be evaluated twice, once for the real part, and once for
			# the imaginary part
			u = meh.fun( "u", expr.getArgument(0), expr.getArgument(1), meh.var("p"), expr.getArgument(2), arglevel=1 )
			u.setLatexArgumentPosition(3, 0)
			return u
		return expr




def generate( complex_valued_terms, order, filename, staggered = True ):
	print("generating {}".format(filename))

	pni = PNInfo2D(order, staggered)

	# this is for checking correctness of the implementation
	#for i in range(pni.num_coeffs()):
	#	(l,m,p) = pni.lmp_index(i)
	#	print( "i={} u: {} l={} m={} check: i={}".format(i, "RI"[p], l, m, pni.u_index(l,m,p)) )


	# source file
	file = open(filename, "w")

	file.write("// This file was generated by {}\n\n".format(__file__))
	file.write("#include <PNSystem.h>\n\n")

	file.write("// truncation order is directly linked to the generated stencil\n")
	file.write("int PNSystem::g_order = {};\n\n".format(order))
	

	arg_sys = "PNSystem::VoxelSystem& sys"
	arg_fields = "PNSystem::Fields& fields"
	prototype = "void set_system_row({},\n\t\t\t\t\t{})\n".format(arg_sys, arg_fields)
	file.write( prototype )
	file.write( "{\n" )
	file.write( "\tV2i vi = sys.getVoxel();\n" )
	file.write( "\tV2d vd = sys.getVoxel().cast<double>();\n" )
	file.write( "\tV2d h_inv( 1.0/({}*sys.getVoxelSize()[0]), 1.0/({}*sys.getVoxelSize()[1]) );\n".format(pni.stencil_half_steps, pni.stencil_half_steps) )

	terms = []
	original_terms = []
	#complex_valued_terms
	#for term_index in range(1):
	for term_index in range(len(complex_valued_terms)):

		#if term_index != 4:
		#	continue

		term = complex_valued_terms[term_index].deep_copy()

		# in each term we replace Llm with Rlm + Ilm which are the real valued coefficients u_0, u_1, ...
		# this is done to create the real valued matrix directly, without requiring a transformation step
		# We replace Llm but have to keep the argument expressions (which could be l+1, etc.).
		substitute = SubstituteL()
		term = meh.apply_recursive(term, substitute)

		#meh.print_expr(term)

		terms.append(term)
		original_terms.append(term_index)

		'''
		if substitute.success == True:
			term = meh.apply_recursive(term, meh.SplitDerivatives())
			term = meh.apply_recursive(term, meh.DistributiveLaw())
			term = meh.apply_recursive(term, meh.CleanupSigns())

			meh.print_expr(term)

			# we should have multiple entries
			if term.__class__ == meh.Addition:
				ops = term.getOperands()
				for op in ops:
					terms.append(op)
					original_terms.append(term_index)
			else:
				raise ValueError( "expected addition after successfull substitution" )
		else:
			terms.append(term)
			original_terms.append(term_index)
		'''

	#print(original_terms[8])

	# for each term
	for term_index in range(len(terms)):
	#for term_index in range(1):
		term = terms[term_index].deep_copy()


		file.write("\n\t// term {} ----------------\n".format(term_index))

		# for each coefficient->l,m
		for coeff_index in range(pni.num_coeffs()):
			(l,m,p) = pni.lmp_index(coeff_index)

			info = EvalInfo()
			#info.debug = True
			info.unknown_symbol = "u"
			info.pni = pni
			info.location = pni.getLocation( np.array([0, 0]), coeff_index )
			info.vars = {}
			# TODO: get rid of info.location alltogether and just use vec{x} ?
			info.vars["\\vec{x}"] = info.location

			# the following snippet will replace l and m with its concrete numbers
			# then it will apply constant folding
			# this is not done for the location variable, as it will change when there
			# are derivatives in the expression
			current_term = term.deep_copy()

			current_term = meh.apply_recursive(current_term, meh.Substitute(meh.var("l'"), meh.num(l)))
			current_term = meh.apply_recursive(current_term, meh.Substitute(meh.var("m'"), meh.num(m)))
			current_term = meh.apply_recursive(current_term, meh.Substitute(meh.var("p"), meh.num(p)))
			current_term = meh.apply_recursive(current_term, meh.FoldConstants())

			# eval_term will parse the term and replace all occurences of u^{lmp} with a special
			# unknown token, holding the information about spatial offset and l,m,p offset
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
			result = meh.apply_recursive(result, FactorizeUnknowns())

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
					if expr.canEvaluate() and np.abs(expr.evaluate()) < 1.0e-8:
						# coefficient is vanishes
						continue

					# the expression consists of function calls or numerical values
					# these are very easily translated into cpp code
					# the t_cpp function parses the expression tree and replaces all meh expression nodes
					# with a string of cpp code
					#print(meh.tree_str(expr))
					expr_cpp = to_cpp(expr)

					file.write( "\tsys.A( {}, vi + V2i({},{}), {} ) += {};".format(coeff_index, u.voxel[0], u.voxel[1], pni.u_index(u.l,u.m,u.p), expr_cpp) )
					file.write("\n")


			else:
				# the result is a simple expression which does not contain any unknown
				# therefore it goes straight to the RHS
				# NB: we do not take care of the sign here. That should be considered
				# by the user
				result_cpp = to_cpp(result)
				file.write( "\tsys.b({}) += {};".format(coeff_index, result_cpp) )
				file.write("\n")

	file.write( "}\n" )

	file.close()



def generate2( order, filename, staggered = True ):
	print("generating {}".format(filename))

	pni = PNInfo2D(order, staggered)

	# source file
	file = open(filename, "w")

	file.write("// This file was generated by {}\n\n".format(__file__))
	file.write("#include <PNSystem.h>\n\n")

	file.write("// truncation order is directly linked to the generated stencil\n")
	file.write("int PNSystem::g_order = {};\n\n".format(order))
	

	arg_sys = "PNSystem::VoxelSystem& sys"
	arg_fields = "PNSystem::Fields& fields"
	prototype = "void set_system_row({},\n\t\t\t\t\t{})\n".format(arg_sys, arg_fields)
	file.write( prototype )
	file.write( "{\n" )
	file.write( "\tV2i vi = sys.getVoxel();\n" )
	file.write( "\tV2d vd = sys.getVoxel().cast<double>();\n" )
	file.write( "\tV2d h_inv( 1.0/({}*sys.getVoxelSize()[0]), 1.0/({}*sys.getVoxelSize()[1]) );\n".format(pni.stencil_half_steps, pni.stencil_half_steps) )
	file.write( "\n" )


	# insert code for setting Mx, Dx, My, Dy, C and q here

	# Currently we handcopied the Mx and My matrices from old simplepn code, which we know is correct.
	# Mx and My do not depend on x (spatial position), nor any RTE variables as the come from
	# discretizing the transport term.
	local_Mx = np.zeros((3, 3))
	local_Mx[1, 0] = 0.57735027
	local_Mx[0, 1] = 0.57735027
	local_My = np.zeros((3, 3))
	local_My[2, 0] = 0.57735027
	local_My[0, 2] = 0.57735027


	file.write( "\n" )


	# use local matrices Mx, My, C and q to build global counterpart
	for i in range(pni.num_coeffs()):
	#for i in range(0, 1):

		location = pni.getLocation( np.array([0,0]), i )
		pVS = location.pVS()


		if i == 0:
			# we exploit the fact that the source is isotropic and therefore has only its 0,0 moment defined
			# usually we would have to apply the transform to the local_q vector: S*local_q in order to
			# transform the coefficients of q from L-space into u-space
			# However, the transform S has no effect on the 0,0 moment
			# we also would have to take the coefficient location into account which we hardcode to voxel center for now
			file.write( "\tsys.coeff_b(0) += fields.q->eval(0,0, sys.voxelToWorld(vd+V2d({}, {}))).real();\n".format(0.5, 0.5) )

			# same is true for the phase function
			file.write( "\tsys.coeff_A(0, vi, 0) += -fields.sigma_s->eval( sys.voxelToWorld(vd+V2d({}, {})) ).real()*fields.f_p->eval( 0, 0, sys.voxelToWorld(vd+V2d({}, {})) ).real();\n".format(pVS[0], pVS[1], pVS[0], pVS[1]) )

		file.write( "\tsys.coeff_A({}, vi, {}) += fields.sigma_t->eval( sys.voxelToWorld(vd+V2d({}, {})) ).real();\n".format(i, i, pVS[0], pVS[1]) )

		#sys.coeff_C( i, vi, i ) = local_C[i,i]
		#sys.coeff_q( i, vi, 0 ) = local_q[i,0]
		for j in range(pni.num_coeffs()):
			local_Ms = [local_Mx, local_My]
			for derivation in range(2):
				# we create the expression M[i, j]*dx(u_j) and let the eval function take care of
				# discretization
				M_ij = meh.fun("M", meh.var("i"), meh.var("j"), body=lambda i,j:local_Ms[derivation][i,j])
				M_ij.setAllSubScripts()
				x = meh.vector( "\\vec{x}", "x", "y", "z" )
				u_j = meh.fun("u", meh.var("j"), x)
				u_j.setLatexArgumentPosition(0, -1)
				current_term = meh.mul( M_ij, meh.deriv( u_j, x.getComponent(derivation), is_partial=True) )


				# replace i, j
				current_term = meh.apply_recursive(current_term, meh.Substitute(meh.var("i"), meh.num(i)))
				current_term = meh.apply_recursive(current_term, meh.Substitute(meh.var("j"), meh.num(j)))
				current_term = meh.apply_recursive(current_term, meh.FoldConstants())

				# now evaluate expression. This is mainly done to apply discretization stencil.
				info = EvalInfo()
				#info.debug = True
				info.unknown_symbol = "u"
				def getUnknownIndex( args ):
					# args[0] -> j (== coefficient index)
					coeff_index_expr = args[0]

					if not coeff_index_expr.__class__ == meh.Number and not coeff_index_expr.__class__ == meh.Negate:
						raise ValueError("expected j to be of expression type meh.Number or meh.Negate")

					# TODO: do we need sure to have ints here?
					coeff_index = coeff_index_expr.evaluate()

					# check bounds
					if coeff_index == None or coeff_index >=info.pni.num_coeffs():
						return None

					return coeff_index

				info.getUnknownIndex = getUnknownIndex
				info.pni = pni
				info.location = location
				info.vars = {}
				# TODO: get rid of info.location alltogether and just use vec{x} ?
				info.vars["\\vec{x}"] = info.location


				# eval_term will parse the term and replace all occurences of u^{lmp} with a special
				# unknown token, holding the information about spatial offset and l,m,p offset
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
				result = meh.apply_recursive(result, FactorizeUnknowns())

				#meh.print_tree(result)
				#print(result)

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
						if expr.canEvaluate() and np.abs(expr.evaluate()) < 1.0e-8:
							# coefficient is vanishes
							continue


						# the expression consists of function calls or numerical values
						# these are very easily translated into cpp code
						# the t_cpp function parses the expression tree and replaces all meh expression nodes
						# with a string of cpp code
						#print(meh.tree_str(expr))
						expr_cpp = to_cpp(expr)

						file.write( "\tsys.coeff_A( {}, vi + V2i({},{}), {} ) += {};\n".format(i, u.voxel[0], u.voxel[1], u.coeff_index, expr_cpp) )
						#file.write("\n")


				else:
					# the result is a simple expression which does not contain any unknown
					# therefore it goes straight to the RHS
					# NB: we do not take care of the sign here. That should be considered
					# by the user
					#result_cpp = to_cpp(result)
					#file.write( "\tsys.b({}) += {};".format(coeff_index, result_cpp) )
					#file.write("\n")
					raise ValueError("rhs?!")





	file.write( "}\n" )

	file.close()



if __name__ == "__main__":
	order = 1
	staggered = True
	filename = "c:/projects/epfl/epfl17/cpp/pnsolver/src/stencil.cpp"

	#terms = []
	'''
	terms.append(rte_terms.lspn.term0_projected_expr())
	terms.append(rte_terms.lspn.term1_projected_expr())
	terms.append(rte_terms.lspn.term2_projected_expr())
	terms.append(rte_terms.lspn.term3_projected_expr())
	terms.append(rte_terms.lspn.term4_projected_expr())
	terms.append(rte_terms.lspn.term5_projected_expr())
	terms.append(rte_terms.lspn.term6_projected_expr())
	'''

	#'''
	#terms.append(rte_terms.fopn.transport_term())
	#terms.append(rte_terms.fopn.collision_term())
	#terms.append(rte_terms.fopn.scattering_term())
	#terms.append(rte_terms.fopn.source_term())
	#'''

	#meh.print_expr(rte_terms.fopn.source_term())
	#meh.print_expr( rte_terms.fopn.transport_term() )


	#terms = rte_terms.splitAddition( terms )


	generate2( order, filename, staggered )
