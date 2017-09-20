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
				return meh.num(0)

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
					return meh.num(0)

			# else keep and return the function expression
			new_fun = meh.fun( expr.getSymbol(), *args)
			for i in range(numArgs):
				new_fun.setLatexArgumentPosition( i, expr.getLatexArgumentPosition(i) )
			return new_fun

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
			return meh.num(0)

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
		term_vanishes = info.term_vanishes
		if info.debug == True:
			print("{}eval_term::Addition".format(meh.indent_string(level)))
		numOperands = expr.numOperands()
		result = []
		for i in range(numOperands):
			info.term_vanishes = False
			evaluated = eval_term(expr.getOperand(i), info, level+1)
			if evaluated.canEvaluate() and np.abs(evaluated.evaluate())<1.0e-8:
				# term is zero, we dont add it to the result list then
				pass
			else:
				result.append(evaluated)
		if info.debug == True:
			print("{}result={}".format(meh.indent_string(level), str(result)))
		info.term_vanishes = term_vanishes
		if len(result) == 0:
			info.term_vanishes = True
			return meh.num(0)
		elif len(result) == 1:
			return result[0]
		else:
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


def to_cpp( expr, info, level=0 ):
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
		return str( "-{}".format(to_cpp(expr.getExpr(), info, level+1) ))
	elif expr.__class__ == meh.Variable:
		# this variable could not be resolved at stencil generation time
		# therefore we directly pass it along to the cpp code and
		# expect the user/client code to deal with it
		return expr.getSymbol()
	elif isinstance(expr, meh.Function):
		return info.fun_to_cpp(expr, info, level)
	elif expr.__class__ == GridLocation2D:
		return "sys.voxelToWorld(vd+V2d({}, {}))".format(expr.getVoxel()[0]+expr.getOffset()[0]*0.5, expr.getVoxel()[1]+expr.getOffset()[1]*0.5 )
	elif expr.__class__ == meh.Derivation:
		raise ValueError("didnt expect any object of type Derivation")
	elif expr.__class__ == meh.Addition:
		numOperands = expr.numOperands()
		cpp_ops = []
		max_str_le = 0
		for i in range(numOperands):
			cpp_op = to_cpp(expr.getOperand(i), info, level+1)
			cpp_ops.append(cpp_op)
			max_str_le = max(max_str_le, len(cpp_op))

		result = "("
		for i in range(numOperands):
			result += cpp_ops[i]
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
			op_cpp = to_cpp(op, info, level+1)
			result += op_cpp
			if i<numOperands-1:
				result+= "*"
		result += ")"
		return result
	elif expr.__class__ == meh.Power:
		base_cpp = to_cpp(expr.getBase(), info)
		exponent_cpp = to_cpp(expr.getExponent(), info)
		result = "std::pow({}, {})".format(base_cpp, exponent_cpp)
		#result = "std::pow({}, 1.0)".format(base_cpp, exponent_cpp)
		return result
	#elif expr.__class__ == meh.Quotient:
	#	num_cpp = to_cpp(expr.getNumerator(), info)
	#	denom_cpp = to_cpp(expr.getDenominator(), info)
	#	result = "({}/{})".format(num_cpp, denom_cpp)
	#	return result
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

		self.build_S()

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


	def build_S(self):
		'''builds the S matrix, which converts from complex-valued to real valued coefficients'''
		# build S matrix ( we iterate over l, m to make sure that the order is correct)

		self.S = np.zeros((self.numCoeffs, self.numCoeffs),dtype=complex)
		count = 0
		for l in range(0, self.order+1):
			for m in range(l, -1, -1):
				# in 2D, we skip coefficients for which l+m is odd
				if (l+m) % 2 != 0:
					continue
				
				# build S matrix, which converts solution from complex to real values

				# computes the real part coefficients for a row (defined by l,m) in the S matrix
				# (see bottom of p.5 in the starmap paper)
				if m == 0:
					self.S[count, self.lm_to_index[(l,m)]] = 1.0
				else:
					self.S[count, self.lm_to_index[(l,m)]] = np.power(-1.0, m)/np.sqrt(2)
					if (l,-m) in self.lm_to_index:
						self.S[count, self.lm_to_index[(l,-m)]] = np.power(-1.0, 2.0*m)/np.sqrt(2)
				count+=1

				# computes the imaginary part coefficients for a row (defined by l,m) in the S matrix
				# (see bottom of p.5 in the starmap paper)
				if m > 0:
					self.S[count, self.lm_to_index[(l,m)]] = np.power(-1.0, m)/np.sqrt(2)*1j
					if (l,-m) in self.lm_to_index:
						self.S[count, self.lm_to_index[(l,-m)]] = -np.power(-1.0, 2*m)/np.sqrt(2)*1j
					count+=1
		self.S_inv = np.linalg.inv(self.S)

	def getS(self):
		return self.S

	def getSInv(self):
		return self.S_inv

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

	def sh_index(self, l, m):
		key = (l,m)
		#NB: we dont use l(l+1)+m because in 2d, we skipp odd (l+m) entries, which changes the sequence.
		if key in self.lm_to_index:
			return self.lm_to_index[key]
		return None

	def lm_index(self, coeff_index):
		return self.index_to_lm[coeff_index]


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



class MatrixExpression(object):
	pass


class UnknownVector(MatrixExpression):
	def __init__(self, expr):
		self.expr = expr
	def getExpr(self):
		return self.expr
	#def setDimensions( shape ):
	#	pass
	def getClassId(self):
		return self.expr.getSymbol()

class DerivationOperator(MatrixExpression):
	def __init__(self, child, variable_symbol):
		self.variable_symbol = variable_symbol
		self.child = child
	#def setDimensions( shape ):
	#	pass
	def getSymbol(self):
		return self.variable_symbol

	def getClassId(self):
		return "d{}{}".format(self.variable_symbol, self.child.getId())

class CoefficientMatrix(MatrixExpression):
	def __init__(self, expr, symbol = "C"):
		self.symbol = symbol
		self.expr = expr

		#self.shape = None
		self.coeff_terms = {}
		#self.coeff_terms = {} # maps coefficient indices (l,m) to expressions
		#self.terms = []

	#def setDimensions( shape ):
	#	self.shape = shape
	#	#self.coeff_terms = [[None for j in shape[1] ] for i in shape[0]]


	def addCoefficientTerm(self, i, j, term ):
		key = (i,j)
		new_term = None
		if not key in self.coeff_terms:
			new_term = term
		else:
			current_term = self.coeff_terms[key]
			if current_term.__class__ != meh.Addition:
				new_term = meh.add(current_term, term)
			else:
				new_term = meh.Addition( current_term.getOperands() + [term] )
		self.coeff_terms[key] = new_term

	def getCoefficientTerm( self, i, j ):
		key = (i,j)
		if not key in self.coeff_terms:
			return None
		return self.coeff_terms[key]

	def derive(self, variable_symbol):
		print("warning: CoefficientMatrix::derive not implemented yet")
		return self
	def getClassId(self):
		return self.symbol
	def __str__(self):
		result = self.getClassId() + "\n"
		for key, expr in self.coeff_terms.items():
			result += "i={} j={} expr=\n{}\n".format(key[0], key[1], meh.latex(expr))
		return result

class MatrixVectorProduct(MatrixExpression):
	def __init__(self, A, B):
		self.A = A
		self.B = B
	def getClassId(self):
		return "{}*{}".format(str(self.A), str(self.B))

class TermClass(object):
	def __init__( self, root ):
		self.root = root
		self.factors = None
		self.pni = None

	def merge( self, other, debug = False ):
		# first as a sanity check make sure factors are identical
		numFactors = len(self.factors)
		if numFactors != len(other.factors):
			raise ValueError("factors should be identical")
		for i in range(numFactors):
			if self.factors[i].__class__ != other.factors[i].__class__:
				raise ValueError("factors should be identical")
			if self.factors[i].__class__ == DerivationOperator and self.factors[i].variable_symbol !=  other.factors[i].variable_symbol:
				raise ValueError("factors should be identical")
		if not other.isLHS() and numFactors != 1:
			raise ValueError("merging RHS term with more than one factor")

		# now iterate all factors and merge coefficient matrix expressions into self
		for i in range(numFactors):
			if self.factors[i].__class__ == CoefficientMatrix:
				# Now we instantiate the unknown and its coefficient term for every possible l,m,
				# therefore iterating over all rows of the coefficient matrix.
				# The unknown index will change and place the term at some column in this matrix
				for coeff_index in range(self.pni.num_coeffs()):

					# find the l, m which are associated with the current coefficient index
					(l,m) = pni.lm_index(coeff_index)

					# we take the term from the other coefficientmatrix
					term = other.factors[i].expr.deep_copy()

					# replace l,m in that given term to instantiate it for the current row in the system
					term = meh.apply_recursive(term, meh.Substitute(meh.var("l'"), meh.num(l)))
					term = meh.apply_recursive(term, meh.Substitute(meh.var("m'"), meh.num(m)))
					term = meh.apply_recursive(term, meh.FoldConstants())


					# if the other term class belongs to the LHS (by having an unknown as a factor),
					# then we derive the column from the coefficient index of that particular unknown
					col = None
					if other.isLHS():
						# We take the unknown with which the other coefficient matrix is associated.
						# This is something like L^(l+1, m-1) and it determines the column at which
						# the term will be placed in the coefficient matrix.
						u_expr = other.getUnknownExpr().deep_copy()

						# replace l,m in the unknown
						u_expr = meh.apply_recursive(u_expr, meh.Substitute(meh.var("l'"), meh.num(l)))
						u_expr = meh.apply_recursive(u_expr, meh.Substitute(meh.var("m'"), meh.num(m)))
						u_expr = meh.apply_recursive(u_expr, meh.FoldConstants())

						# Extract l,m indices from the unknownexpression. This is needed to find the
						# column in the coefficient matrix

						# currently we assume arg[0]->l arg[1]->m
						args = u_expr.getArguments()
						l_expr = args[0]
						m_expr = args[1]

						if not l_expr.canEvaluate():
							raise ValueError("expected expression for l to be able to evaluate")
						if not m_expr.canEvaluate():
							raise ValueError("expected expression for m to be able to evaluate")

						#print( "L l={} m={}".format(coeff_index, col) )

						col = self.pni.sh_index(l_expr.evaluate(), m_expr.evaluate())

						#print( "matrix i={} j={}".format(coeff_index, col) )

						# check bounds
						if col == None or col>=self.pni.num_coeffs():
							# the current term can not be placed within the coefficient matrix
							continue
					else:
						# term is a RHS term, there is just a single column
						col = 0

						# just a quick sanity check



					# finally we can place the term (in which l,m have been replaced) in the matrix
					# if the matrix already contains a term at this position, it will be added.
					if debug == True:
						print("TermClass::merge adding term at {} {}:".format(coeff_index, col))
						meh.print_expr(term)

					self.factors[i].addCoefficientTerm( coeff_index, col, term )





	def getUnknownExpr(self):
		if self.factors is None:
			raise ValueError("must first call finalize")
		for f in self.factors:
			if f.__class__ == UnknownVector:
				return f.getExpr()
		return None

	def isLHS(self):
		if self.getUnknownExpr() is None:
			return False
		return True

	def finalize(self, pni):
		self.factors = self.flatten(self.root)
		self.pni = pni

		if not self.isLHS():
			# sanity check
			if len(self.factors) != 1 or self.factors[0].__class__ != CoefficientMatrix:
				raise ValueError("did expect a single factor for RHS term of type CoefficientMatrix")
			#self.factors[0].setDimensions((pni.num_coeffs(), 1))
		else:
			pass
			#for i in range(numFactors):
			#	if self.factors[i].__class__ == CoefficientMatrix:
			#		self.factors[i].setDimensions((pni.num_coeffs(), pni.num_coeffs()))

	def flatten(self, ct):
		if ct.__class__ == MatrixVectorProduct:
			return self.flatten(ct.A) + self.flatten(ct.B)
		elif ct.__class__ == CoefficientMatrix:
			return [ct]
		elif ct.__class__ == DerivationOperator:
			return [ct] + self.flatten(ct.child)
		elif ct.__class__ == UnknownVector:
			return [ct]
		else:
			raise ValueError("unable to handle class {}".format(ct.__class__.__name__))

	def getId(self):
		if self.factors is None:
			items = self.flatten(self.root)
		else:
			items = self.factors
		result = ""
		for f in items :
			f_str = None
			if f.__class__ == MatrixVectorProduct:
				raise ValueError("unexpected MatrixVectorProduct")
			elif f.__class__ == CoefficientMatrix:
				f_str = f.getClassId()
			elif f.__class__ == DerivationOperator:
				f_str = "d{}".format(f.variable_symbol)
			elif f.__class__ == UnknownVector:
				f_str = f.getClassId()
			else:
				raise ValueError("unable to handle class {}".format(ct.__class__.__name__))

			result += f_str
		return "{}".format(result)

	def __str__(self):
		if self.factors is None:
			items = self.flatten(self.root)
		else:
			items = self.factors

		result = ""
		for item in items:
			result += str(item)


		return result

# classify_term will parse the term and create a string which identifies the term structure
#
def classify_term_recursive( expr, unknown_symbol, level=0 ):
	istr = meh.indent_string(level)

	if expr.__class__ == meh.Number:
		if level == 0:
			return TermClass(CoefficientMatrix(expr))
		else:
			return CoefficientMatrix(expr)
	elif expr.__class__ == meh.Negate:
		ct = classify_term_recursive(expr.getExpr(), unknown_symbol, level+1)
		if ct.__class__ == CoefficientMatrix:
			result = CoefficientMatrix(meh.neg(ct.expr))
		elif ct.__class__ == MatrixVectorProduct:
			if ct.A.__class__ == CoefficientMatrix and ct.B.__class__ != CoefficientMatrix:
				ct.A = CoefficientMatrix(meh.neg(ct.A.expr))
				result = ct
			else:
				raise ValueError("we currently assume that A is a CoefficientMatrix and the unknown vector is always B")
		else:
			result = MatrixVectorProduct( CoefficientMatrix(meh.Number(-1)), ct )
		if level == 0:
			return TermClass(result)
		else:
			return result
	elif isinstance(expr, meh.Variable):
		#print(expr.getSymbol())
		#raise ValueError("work out negate")
		return CoefficientMatrix(expr)
	elif isinstance(expr, meh.Function):
		fun_class = None
		if expr.getSymbol() == unknown_symbol:
			fun_class = UnknownVector(expr)
		else:
			fun_class = CoefficientMatrix( expr )
		if level == 0:
			return TermClass(fun_class)
		else:
			return fun_class

	elif expr.__class__ == meh.Derivation:
		ct = classify_term_recursive(expr.getExpr(), unknown_symbol, level+1)
		if ct.__class__ == CoefficientMatrix:
			return CoefficientMatrix( meh.deriv(ct.expr, expr.getVariable(), is_partial=True)  )
		else:
			# assuming unknownvector or a DerivationOperator aplied to an unknown vector
			return DerivationOperator(ct, expr.getVariable().getSymbol() )
	elif expr.__class__ == meh.Multiplication:
		numOperands = expr.numOperands()
		constant = None
		other = None
		for i in range(numOperands):
			ct = classify_term_recursive(expr.getOperand(i), unknown_symbol, level+1)
			if ct.__class__ == CoefficientMatrix:
				if constant == None:
					constant = ct
				else:
					if ct.expr.__class__ == meh.Multiplication:
						constant.expr = meh.Multiplication( [constant.expr] + ct.expr.getOperands() )
					else:
						constant.expr = meh.mul( constant.expr, ct.expr )

			elif ct.__class__ != CoefficientMatrix and other is None:
				other = ct
			elif ct.__class__ != CoefficientMatrix:
				raise ValueError("expected only one operand which is not a coefficientmatrix")

		result = None
		if other is None and not constant is None:
			result = constant
		elif constant is None and not other is None:
			result = other
		else:
			result = MatrixVectorProduct(constant, other)
		
		if level == 0:
			return TermClass(result)
		else:
			return result
	elif expr.__class__ == meh.Addition:
		numOperands = expr.numOperands()
		result = []
		for i in range(numOperands):
			result.append(eval_term(expr.getOperand(i), unknown_symbol, level+1))
		if len(result) <= 1:
			raise ValueError("expected at least two operands")
		return meh.Addition(result)
	elif expr.__class__ == meh.Power:
		base_ct = classify_term_recursive(expr.getBase(), unknown_symbol, level+1)
		exponent_ct = classify_term_recursive(expr.getExponent(), unknown_symbol, level+1)
		if base_ct.__class__ == CoefficientMatrix and exponent_ct.__class__ == CoefficientMatrix:
			return CoefficientMatrix( meh.pow(base_ct.expr, exponent_ct.expr) )
		else:
			raise ValueError("handling of unknowns or derivatives in nominator or denominator not implemented yet.")
	elif expr.__class__ == meh.Quotient:
		num_ct = classify_term_recursive(expr.getNumerator(), unknown_symbol, level+1)
		denom_ct = classify_term_recursive(expr.getDenominator(), unknown_symbol, level+1)
		if num_ct.__class__ == CoefficientMatrix and denom_ct.__class__ == CoefficientMatrix:
			return CoefficientMatrix( meh.frac(num_ct.expr, denom_ct.expr) )
		else:
			raise ValueError("sadsadsad")
	else:
		raise ValueError("unable to handle expression of type {}".format(expr.__class__.__name__))

def a_lm( l, m ):
	base = (l-m+1)*(l+m+1)/((2*l+3)*(2*l+1))
	if base < 0.0:
		return np.nan
	return np.sqrt(base)
def b_lm( l, m ):
	base = (l-m)*(l+m)/((2*l+1)*(2*l-1))
	if base < 0.0:
		return np.nan
	return np.sqrt(base)
def c_lm( l, m ):
	base = (l+m+1)*(l+m+2)/((2*l+3)*(2*l+1))
	if base < 0.0:
		return np.nan
	return np.sqrt(base)
def d_lm( l, m ):
	base = (l-m)*(l-m-1)/((2*l+1)*(2*l-1))
	if base < 0.0:
		return np.nan
	return np.sqrt(base)
def e_lm( l, m ):
	base = (l-m+1)*(l-m+2)/((2*l+3)*(2*l+1))
	if base < 0.0:
		return np.nan
	return np.sqrt(base)
def f_lm( l, m ):
	base = (l+m)*(l+m-1)/((2*l+1)*(2*l-1))
	if base < 0.0:
		return np.nan
	return np.sqrt(base)





class ReplaceCoefficientMatrixComponents(object):
	def __init__(self, coefficient_matrix_register):
		self.coefficient_matrix_register = coefficient_matrix_register
	def visit_Function( self, expr ):
		if expr.getSymbol() == "M":
			global_index = expr.getArgument(0).evaluate()
			i = expr.getArgument(1).evaluate()
			j = expr.getArgument(2).evaluate()

			f = self.coefficient_matrix_register[global_index]
			if hasattr(f, "M_real"):
				return meh.num(f.M_real[i, j])
			elif f.getCoefficientTerm(i, j) == None:
				# if there is no coefficient term for the given i,j index
				# then we return zero (which allows better pruning)
				return meh.num(0)
			else:
				return expr
		else:
			return expr



def fun_to_cpp( expr, info, level ):
	numArgs = expr.numArguments()
	# turn all arguments into cpp
	arg_str = ""
	for i in range(numArgs):
		arg_str += to_cpp(expr.getArgument(i), info, level+1)
		if i < numArgs-1:
			arg_str += ", "

	symbol = expr.getSymbol()
	if symbol == "M":
		# we reference a coefficient matrix
		# the first argument is the global matrix index
		global_index = expr.getArgument(0).evaluate()
		# second and third argument are the row and column indices
		i = expr.getArgument(1).evaluate()
		j = expr.getArgument(2).evaluate()

		# now get the CoefficientMatrix object behind
		f = info.coefficient_matrix_register[global_index]

		# if the real matrix exists, we can use the constant value diretly
		if hasattr(f, "M_real"):
			return str(f.M_real[i, j])
		else:
			# we assume that a matrix f.id+"_real" has been created in code
			return "{}_real.coeffRef({}, {})".format(f.id, i, j)
	elif symbol == "\\sigma_t":
		return "fields.{}->eval({})".format("sigma_t", arg_str)
	elif symbol == "\\sigma_a":
		return "fields.{}->eval({})".format("sigma_a", arg_str)
	elif symbol == "\\sigma_s":
		return "fields.{}->eval({})".format("sigma_s", arg_str)
	elif symbol == "q":
		return "fields.{}->eval({})".format("q", arg_str)
	elif symbol == "f_p":
		return "fields.{}->eval({})".format("f_p", arg_str)
	else:
		raise ValueError("unable to convert function {} to c++ code.".format(symbol))


if __name__ == "__main__":
	order = 1
	staggered = True
	filename = "c:/projects/epfl/epfl17/cpp/pnsolver/src/stencil.cpp"


	pni = PNInfo2D(order, staggered)



	terms = []
	#'''
	terms.append(rte_terms.lspn.term0_projected_expr()) #check
	terms.append(rte_terms.lspn.term1_projected_expr()) #check
	terms.append(rte_terms.lspn.term2_projected_expr()) #check
	terms.append(rte_terms.lspn.term3_projected_expr()) #check
	terms.append(rte_terms.lspn.term4_projected_expr()) #check
	terms.append(rte_terms.lspn.term5_projected_expr()) #check
	terms.append(rte_terms.lspn.term6_projected_expr()) #check
	#'''

	'''
	terms.append(rte_terms.fopn.transport_term())
	terms.append(rte_terms.fopn.collision_term())
	terms.append(rte_terms.fopn.scattering_term())
	terms.append(rte_terms.fopn.source_term())
	'''

	#meh.print_expr(rte_terms.fopn.source_term())
	#meh.print_expr( rte_terms.fopn.transport_term() )


	terms = rte_terms.splitAddition( terms )


	#generate2( order, filename, staggered )



	'''
	omega = meh.tensor("\\omega", rank=1, dimension=3)
	omega_x = omega.getComponent(0)
	omega_y = omega.getComponent(1)
	omega_z = omega.getComponent(2)

	x = meh.tensor("\\vec{x}", rank=1, dimension=3)
	x.setComponent(0, meh.var("x"))
	x.setComponent(1, meh.var("y"))
	x.setComponent(2, meh.var("z"))
	x.collapsed = True

	u = meh.fun( "L", meh.var("i"), x)
	dx_u = meh.deriv( u, x.getComponent(0), is_partial=True )
	sigma_t = meh.fun("sigma_t")
	dx_sigma_t = meh.deriv( sigma_t, x.getComponent(0), is_partial=True )

	#term = u
	#term = sigma_t
	#term = dx_sigma_t
	#term = meh.num(2)
	#term = meh.mul(meh.num(1), meh.num(2), meh.num(2), u)
	#term = meh.mul(dx_sigma_t, meh.num(3), meh.num(2), u)
	#term = dx_u


	#info = EvalInfo()
	#info.unknown_symbol = "L"
	#term_class = classify_term(term, info)
	#print(term_class)


	'''

	#term = terms[0]

	#meh.print_expr(term)
	#meh.print_tree(term)






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

	file.write( "\tEigen::Matrix<std::complex<double>, {}, {}> S;\n".format(pni.getS().shape[0], pni.getS().shape[1]) )
	for i in range(pni.getS().shape[0]):
		for j in range(pni.getS().shape[1]):
			value = pni.getS()[i,j]
			file.write("\tS.coeffRef({}, {}) = std::complex<double>({}, {});\n".format(i, j, np.real(value), np.imag(value)))
	file.write( "\tEigen::Matrix<std::complex<double>, {}, {}> SInv;\n".format(pni.getS().shape[0], pni.getS().shape[1]) )
	for i in range(pni.getS().shape[0]):
		for j in range(pni.getS().shape[1]):
			value = pni.getSInv()[i,j]
			file.write("\tSInv.coeffRef({}, {}) = std::complex<double>({}, {});\n".format(i, j, np.real(value), np.imag(value)))

	file.write( "\n" )



	term_classes = {}

	#print(len(terms))
	#for term in terms:
	#	meh.print_expr(term)

	for term in terms:
	#for term in [terms[i] for i in range(2)]:
	#for term in [terms[0]]:
		
		#meh.print_expr(term)
		#meh.print_tree(term)
		term_class = classify_term_recursive(term, "L")
		term_class.finalize(pni)

		key = term_class.getId()

		if not key in term_classes:
			term_classes[key] = term_class

		#print("merging in:")
		#meh.print_expr(term_class.factors[0].expr)

		# now accumulate the expression in term_class 
		#print("merging key={} expr={}".format(key, meh.latex(term_class.factors[0].expr)))
		term_classes[key].merge(term_class)

		#print("after merge:")
		#meh.print_expr(term_classes["C"].factors[0].expr)


		#print(term_class)
		#meh.print_expr(term_class.root.A.expr)

	#print(term_classes)
	#meh.print_expr(term_classes["C"].factors[0].getCoefficientTerm(0,0))
	#meh.print_expr(term_classes["C"].factors[0].getCoefficientTerm(1,0))
	#meh.print_expr(term_classes["C"].factors[0].getCoefficientTerm(2,0))
	#exit(1)
	# Now we have factorized the whole discretized RTE into matrix form. The next step is to generate
	# the c++ stencil code for setting the different (complex-valued) matrix components.


	# Produce stencil code for generating complex-valued coefficient matrices =========================

	# To do this, we iterate over all classes and for each class we iterate over all
	# coefficient matrices. If all its terms can be evaluated (and therefore dont depend
	# on any RTE parameter), we will evaluate them immediately and just set the constant value in c++.
	# If they cant be evaluated, we will generate stencil code to construct the
	# complex-valued coefficient matrix per voxel (using RTE parameter fields)
	coefficient_matrix_register = []
	counter = 0
	for key, value in term_classes.items():
		for f in value.factors:
			if f.__class__ == CoefficientMatrix:

				# the variable name, which will be used in generated c++ code
				f.global_index = counter
				f.id = "M_{}".format(counter)
				coefficient_matrix_register.append(f)
				counter += 1

				file.write( "\t//{} =============\n".format(f.id) )

				# find out if the matrix is constant
				is_constant = True
				for key, expr in f.coeff_terms.items():
					if not expr.canEvaluate():
						is_constant = False
						break
						

				f.numRows = pni.num_coeffs()
				f.numCols = None
				if value.isLHS():
					f.numCols = pni.num_coeffs()
				else:
					f.numCols = 1


				if is_constant and f.coeff_terms:

					# evaluate complex valued matrix straight away
					f.M_complex = np.zeros((f.numRows, f.numCols), dtype=complex)
					for key, expr in f.coeff_terms.items():
						f.M_complex[key[0], key[1]] = expr.evaluate()

					# If M_complex is constant and has already been computed, then we perform the transformation straight away
					# and use the resulting constants directly.
					f.M_real = np.real(pni.getS().dot(f.M_complex.dot(pni.getSInv())))
				else:
					matrix_written = False

					# We later will generate c++ code for converting complex-valued matrices to real-valued ones.
					# For 2d, some of these matrices will come out as zero-matrices in which case we
					# can skip the conversion step. This makes our code easier to read and avoids wasting compute.
					skip_conversion_step = True

					# now generate c++ code for generating the coefficient matrix
					# this is done by calling eval_term on each coefficient term (in order to get constant folding
					# and discretization of differential operators) and converting the result to cpp code expressions.
					for i in range(f.numRows):
						for j in range(f.numCols):
							term = f.getCoefficientTerm(i,j)
							if term is None:
								continue

							#meh.print_expr(term)

							# now we run eval_term in order to discretize the differential operators
							info = EvalInfo()
							#info.debug = True
							info.unknown_symbol = "u"
							info.pni = pni
							info.location = pni.getLocation( np.array([0,0]), i )
							info.vars = {}
							# TODO: get rid of info.location alltogether and just use vec{x} ?
							info.vars["\\vec{x}"] = info.location

							#print("-----------------")
							#meh.print_expr(term)
							result = eval_term(term, info)

							# terms vanishes in 2d if (l+m)%2 == 1
							# or if l,m indices are out of bounds
							if info.term_vanishes == True:
								#print("term vanishes at matrix component={} {}".format(i,j))
								continue
							else:
								# At least one term does not vanish and therefore will have to do 
								# the conversion from complex-valued to real-valued matrices
								skip_conversion_step = False

							#meh.print_expr(result)

							# sanity check
							if result.__class__ == UnknownSet:
								raise ValueError( "no unknown expected in coefficient term" )

							# Here we check if the coefficient is reduced to zero. This happens alot.
							if result.canEvaluate() and np.abs(result.evaluate()) < 1.0e-8:
								# coefficient is vanishes
								continue

							#print(result)
							#meh.print_tree(result)
							#meh.print_expr(result)

							# convert to cpp
							info = EvalInfo()
							info.fun_to_cpp = fun_to_cpp
							expr_cpp = to_cpp(result, info)

							if not matrix_written:
								matrix_written = True
								# generate c++ code for building the matrix in our stencil function per voxel
								file.write( "\tEigen::Matrix<std::complex<double>, {}, {}> {};\n".format(f.numRows, f.numCols, f.id) )


							#
							file.write( "\t{}({}, {}) = {};\n".format( f.id, i, j, expr_cpp ) )



					if skip_conversion_step == False:
						# Produce stencil code for converting from complex-valued system to real valued system =================
						# convert the matrix to real valued matrix
						if f.numCols == 1:
							file.write( "\tEigen::Matrix<double, {}, {}> {}_real = (S*{}).real();\n".format(f.numRows, f.numCols, f.id, f.id) )
						else:
							file.write( "\tEigen::Matrix<double, {}, {}> {}_real = (S*{}*SInv).real();\n".format(f.numRows, f.numCols, f.id, f.id) )

				file.write("\n")




	


	# Now we have the system in matrix form where all terms are real-valued. What remains to be done
	# is to discretize the differential operators dx, dy, dz. This is done using our expression evaluation
	# function. It is able to create correct stencils for first and second order derivatives including
	# RTE parameters (or their gradients) as factors.

	# Produce stencil code for computing the final coefficients of the global matrix A =====
	# TODO: do matrix product and iterate over all occuring terms (including discretization). build terms using meh and run eval_term
	# note: we now have u as unknowns

	file.write( "\t// Assembling global system =============\n" )
	# for every term class
	# 
	x = meh.tensor("\\vec{x}", rank=1, dimension=3)
	x.setComponent(0, meh.var("x"))
	x.setComponent(1, meh.var("y"))
	x.setComponent(2, meh.var("z"))
	x.collapsed = True
	for key, tc in term_classes.items():
		print( "term {} ------------".format(key) )

		#if key != "CdyL":
		#	continue

		numFactors = len(tc.factors)
		variable_names = ["i", "j", "k", "l", "m"]
		next_free_variable = 0

		# We generate terms which represent the Matrix-vector products in einstein summation convention.
		# For staggered grids, we need to know to which row of the global system each terms belongs,
		# this is identified in einstein summation convention by the free independent variable which
		# can not be contracted.
		unbound_variable_index = None

		# just make sure the last factor is the unknown if there is one
		if tc.isLHS() and not tc.factors[numFactors-1].__class__ == UnknownVector:
			raise ValueError("currently we expect the unknown to be the last item in a LHS term classification")

		generic_term = None
		# we iterate over all factors backwards and construct a term which we can use
		# to generate code for the whole matrix vector product (including differntial operators)
		for f in reversed(tc.factors):
			if f.__class__ == UnknownVector:
				generic_term = meh.fun( "u", meh.var(variable_names[next_free_variable]) )
				generic_term.setAllSubScripts()
				unbound_variable_index = next_free_variable
				next_free_variable += 1
			elif f.__class__ == CoefficientMatrix:
				#M_ij = meh.fun("M", meh.num(f.global_index), meh.var(variable_names[next_free_variable-1]), meh.var(variable_names[next_free_variable]))
				row_index = meh.var(variable_names[next_free_variable])
				unbound_variable_index = next_free_variable
				if tc.isLHS():
					col_index = meh.var(variable_names[next_free_variable-1])
				else:
					col_index = meh.num(0)
				M_ij = meh.fun("M", meh.num(f.global_index), row_index, col_index )
				M_ij.setAllSubScripts()
				if generic_term is None:
					# happens for RHS terms
					generic_term = M_ij
				else:
					generic_term = meh.mul( M_ij, generic_term )
				next_free_variable += 1
			elif f.__class__ == DerivationOperator:
				component = {"x":0, "y":1, "z":2}[f.getSymbol()]
				generic_term = meh.deriv( generic_term, x.getComponent(component), is_partial=True )

		# The factor is given in einstein summation convention. We can generate all the different
		# terms by iterating over all combinations of free variables. From these individual terms,
		# we generate the code to set the global coefficient matrix.
		#
		# We take into account of purely constant matrices have already been computed and use their
		# values directly, instead of refering to a variable in c++
		numVariables = next_free_variable
		offset_combinations = itertools.product(*[range(pni.num_coeffs()) for d in range(numVariables)])
		for var_values in offset_combinations:

			#if var_values != (1,0):
			#	continue

			# For staggered grids and for assignment in the global matrix A, we need to know
			# to which row of the local(voxel) system each terms belongs to.
			# This also determines the location at which the term is being evaluated for staggered
			# grids.
			coeff_index = var_values[unbound_variable_index]

			# instantiate the term
			term = generic_term.deep_copy()

			# and replace all the indices by their numbers
			for var_index in range(len(var_values)):
				var_value = var_values[var_index]
				term = meh.apply_recursive(term, meh.Substitute(meh.var(variable_names[var_index]), meh.num(var_value)))

			location = pni.getLocation( np.array([0,0]), coeff_index )

			# now we run eval_term in order to discretize the differential operators
			info = EvalInfo()
			#info.debug = True
			info.unknown_symbol = "u"
			def getUnknownIndex( args ):
				# args[0] -> j (== coefficient index)
				coeff_index_expr = args[0]

				if not coeff_index_expr.__class__ == meh.Number and not coeff_index_expr.__class__ == meh.Negate:
					raise ValueError("expected i to be of expression type meh.Number or meh.Negate")

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

			#meh.print_expr(term)

			result = eval_term( term, info )

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

					# here we make sure that we collapse numbers 
					expr = meh.apply_recursive(expr, meh.CleanupSigns())
					expr = meh.apply_recursive(expr, meh.FoldConstants())
					# TODO: explain...
					expr = meh.apply_recursive(expr, ReplaceCoefficientMatrixComponents(coefficient_matrix_register))

					# Here we check if the coefficient is reduced to zero. This happens alot.
					if expr.canEvaluate() and np.abs(expr.evaluate()) < 1.0e-8:
						# coefficient is vanishes
						continue

					#meh.print_expr(expr)

					# the expression consists of function calls or numerical values
					# these are very easily translated into cpp code
					# the to_cpp function parses the expression tree and replaces all meh expression nodes
					# with a string of cpp code
					#print(meh.tree_str(expr))
					info = EvalInfo()
					info.coefficient_matrix_register = coefficient_matrix_register
					info.fun_to_cpp = fun_to_cpp
					expr_cpp = to_cpp(expr, info)

					file.write( "\tsys.coeff_A( {}, vi + V2i({},{}), {} ) += {};\n".format(coeff_index, u.voxel[0], u.voxel[1], u.coeff_index, expr_cpp) )
					#file.write("\n")
			else:
				expr = meh.apply_recursive(result, ReplaceCoefficientMatrixComponents(coefficient_matrix_register))

				# Here we check if the coefficient is reduced to zero. This happens alot.
				if expr.canEvaluate() and np.abs(expr.evaluate()) < 1.0e-8:
					# coefficient is vanishes
					continue

				#meh.print_expr(expr)

				# the expression consists of function calls or numerical values
				# these are very easily translated into cpp code
				# the to_cpp function parses the expression tree and replaces all meh expression nodes
				# with a string of cpp code
				#print(meh.tree_str(expr))
				info = EvalInfo()
				info.coefficient_matrix_register = coefficient_matrix_register
				info.fun_to_cpp = fun_to_cpp
				expr_cpp = to_cpp(expr, info)

				file.write( "\tsys.coeff_b( {} ) += {};\n".format(coeff_index, expr_cpp) )

				
			#meh.print_expr(generic_term)

	file.write( "}\n" )

	file.close()





