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


class StaggeredGridLocation(object):
	def __init__( self, voxel = np.array([0,0,0]), offset = np.array([0,0,0]) ):
		voxel_i = voxel[0]
		voxel_j = voxel[1]
		voxel_k = voxel[2]

		(voxel_offset_i, offset_i) = divmod( offset[0], 2 )
		(voxel_offset_j, offset_j) = divmod( offset[1], 2 )
		(voxel_offset_k, offset_k) = divmod( offset[2], 2 )

		self.voxel = np.array([voxel_i + voxel_offset_i,voxel_j + voxel_offset_j,voxel_k + voxel_offset_k])
		self.offset = np.array([offset_i,offset_j,offset_k])
	def getVoxel(self):
		return self.voxel
	def getOffset(self):
		return self.offset
	def isZero(self):
		for i in range(3):
			if self.voxel[i] != 0 and self.offset[i] != 0:
				return False
		return True

	def __add__(self, other):
		return StaggeredGridLocation( self.voxel+other.voxel, self.offset + other.offset )



class GridLocation3D(object):
	def __init__(self, start, shift, test):
		self.start = start
		self.shift = shift
		self.test = test


	def getOffset(self):
		return (self.start+self.shift).getOffset()
	def getVoxel(self):
		return (self.start+self.shift).getVoxel()
	def getShiftedLocation(self, step):
		return GridLocation3D(self.start, self.shift+StaggeredGridLocation(np.array([0,0,0]), step), "check")

	'''
	def pVS(self):
		return self.voxel + self.offset*0.5
	def __str__(self):
		return "voxel={} {} {} offset={} {} {}".format(self.voxel[0], self.voxel[1], self.voxel[2], self.offset[0], self.offset[1], self.offset[2])
	def toLatex(self):
		return self.__class__.__name__
	def deep_copy(self):
		return GridLocation3D(self.voxel, self.offset)
		'''


'''
class GridLocation3D(object):
	def __init__(self, voxel, offset):
		voxel_i = voxel[0]
		voxel_j = voxel[1]
		voxel_k = voxel[2]

		(voxel_offset_i, offset_i) = divmod( offset[0], 2 )
		(voxel_offset_j, offset_j) = divmod( offset[1], 2 )
		(voxel_offset_k, offset_k) = divmod( offset[2], 2 )

		self.voxel_i = voxel_i + voxel_offset_i
		self.voxel_j = voxel_j + voxel_offset_j
		self.voxel_k = voxel_k + voxel_offset_k
		self.voxel = np.array([voxel_i + voxel_offset_i,voxel_j + voxel_offset_j,voxel_k + voxel_offset_k])
		self.offset = np.array([offset_i,offset_j,offset_k])
	def getOffset(self):
		return self.offset
	def getVoxel(self):
		return self.voxel
	def pVS(self):
		return self.voxel + self.offset*0.5
	def getShiftedLocation(self, offset):
		return GridLocation3D(self.voxel, self.offset+offset)
	def __str__(self):
		return "voxel={} {} {} offset={} {} {}".format(self.voxel[0], self.voxel[1], self.voxel[2], self.offset[0], self.offset[1], self.offset[2])
	def toLatex(self):
		return self.__class__.__name__
	def deep_copy(self):
		return GridLocation3D(self.voxel, self.offset)
'''



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
		return "coeff_index={} voxel={} {} {} weight_expr={}\n".format(self.coeff_index, self.voxel[0], self.voxel[1], self.voxel[2], weight_expr_str)



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


			# location_offset is the staggered grid location at which to evaluate the unknown
			location_offset = info.location.getOffset()
			# unknown_offset is the staggered grid location at which the unknown is defined
			unknown_offset = info.pni.getOffset(coeff_index)

			# now find in which dimension the staggered grid locations match
			same_axis = np.array([0 for i in range(3)])
			for i in range(3):
				if location_offset[i] == unknown_offset[i]:
					same_axis[i] = 1

			# number of axes in which evaluation location and coefficient location are equal
			numSameAxes = np.sum(same_axis)

			# now we build a list of local coordinate axes from those dimensions, where the
			# staggered grid locations of unknown and evaluation location didnt match
			all_axes = [ np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1]) ]
			local_axes = []
			for i in range(3):
				if same_axis[i] == 0:
					local_axes.append(all_axes[i])
			numAxes = len(local_axes)

			if numAxes == 0:
				# unknown location and eval location are at the same spot---no interpolation needed
				u = Unknown(coeff_index, info.location.getVoxel())
				u.interpolated = False
				result = UnknownSet([u])
			else:
				# now we interpolate from the neighbous along the defined local axes
				offset_combinations = itertools.product(*[[-1, 1] for d in range(numAxes)])
				num_offset_combinations = 2**numAxes
				weight =  1.0/num_offset_combinations
				unknowns = []
				for o in offset_combinations:
					offset = np.array([0, 0, 0])
					for i in range(numAxes):
						offset += o[i]*local_axes[i]
					u = Unknown(coeff_index, info.location.getShiftedLocation(offset).getVoxel(), meh.num(weight))
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

		if info.pni.is2D() and dimension == 2:
			# we have a derivative in z, although we only work in 2d domain
			info.term_vanishes = True
			return meh.num(0)

		# stepsize determines the stepsize of the stencil in number of half-voxels
		stepsize = info.pni.stencil_half_steps
		step = np.zeros(3, dtype=int)
		step[dimension] = stepsize

		# 
		location = info.location

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

		result = meh.add(a, b)
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
	elif expr.__class__ == GridLocation3D:
		loc = expr
		#return "TODO"
		#return "domain.voxelToWorld(vd+V3d({}, {}, {}))".format(expr.getVoxel()[0]+expr.getOffset()[0]*0.5, expr.getVoxel()[1]+expr.getOffset()[1]*0.5, expr.getVoxel()[2]+expr.getOffset()[2]*0.5 )

		#return "domain.voxelToWorld(vd+V3d({}, {}, {}))".format(expr.getVoxel()[0]+expr.getOffset()[0]*0.5, expr.getVoxel()[1]+expr.getOffset()[1]*0.5, expr.getVoxel()[2]+expr.getOffset()[2]*0.5 )

		#'''
		if loc.shift.isZero():
			grid_index = info.pni.getGridIndex(loc.getOffset())
			start_code = "+ctx.getVoxelSpaceOffsetFromGrid2({})".format(grid_index)
			if hasattr(loc.start, "to_cpp_str"):
				start_code = loc.start.to_cpp_str
			return "domain.voxelToWorld(vd{})".format(start_code)
		else:
			grid_index = info.pni.getGridIndex(loc.start.getOffset())
			shift = loc.shift.getVoxel()+loc.shift.getOffset()*0.5
			return "domain.voxelToWorld(vd+ctx.getVoxelSpaceOffsetFromGrid2({})+V3d({}, {}, {}))".format(grid_index, shift[0], shift[1], shift[2])		#'''
		
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

		# staggered grid locations
		self.grid_offsets = []
		self.grid_offsets.append( (0, 0, 1) )
		self.grid_offsets.append( (1, 0, 1) )
		self.grid_offsets.append( (1, 1, 1) )
		self.grid_offsets.append( (0, 1, 1) )
		self.grid_offsets.append( (0, 0, 0) )
		self.grid_offsets.append( (1, 0, 0) )
		self.grid_offsets.append( (1, 1, 0) )
		self.grid_offsets.append( (0, 1, 0) )
		self.grid_from_offset = {}
		for i in range(len(self.grid_offsets)):
			self.grid_from_offset[self.grid_offsets[i]] = i


		# we keep track of where the unknowns are placed
		self.unknown_info = [ {} for i in range(self.numCoeffs)]

		# by default we place all unknowns at the cell centers
		for i in range(self.numCoeffs):
			self.place_unknown(i, (1, 1, 1))
		# and we use full voxel central differences
		self.stencil_half_steps = 2

		if staggered == True:
			if self.order >= 1:
				#pass
				self.place_unknown( 1, (0, 1, 1) )
				self.place_unknown( 2, (1, 0, 1) )
			#'''
			if self.order >= 2:
				self.place_unknown( 3, (1, 1, 1) )
				self.place_unknown( 5, (1, 1, 1) )
				self.place_unknown( 4, (0, 0, 1) )
			if self.order >= 3:
				self.place_unknown( 7, (1, 0, 1) )
				self.place_unknown( 9, (1, 0, 1) )
				self.place_unknown( 6, (0, 1, 1) )
				self.place_unknown( 8, (0, 1, 1) )
			if self.order >= 4:
				self.place_unknown( 10, (1, 1, 1) )
				self.place_unknown( 12, (1, 1, 1) )
				self.place_unknown( 14, (1, 1, 1) )
				self.place_unknown( 11, (0, 0, 1) )
				self.place_unknown( 13, (0, 0, 1) )
			if self.order >= 5:
				self.place_unknown( 16, (1, 0, 1) )
				self.place_unknown( 18, (1, 0, 1) )
				self.place_unknown( 20, (1, 0, 1) )
				self.place_unknown( 15, (0, 1, 1) )
				self.place_unknown( 17, (0, 1, 1) )
				self.place_unknown( 19, (0, 1, 1) )
			#'''
			if self.order > 5:
				raise ValueError("CHECK!")
			self.stencil_half_steps = 1

	def getGridOffset(self, grid_index):
		return self.grid_offsets[grid_index]

	def getGridIndex(self, offset):
		key = (offset[0], offset[1], offset[2])
		return self.grid_from_offset[key]


	def getGridIndexFromCoeff(self, coeff):
		offset = self.getOffset(coeff)
		return self.getGridIndex(offset)


	def getNumGrids(self):
		return 4

	def is2D(self):
		return True
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

	def to_complex( self, x_real):
	    # use this to convert the solution from complex valued to real valued
	    numVoxels = int(x_real.shape[0]/self.numCoeffs)
	    x_complex = np.zeros( (numVoxels*self.numCoeffs, 1), dtype=complex )
	    for i in range(numVoxels):
	        block_i = i*self.numCoeffs
	        x_complex[block_i:block_i + self.numCoeffs, :] = self.S_inv @ x_real[block_i:block_i + self.numCoeffs, :]

	    return x_complex

	'''
	def to_real(self, x_complex):
		# use this to convert the solution from real valued to complex valued
		numVoxels = self.domain.res_x*self.domain.res_y
		x_real = np.zeros( (numVoxels*self.numCoeffs), dtype=float )
		for voxel_x in range(self.domain.res_x):
			for voxel_y in range(self.domain.res_y):
				block_i = self.get_global_index(voxel_x, voxel_y, 0)
				#if block_i <= 6628 and block_i+self.numCoeffs >= 6628:
				#if block_i == 6627:
				#	print(block_i)
				#	print(x_complex[block_i:block_i + self.numCoeffs])
				#	print(np.real(self.S.dot(x_complex[block_i:block_i + self.numCoeffs])))
				x_real[block_i:block_i + self.numCoeffs] = np.real(self.S.dot(x_complex[block_i:block_i + self.numCoeffs]))
		return x_real
	'''

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
		self.unknown_info[coeff_index]['offset'] = np.array( [grid_id[0], grid_id[1], grid_id[2]] , dtype=int)

	def getOffset(self, coeff_index):
		return self.unknown_info[coeff_index]['offset']
	'''
	def getLocation(self, voxel, coeff_index):
		#return GridLocation3D(voxel, self.unknown_offset(coeff_index)) 
		start = StaggeredGridLocation(voxel, self.getOffset(coeff_index))
		offset = StaggeredGridLocation()
		return GridLocation3D(start, offset, "check") 
	'''



class PNInfo3D(object):
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
				self.index_to_lm.append( (l, m) )
				sh_index = len(self.index_to_lm)-1
				self.lm_to_index[(l,m)] = sh_index
				#print( "check: {} {}".format(sh_index, util.sh_index(l,m)) )
		self.numCoeffs = len(self.index_to_lm)
		#print( "check: {} {}".format(self.numCoeffs, util.num_sh_coeffs(self.order)) )

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


		# staggered grid locations
		self.grid_offsets = []
		self.grid_offsets.append( (0, 0, 1) )
		self.grid_offsets.append( (1, 0, 1) )
		self.grid_offsets.append( (1, 1, 1) )
		self.grid_offsets.append( (0, 1, 1) )
		self.grid_offsets.append( (0, 0, 0) )
		self.grid_offsets.append( (1, 0, 0) )
		self.grid_offsets.append( (1, 1, 0) )
		self.grid_offsets.append( (0, 1, 0) )
		self.grid_from_offset = {}
		for i in range(len(self.grid_offsets)):
			self.grid_from_offset[self.grid_offsets[i]] = i

		# we keep track of where the unknowns are placed
		self.unknown_info = [ {} for i in range(self.numCoeffs)]

		# by default we place all unknowns at the cell centers
		for i in range(self.numCoeffs):
			self.place_unknown(i, (1, 1, 1))
		# and we use full voxel central differences
		self.stencil_half_steps = 2

		if staggered == True:
			if self.order >= 1:
				self.place_unknown( 0, (1, 1, 1) )
				self.place_unknown( 3, (1, 1, 0) )
				self.place_unknown( 2, (1, 0, 1) )
				self.place_unknown( 1, (0, 1, 1) )
			if self.order >= 2:
				self.place_unknown( 4, (1, 1, 1) )
				self.place_unknown( 8, (1, 1, 1) )
				self.place_unknown( 7, (1, 0, 0) )
				self.place_unknown( 6, (0, 1, 0) )
				self.place_unknown( 5, (0, 0, 1) )
			if self.order >= 3:
				self.place_unknown( 11, (1, 1, 0) )
				self.place_unknown( 15, (1, 1, 0) )
				self.place_unknown( 10, (1, 0, 1) )
				self.place_unknown( 14, (1, 0, 1) )
				self.place_unknown( 9, (0, 1, 1) )
				self.place_unknown( 13, (0, 1, 1) )
				self.place_unknown( 12, (0, 0, 0) )
			if self.order >= 4:
				self.place_unknown( 16, (1, 1, 1) )
				self.place_unknown( 20, (1, 1, 1) )
				self.place_unknown( 24, (1, 1, 1) )
				self.place_unknown( 19, (1, 0, 0) )
				self.place_unknown( 23, (1, 0, 0) )
				self.place_unknown( 18, (0, 1, 0) )
				self.place_unknown( 22, (0, 1, 0) )
				self.place_unknown( 17, (0, 0, 1) )
				self.place_unknown( 21, (0, 0, 1) )
			if self.order >= 5:
				self.place_unknown( 27, (1, 1, 0) )
				self.place_unknown( 31, (1, 1, 0) )
				self.place_unknown( 35, (1, 1, 0) )
				self.place_unknown( 26, (1, 0, 1) )
				self.place_unknown( 30, (1, 0, 1) )
				self.place_unknown( 34, (1, 0, 1) )
				self.place_unknown( 25, (0, 1, 1) )
				self.place_unknown( 29, (0, 1, 1) )
				self.place_unknown( 33, (0, 1, 1) )
				self.place_unknown( 28, (0, 0, 0) )
				self.place_unknown( 32, (0, 0, 0) )
			#'''
			if self.order > 5:
				raise ValueError("CHECK!")
			self.stencil_half_steps = 1

	def getGridOffset(self, grid_index):
		return self.grid_offsets[grid_index]

	def getGridIndex(self, offset):
		key = (offset[0], offset[1], offset[2])
		return self.grid_from_offset[key]


	def getGridIndexFromCoeff(self, coeff):
		offset = self.getOffset(coeff)
		return self.getGridIndex(offset)


	def getNumGrids(self):
		return 8

	def is2D(self):
		return False
	def build_S(self):
		'''builds the S matrix, which converts from complex-valued to real valued coefficients'''
		# build S matrix ( we iterate over l, m to make sure that the order is correct)

		self.S = np.zeros((self.numCoeffs, self.numCoeffs),dtype=complex)
		count = 0
		for l in range(0, self.order+1):
			for m in range(l, -1, -1):

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

	def to_complex( self, x_real):
		# use this to convert the solution from complex valued to real valued
		numVoxels = int(x_real.shape[0]/self.numCoeffs)
		x_complex = np.zeros( (numVoxels*self.numCoeffs, 1), dtype=complex )
		for i in range(numVoxels):
			block_i = i*self.numCoeffs
			x_complex[block_i:block_i + self.numCoeffs, :] = self.S_inv @ x_real[block_i:block_i + self.numCoeffs, :]

		return x_complex

	'''
	def to_real(self, x_complex):
		# use this to convert the solution from real valued to complex valued
		numVoxels = self.domain.res_x*self.domain.res_y*self.domain.res_z
		x_real = np.zeros( (numVoxels*self.numCoeffs), dtype=float )
		for voxel_x in range(self.domain.res_x):
			for voxel_y in range(self.domain.res_y):
				block_i = self.get_global_index(voxel_x, voxel_y, 0)
				#if block_i <= 6628 and block_i+self.numCoeffs >= 6628:
				#if block_i == 6627:
				#	print(block_i)
				#	print(x_complex[block_i:block_i + self.numCoeffs])
				#	print(np.real(self.S.dot(x_complex[block_i:block_i + self.numCoeffs])))
				x_real[block_i:block_i + self.numCoeffs] = np.real(self.S.dot(x_complex[block_i:block_i + self.numCoeffs]))
		return x_real
	'''

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
		self.unknown_info[coeff_index]['offset'] = np.array( [grid_id[0], grid_id[1], grid_id[2]] , dtype=int)

	def getOffset(self, coeff_index):
		return self.unknown_info[coeff_index]['offset']

	'''
	def getLocation(self, voxel, coeff_index):
		#return GridLocation3D(voxel, self.unknown_offset(coeff_index)) 
		start = StaggeredGridLocation(voxel, self.getOffset(coeff_index))
		offset = StaggeredGridLocation()
		return GridLocation3D(start, offset, "check") 
	'''




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





class MatrixExpression(object):
	def __init__(self, id_str = ""):
		self.id = id_str
	def setId(self, id_str):
		self.id = id_str
	def getId( self ):
		return self.id



class UnknownVector(MatrixExpression):
	def __init__(self, expr):
		self.expr = expr
		super(UnknownVector, self).__init__()
	def getExpr(self):
		return self.expr
	def getId( self ):
		return self.expr.getSymbol()
	def getClassId(self):
		return self.expr.getSymbol()

class DerivationOperator(MatrixExpression):
	def __init__(self, child, variable_symbol):
		self.variable_symbol = variable_symbol
		self.child = child
		super(DerivationOperator, self).__init__()
	#def setDimensions( shape ):
	#	pass
	def getSymbol(self):
		return self.variable_symbol
	def getId( self ):
		return "d{}".format(self.variable_symbol)
	def getClassId(self):
		return "d{}{}".format(self.variable_symbol, self.child.getId())

class CoefficientMatrix(MatrixExpression):
	def __init__(self, expr, symbol = "C"):
		self.symbol = symbol
		self.expr = expr
		self.coeff_terms = {}
		super(CoefficientMatrix, self).__init__(symbol)


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
		super(MatrixVectorProduct, self).__init__()

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
					(l,m) = self.pni.lm_index(coeff_index)

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
				f_str = f.getId()
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
		#return "fields.{}->eval({})".format("sigma_t", arg_str)
		return "problem.evalExtinction({})[color_channel]".format(arg_str)
	elif symbol == "\\sigma_a":
		#return "fields.{}->eval({})".format("sigma_a", arg_str)
		return "problem.evalAbsorption({})[color_channel]".format(arg_str)
	elif symbol == "\\sigma_s":
		#return "fields.{}->eval({})".format("sigma_s", arg_str)
		return "problem.evalScattering({})[color_channel]".format(arg_str)
	elif symbol == "q":
		#return "fields.{}->eval({})".format("q", arg_str)
		return "problem.evalEmission({})[color_channel]".format(arg_str)
	elif symbol == "f_p":
		#return "fields.{}->eval({})".format("f_p", arg_str)
		return "problem.evalPhase({})[color_channel]".format(arg_str)

	else:
		raise ValueError("unable to convert function {} to c++ code.".format(symbol))


def generate_stencil_code( stencil_name, filename, terms, pni ):

	print( "generating stencil code {} order={} staggered={}\nfilename={}".format(stencil_name, order, staggered, filename) )

	# the width of the stencil (the number of neighbouring voxels the stencil touches)
	# this values is determined later during construction the code for building A
	stencil_width = 0

	# source file
	file = open(filename, "w")

	file.write("// This file was generated by {}\n\n".format(__file__))
	file.write("#include <PNSystem.h>\n\n")

	#file.write("// truncation order is directly linked to the generated stencil\n")
	#file.write("int PNSystem::g_order = {};\n\n".format(order))
	

	arg_sys = "PNSystem::Stencil::Context& ctx"
	#arg_voxel = "const V2i& voxel"
	#prototype = "void {}({},\n\t\t\t\t\t{})\n".format(stencil_name, arg_sys, arg_voxel)
	prototype = "void {}({})\n".format(stencil_name, arg_sys)
	file.write( prototype )
	file.write( "{\n" )
	file.write( "\tV3i vi = ctx.getVoxelCoord();\n" )
	file.write( "\tV3d vd = vi.cast<double>();\n" )
	file.write( "\tconst Domain& domain = ctx.getDomain();\n" )
	file.write( "\tconst PNVolume& problem = ctx.getProblem();\n" )
	file.write( "\tV3d h_inv( 1.0/({}*domain.getVoxelSize()[0]), 1.0/({}*domain.getVoxelSize()[1]), 1.0/({}*domain.getVoxelSize()[2]) );\n".format(pni.stencil_half_steps, pni.stencil_half_steps, pni.stencil_half_steps) )
	file.write( "\tint color_channel = 0;\n" )
	file.write( "\n" )

	file.write( "\tEigen::Matrix<std::complex<double>, {}, {}> S;\n".format(pni.getS().shape[0], pni.getS().shape[1]) )
	for i in range(pni.getS().shape[0]):
		for j in range(pni.getS().shape[1]):
			value = pni.getS()[i,j]
			if np.abs(value) > 1.0e-8:
				file.write("\tS.coeffRef({}, {}) = std::complex<double>({}, {});\n".format(i, j, np.real(value), np.imag(value)))
	file.write( "\tEigen::Matrix<std::complex<double>, {}, {}> SInv;\n".format(pni.getS().shape[0], pni.getS().shape[1]) )
	for i in range(pni.getS().shape[0]):
		for j in range(pni.getS().shape[1]):
			value = pni.getSInv()[i,j]
			if np.abs(value) > 1.0e-8:
				file.write("\tSInv.coeffRef({}, {}) = std::complex<double>({}, {});\n".format(i, j, np.real(value), np.imag(value)))

	file.write( "\n" )



	term_classes = {}


	for term in terms:
		term_class = classify_term_recursive(term, "L")
		term_class.finalize(pni)

		key = term_class.getId()

		if not key in term_classes:
			term_classes[key] = term_class

		# now accumulate the expression in term_class 
		term_classes[key].merge(term_class)

	# assign the ids to coefficient matrices ----------------
	equ_str = ""
	coefficient_matrix_register = []
	counter = 0
	rhs = None
	for key, value in term_classes.items():
		if equ_str and value.isLHS():
			equ_str += " + "
		for f in value.factors:
			if f.__class__ == CoefficientMatrix:
				# the variable name, which will be used in generated c++ code
				f.global_index = counter
				if value.isLHS():
					f.setId("M_{}".format(counter))
				else:
					f.setId("b")
					if not len(value.factors) == 1:
						raise ValueError("expected only a single coefficient matrix for RHS term")
					if not rhs is None:
						raise ValueError("expected only a single RHS coefficient matrix")
					rhs = value
				coefficient_matrix_register.append(f)
				counter += 1
			if value.isLHS():
				equ_str += "{}".format(f.getId())
	if rhs:
		equ_str += " = {}".format(rhs.getId())

	# Produce stencil code for generating complex-valued coefficient matrices =========================

	# To do this, we iterate over all classes and for each class we iterate over all
	# coefficient matrices. If all its terms can be evaluated (and therefore dont depend
	# on any RTE parameter), we will evaluate them immediately and just set the constant value in c++.
	# If they cant be evaluated, we will generate stencil code to construct the
	# complex-valued coefficient matrix per voxel (using RTE parameter fields)

	file.write( "\t//Producing complex-valued matrices =============\n" )
	file.write( "\t//{}\n".format(equ_str) )
	file.write( "\n" )

	for key, value in term_classes.items():
		for f in value.factors:
			if f.__class__ == CoefficientMatrix:

				file.write( "\t//{} ---\n".format(f.id) )

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
					# NB: I hope we are right in the assumption that for constant
					# coefficient matrices (which dont depend on any RTE variable or position)
					# it will be the same for every grid point
					f.M_complex = np.zeros((f.numRows, f.numCols), dtype=complex)
					for key, expr in f.coeff_terms.items():
						f.M_complex[key[0], key[1]] = expr.evaluate()

					# If M_complex is constant and has already been computed, then we perform the transformation straight away
					# and use the resulting constants directly.
					f.M_real = np.real(pni.getS().dot(f.M_complex.dot(pni.getSInv())))

					file.write( "\t// is constant\n" )
					#file.write( "\tEigen::Matrix<double, {}, {}> {}_real;\n".format(f.numRows, f.numCols, f.id) )
					#for i in range(f.numRows):
					#	for j in range(f.numCols):
					#		file.write( "\t{}_real({}, {}) = {};\n".format( f.id, i, j, f.M_real[i, j] ) )
				else:

					# now generate c++ code for generating the coefficient matrix
					# this is done by calling eval_term on each coefficient term (in order to get constant folding
					# and discretization of differential operators) and converting the result to cpp code expressions.
					code = ""
					for i in range(f.numRows):
						for j in range(f.numCols):
							term = f.getCoefficientTerm(i,j)
							if term is None:
								continue

							# now we run eval_term in order to discretize the differential operators
							info = EvalInfo()
							#info.debug = True
							info.unknown_symbol = "u"
							info.pni = pni
							#info.location = pni.getLocation( np.array([0,0,0]), i )
							info.location = GridLocation3D(StaggeredGridLocation(np.array([0,0,0]), pni.getOffset(i)), StaggeredGridLocation(), "check") 
							info.location.start.to_cpp_str = "+ctx.getVoxelSpaceOffsetFromGrid2(i)"

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


							# sanity check
							if result.__class__ == UnknownSet:
								raise ValueError( "no unknown expected in coefficient term" )

							# Here we check if the coefficient is reduced to zero. This happens alot.
							if result.canEvaluate() and np.abs(result.evaluate()) < 1.0e-8:
								# coefficient is vanishes
								continue

							info = EvalInfo()
							info.fun_to_cpp = fun_to_cpp
							info.pni = pni
							expr_cpp = to_cpp(result, info)

							#
							#file.write( "\t{}({}, {}) = {};\n".format( f.id, i, j, expr_cpp ) )
							#code += "\t{}({}, {}) = {};\n".format( f.id, i, j, expr_cpp )
							code += "\t\t{}({}, {}) = {};\n".format( f.id, i, j, expr_cpp )



					if code:
						# generate c++ code for building the matrix in our stencil function per voxel						
						#file.write( "\tEigen::Matrix<std::complex<double>, {}, {}> {};\n".format(f.numRows, f.numCols, f.id) )
						matrix_def = "\tEigen::Matrix<double, {}, {}> {}_real_staggered[{}];\n".format(f.numRows, f.numCols, f.id, pni.getNumGrids())

						loop = "\tfor( int i=0;i<{};++i )\n".format(pni.getNumGrids())
						loop += "\t{\n"
						loop += "\t\tEigen::Matrix<std::complex<double>, {}, {}> {};\n".format(f.numRows, f.numCols, f.id)

						code = matrix_def + loop + code

						# Produce stencil code for converting from complex-valued system to real valued system =================
						# convert the matrix to real valued matrix
						if not value.isLHS():
							#file.write( "\tEigen::Matrix<double, {}, {}> {}_real = (S*{}).real();\n".format(f.numRows, f.numCols, f.id, f.id) )
							#code += "\tEigen::Matrix<double, {}, {}> {}_real = (S*{}).real();\n".format(f.numRows, f.numCols, f.id, f.id)
							code += "\t\t{}_real_staggered[i] = (S*{}).real();\n".format(f.id, f.id)
						else:
							#file.write( "\tEigen::Matrix<double, {}, {}> {}_real = (S*{}*SInv).real();\n".format(f.numRows, f.numCols, f.id, f.id) )
							#code += "\tEigen::Matrix<double, {}, {}> {}_real = (S*{}*SInv).real();\n".format(f.numRows, f.numCols, f.id, f.id)
							code += "\t\t{}_real_staggered[i] = (S*{}*SInv).real();\n".format(f.id, f.id)
						code += "\t}\n"

						# now assemble final real matrix
						code += "\tEigen::Matrix<double, {}, {}> {}_real;\n".format(f.numRows, f.numCols, f.id)
						for i in range(pni.num_coeffs()):
							code += "\t{}_real.row({}) = {}_real_staggered[{}].row({});\n".format(f.id, i, f.id, pni.getGridIndexFromCoeff(i), i)


						file.write( code )
					else:
						file.write( "\t// all components vanish\n" )



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
		#print( "term {} ------------".format(key) )

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

			location = GridLocation3D(StaggeredGridLocation(np.array([0,0,0]), pni.getOffset(coeff_index)), StaggeredGridLocation(), "check") 
			#location = pni.getLocation( np.array([0,0,0]), coeff_index )

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

					#if u.interpolated == True:
					#	print("interpolated unknown!")
					
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

					if np.max(np.abs(u.voxel)) > stencil_width:
						stencil_width = np.max(np.abs(u.voxel))

					file.write( "\tctx.coeff_A( {}, vi + V3i({},{},{}), {} ) += {};\n".format(coeff_index, u.voxel[0], u.voxel[1], u.voxel[2], u.coeff_index, expr_cpp) )
					#file.write("\n")
			else:
				expr = meh.apply_recursive(result, ReplaceCoefficientMatrixComponents(coefficient_matrix_register))

				# Here we check if the coefficient is reduced to zero. This happens alot.
				if expr.canEvaluate() and np.abs(expr.evaluate()) < 1.0e-8:
					# coefficient is vanishes
					continue

				# the expression consists of function calls or numerical values
				# these are very easily translated into cpp code
				# the to_cpp function parses the expression tree and replaces all meh expression nodes
				# with a string of cpp code
				#print(meh.tree_str(expr))
				info = EvalInfo()
				info.coefficient_matrix_register = coefficient_matrix_register
				info.fun_to_cpp = fun_to_cpp
				expr_cpp = to_cpp(expr, info)

				file.write( "\tctx.coeff_b( {} ) += {};\n".format(coeff_index, expr_cpp) )


	file.write( "}\n" )
	file.write("V3i {}_get_offset(int coeff)\n".format(stencil_name))
	file.write("{\n")
	file.write("\tswitch(coeff)\n")
	file.write("\t{\n")
	for i in range(pni.num_coeffs()):
		offset = pni.getOffset(i)
		file.write("\t\tcase {}:return V3i({}, {}, {});break;\n".format(i, offset[0], offset[1], offset[2]))
	file.write("\t\tdefault:throw std::runtime_error(\"unexpected coefficient index\");break;\n")
	file.write("\t};\n")
	file.write("}\n")
	file.write( "REGISTER_STENCIL({}, {}, {}, {})\n".format(stencil_name, order, pni.num_coeffs(), stencil_width) )

	file.close()




if __name__ == "__main__":
	#order = 5
	#staggered = True
	#filename = "c:/projects/epfl/epfl17/cpp/pnsolver/src/stencil.cpp"

	#test = PNInfo3D(4, True)
	#exit(1)



	terms_sopn = []
	terms_sopn.append(rte_terms.lspn.term0_projected_expr())
	terms_sopn.append(rte_terms.lspn.term1_projected_expr())
	terms_sopn.append(rte_terms.lspn.term2_projected_expr())
	terms_sopn.append(rte_terms.lspn.term3_projected_expr())
	terms_sopn.append(rte_terms.lspn.term4_projected_expr())
	terms_sopn.append(rte_terms.lspn.term5_projected_expr())
	terms_sopn.append(rte_terms.lspn.term6_projected_expr())
	terms_sopn = rte_terms.splitAddition( terms_sopn )

	terms_fopn = []
	terms_fopn.append(rte_terms.fopn.transport_term())
	terms_fopn.append(rte_terms.fopn.collision_term())
	terms_fopn.append(rte_terms.fopn.scattering_term())
	terms_fopn.append(rte_terms.fopn.source_term())
	
	
	terms_fopn = rte_terms.splitAddition( terms_fopn )

	#for t in terms_fopn:
	#	print(meh.latex(t) + "+")

	#for term in terms_fopn:
	#	print("\n-------------------\n")
	#	meh.print_expr(term)


	staggered_id = {True:"sg", False:"cg"}
	terms = {"fopn":terms_fopn, "sopn":terms_sopn}

	#rte_forms = ["fopn", "sopn"]
	#staggered = [True, False]
	#order = [0,1,2,3,4,5]
	#order = [0,1,2,3,4]

	rte_forms = ["fopn"]
	#order = [1]
	order = [1,2,3,4,5]
	staggered = [True]

	test = itertools.product(rte_forms, order, staggered)
	for c in test:
		rte_form_name = c[0]
		rte_form_terms = terms[rte_form_name]
		order = c[1]
		is_staggered = c[2]
		#print( "terms={} order={} staggered={}".format(rte_form_name, order, str(staggered_state)))

		path = "c:/projects/epfl/epfl17/cpp/pnsolver/src"
		stencil_name = "stencil_{}_p{}_{}".format(rte_form_name, order, staggered_id[is_staggered] )
		filename = "{}/{}.cpp".format(path, stencil_name)

		#pni = PNInfo2D(order, is_staggered)
		pni = PNInfo3D(order, is_staggered)


		generate_stencil_code( stencil_name, filename, rte_form_terms, pni )

		# clear stencil file
		#file = open(filename, "w")
		#file.close()









