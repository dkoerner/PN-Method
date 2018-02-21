import datetime
import numpy as np
import itertools
import rte_terms
import meh


class PNInfo3D(object):
	def __init__(self, order, staggered=True):
		self.order = order

		# the following dictionaries help mapping lm to shindex and vice versa ---
		self.index_to_lm = [] # this will map equation index to l,m indices
		self.lm_to_index = {} # this dict will map l,m indices to the equation index
		# iterate all sh coefficients for given truncation order
		for l in range(0, self.order+1):
			for m in range(-l, l+1):
				self.index_to_lm.append( (l, m) )
				sh_index = len(self.index_to_lm)-1
				self.lm_to_index[(l,m)] = sh_index
				#print("check {} {}".format(sh_index, self.getIndex(l, m)))
		self.numCoeffs = len(self.index_to_lm)		

		# info about placement of unknowns on staggerd grid locations...
		self.unknown_info = [ {} for i in range(self.numCoeffs)]

		# by default we place all unknowns at the cell centers
		for i in range(self.numCoeffs):
			self.place_unknown(i, (1, 1, 1))

		# and we use full voxel central differences
		self.stencil_half_steps = 2

		# currently we manually place unknowns
		# TODO: come up with an automated approach and some motivation
		# for its solution.
		if staggered == True:
			if self.order >= 0:
				self.place_unknown( 0, (1, 1, 1) )
			if self.order >= 1:
				self.place_unknown( 1, (1, 0, 1) )
				self.place_unknown( 2, (1, 1, 0) )
				self.place_unknown( 3, (0, 1, 1) )
			if self.order >= 2:
				self.place_unknown( 4, (0, 0, 1) )
				self.place_unknown( 5, (1, 0, 0) )
				self.place_unknown( 6, (1, 1, 1) )
				self.place_unknown( 7, (0, 1, 0) )
				self.place_unknown( 8, (1, 1, 1) )
			if self.order >= 3:
				self.place_unknown( 9, (1, 0, 1) )
				self.place_unknown( 10, (0, 0, 0) )
				self.place_unknown( 11, (1, 0, 1) )
				self.place_unknown( 12, (1, 1, 0) )
				self.place_unknown( 13, (0, 1, 1) )
				self.place_unknown( 14, (1, 1, 0) )
				self.place_unknown( 15, (0, 1, 1) )
			if self.order >= 4:
				self.place_unknown( 16, (0, 0, 1) )
				self.place_unknown( 17, (1, 0, 0) )
				self.place_unknown( 18, (0, 0, 1) )
				self.place_unknown( 19, (1, 0, 0) )
				self.place_unknown( 20, (1, 1, 1) )
				self.place_unknown( 21, (0, 1, 0) )
				self.place_unknown( 22, (1, 1, 1) )
				self.place_unknown( 23, (0, 1, 0) )
				self.place_unknown( 24, (1, 1, 1) )
			if self.order >= 5:
				self.place_unknown( 25, (1, 0, 1) )
				self.place_unknown( 26, (0, 0, 0) )
				self.place_unknown( 27, (1, 0, 1) )
				self.place_unknown( 28, (0, 0, 0) )
				self.place_unknown( 29, (1, 0, 1) )
				self.place_unknown( 30, (1, 1, 0) )
				self.place_unknown( 31, (0, 1, 1) )
				self.place_unknown( 32, (1, 1, 0) )
				self.place_unknown( 33, (0, 1, 1) )
				self.place_unknown( 34, (1, 1, 0) )
				self.place_unknown( 35, (0, 1, 1) )
			if self.order > 5:
				raise ValueError("CHECK!")
			self.stencil_half_steps = 1

	def is2D(self):
		return False

	def getNumCoeffs(self):
		return self.numCoeffs

	def getIndex( self, l, m ):
		if l > self.order:
			return None
		return l*(l+1)+m 

	def getLMIndex( self, index ):
		return self.index_to_lm[index]

	def getUnknownSymbol(self):
		return "L"

	def place_unknown( self, coeff_index, grid_id ):
		self.unknown_info[coeff_index]['grid_id'] = grid_id
		self.unknown_info[coeff_index]['offset'] = np.array( [grid_id[0], grid_id[1], grid_id[2]] , dtype=int)

	def getOffset(self, coeff_index):
		return self.unknown_info[coeff_index]['offset']

	def getStencilHalfSteps(self):
		return self.stencil_half_steps



class PNInfo2D(object):
	def __init__(self, order, staggered=True):
		self.pni3d = PNInfo3D(order, staggered)

		self.index3d_to_index2d = [None for i in range(self.pni3d.getNumCoeffs())]
		self.index2d_to_index3d = []
		for l in range(0, order+1):
			for m in range(-l, l+1):
				index3d = self.pni3d.getIndex(l, m)
				# in 2d, we only need to solve for moments where l+m is even
				if (l+m) % 2 == 0:
					index2d = len(self.index2d_to_index3d)
					self.index2d_to_index3d.append(index3d)
				else:
					index2d = None
				self.index3d_to_index2d[index3d] = index2d
		self.numCoeffs = len(self.index2d_to_index3d)

	def getUnknownSymbol(self):
		return "L"

	def is2D(self):
		return True

	def getNumCoeffs(self):
		return self.numCoeffs

	def getIndex( self, l, m ):
		#print( "getIndex l={} m{}".format(l,m) )
		index3d = self.pni3d.getIndex(l,m)
		if not index3d is None:
			return self.index3d_to_index2d[index3d]
		return None

	def getLMIndex( self, index2d ):
		index3d = self.index2d_to_index3d[index2d]
		if not index3d is None:
			return self.pni3d.getLMIndex(index3d)
		return None, None

	def getOffset(self, index2d):
		index3d = self.index2d_to_index3d[index2d]
		if not index3d is None:
			return self.pni3d.getOffset(index3d)
		return None

	def getStencilHalfSteps(self):
		return self.pni3d.getStencilHalfSteps()





class Coefficient(object):
	def __init__(self, expr):
		self.expr = expr
	def getExpr(self):
		return self.expr
	def merge(self, b):
		self.expr = meh.add(self.expr, b.expr)

class Unknown(object):
	def __init__(self, expr):
		self.expr = expr
		#self.id_str = id_str
	#def getId(self):
	#	#return "{}[{}]".format(self.getSymbol(), self.getIndex())
	#	return self.id_str
	def getExpr(self):
		return self.expr
	def getSymbol(self):
		return self.expr.getSymbol()
	#def getIndex(self):
	#	return self.index

class Derivation(object):
	def __init__(self, child, variable):
		self.variable = variable
		self.child = child
	def getVariable(self):
		return self.variable
	def getChild(self):
		return self.child
	#def getId(self):
	#	return "d{}({})".format(self.getVariable().getSymbol(), self.child.getId())
class Product(object):
	def __init__( self, a, b ):
		self.a = a
		self.b = b
	#def getId(self):
	#	return "({})*({})".format(self.a.getId(), self.b.getId())

class TermStructure(object):
	def __init__( self, expr, unknown_symbol ):
		self.unknown_symbol = unknown_symbol
		self.root = self.build_recursive(expr)
		#self.id = self.root.getId()
		#print(r.getId())

	def getRoot(self):
		return self.root

	def merge( self, other ):
		#if self.id != other.getId():
		#	raise ValueError("unable to merge terms of different structure")
		self.merge_recursive(self.root, other.root)

	def merge_recursive(self, a, b):
		if a.__class__ == Coefficient and b.__class__ == Coefficient:
			# merge coefficients
			#a.expr = meh.add(a.expr, b.expr)
			a.merge(b)
		elif a.__class__ == Product and b.__class__ == Product:
			# merge childs
			self.merge_recursive( a.a, b.a )
			self.merge_recursive( a.b, b.b )
		elif a.__class__ == Derivation and b.__class__ == Derivation:
			self.merge_recursive( a.child, b.child )
		elif a.__class__ == Unknown and b.__class__ == Unknown:
			#if a.getId() != b.getId():
			if a.getSymbol() != b.getSymbol():
				raise ValueError("TermStructure::merge_recursive unknowns do not match")
		else:
			raise ValueError("TermStructure::merge_recursive unable to merge")

	#def getId(self):
	#	return self.id

	def print_info(self, node=None, level=0):
		if node == None:
			self.print_info(self.root)
			return
		istr = meh.indent_string(level)
		if node.__class__ == Coefficient:
			print("{}Coefficient expr={}".format(istr, meh.latex(node.expr)))
		elif node.__class__ == Product:
			print("{}Product".format(istr))
			self.print_info(node.a, level+1)
			self.print_info(node.b, level+1)
		elif node.__class__ == Derivation:
			print("{}Derivation".format(istr))
			self.print_info(node.child, level+1)
		elif node.__class__ == Unknown:
			print("{}Unknown id={}".format(istr, node.getSymbol()))
		else:
			raise ValueError("TermStructure::print_info unknown type {}".format(node.__class__.__name__))


	def getExpr(self, node=None, level=0):
		if node == None:
			return self.getExpr(self.root)

		if node.__class__ == Coefficient:
			return node.expr
		elif node.__class__ == Product:
			expr_a = self.getExpr(node.a)
			expr_b = self.getExpr(node.b)
			return meh.mul(expr_a, expr_b)
		elif node.__class__ == Derivation:
			expr_child = self.getExpr(node.child)
			return meh.deriv( expr_child, node.getVariable(), True )
		elif node.__class__ == Unknown:
			return node.expr			
		else:
			raise ValueError("TermStructure::getExpr unknown type {}".format(node.__class__.__name__))


	def build_recursive( self, expr, level=0 ):
		if expr.__class__ == meh.Number:
			return Coefficient(expr)
		elif expr.__class__ == meh.Negate:
			# get structure object from nested expression
			tmp = self.build_recursive(expr.getExpr(), level+1)
			if tmp.__class__ == Coefficient:
				# if the nested expression is a simple coefficient
				# merge the negation into its expression
				return Coefficient( meh.neg(tmp.expr) )
			elif tmp.__class__ == Product and (tmp.a.__class__ == Coefficient or tmp.b.__class__ == Coefficient):
				# if the nested expression is a product involving a coefficient...
				if tmp.a.__class__ == Coefficient:
					coeff = tmp.a
				else:
					coeff = tmp.b
				# then merge the negate into this one and return the product
				coeff.expr = meh.neg(coeff.expr)
				return tmp
			else:
				# otherwise just multiply by -1 coefficient
				return Product( Coefficient(meh.Number(-1)), tmp )
		elif isinstance(expr, meh.Variable):
			return Coefficient(expr)
		elif isinstance(expr, meh.Function):
			if expr.getSymbol() == self.unknown_symbol:
				return Unknown(expr)
			else:
				return Coefficient( expr )
		elif expr.__class__ == meh.Derivation:
			tmp = self.build_recursive(expr.getExpr(), level+1)
			if tmp.__class__ == Coefficient:
				return Coefficient( meh.deriv(tmp.expr, expr.getVariable(), is_partial=True)  )
			else:
				# assuming unknownvector or a DerivationOperator aplied to an unknown vector
				return Derivation(tmp, expr.getVariable() )
		elif expr.__class__ == meh.Multiplication:
			# we merge all coefficients into a single coefficient expression
			# and expect only a single non-coefficient element to multiply with
			numOperands = expr.numOperands()
			merged_coefficients = None
			other = None
			for i in range(numOperands):
				tmp = self.build_recursive(expr.getOperand(i), level+1)
				# if the expression is a coefficient
				if tmp.__class__ == Coefficient:
					# we check if we already came across a coefficient matrix
					if merged_coefficients == None:
						# no, then we use this coefficient matrix into which
						# we may merge other coefficients
						merged_coefficients = tmp
					else:
						# yes, there is a coefficient matrix, so merge other into it
						merged_coefficients.expr = meh.mul( merged_coefficients.expr, tmp.expr )
				elif tmp.__class__ != Coefficient and other is None:
					other = tmp
				elif tmp.__class__ != Coefficient:
					raise ValueError("expected only one operand which is not a coefficientmatrix")

			result = None
			if other is None and not merged_coefficients is None:
				# we have only a single coefficient
				return merged_coefficients
			elif merged_coefficients is None and not other is None:
				# we have only a single (non-coefficient) term
				return other
			else:
				# we have a coefficient term and a non-coefficient
				# which are multiplied
				return Product(merged_coefficients, other)
		elif expr.__class__ == meh.Power:
			tmp_base = self.build_recursive(expr.getBase(), level+1)
			tmp_exp = self.build_recursive(expr.getExponent(), level+1)
			if tmp_base.__class__ == Coefficient and tmp_exp.__class__ == Coefficient:
				return Coefficient( meh.pow(tmp_base.expr, tmp_exp.expr) )
			else:
				raise ValueError("unknowns or derivatives in potentials is power expressions are not supported .")
		elif expr.__class__ == meh.Quotient:
			tmp_num = self.build_recursive(expr.getNumerator(), level+1)
			tmp_denom = self.build_recursive(expr.getDenominator(), level+1)
			if tmp_num.__class__ == Coefficient and tmp_denom.__class__ == Coefficient:
				return Coefficient( meh.frac(tmp_num.expr, tmp_denom.expr) )
			else:
				raise ValueError("unknowns or derivatives in quotient expressions are not supported .")
		elif expr.__class__ == meh.Sqrt:
			tmp = self.build_recursive(expr.getExpr(), level+1)
			if tmp.__class__ == Coefficient:
				return Coefficient( meh.sqrt(tmp.expr) )
			else:
				raise ValueError("unknowns or derivatives in sqrt expressions are not supported .")
		else:
			raise ValueError("TermStructure::build_recursive unable to handle expression of type {}".format(expr.__class__.__name__))

		#meh.print_expr(expr)




def getUnknownIndex( pni, args ):
	# args[0] -> l
	l_expr = args[0]
	# args[1] -> m
	m_expr = args[1]

	if not (l_expr.canEvaluate() and m_expr.canEvaluate()):
		raise ValueError( "getUnknownIndex: unable to evaluate expressions for l or m" )

	return pni.getIndex(l_expr.evaluate(), m_expr.evaluate() )

'''
# getUnknownId is a function used by the factorization. 
# It returns an id for the args of a given unknown
# this is because the unknown actually is changed from a continuous, to a discretized
# version, which changes it parameters. We use factorization for both versions, which
# means that the continuous version should produce the same id for general x, while
# the discretized version should create different ids for different voxel coordinates.
def getUnknownId( pni, unknown_expr ):
	args = unknown_expr.getArguments()
	# the linear index of the given unknown within the
	# vector of SH coefficients. This is used during discretization
	# of the expressions which represent the expanded RTE.
	index = getUnknownIndex(pni, args)
	position_str = ""

	# if we are dealing with a discretized version, we take the voxel
	# position into account to discriminate unknowns
	if unknown_expr.__class__ == DiscretizedFunction:
		# args[2] -> i
		i_expr = args[2]
		# args[3] -> j
		j_expr = args[3]
		# args[4] -> k
		k_expr = args[4]

		if not (i_expr.canEvaluate() and j_expr.canEvaluate() and k_expr.canEvaluate()):
			raise ValueError( "getUnknownId: unable to evaluate expressions for i, j or k" )

		i = i_expr.evaluate()
		j = j_expr.evaluate()
		k = k_expr.evaluate()

		position_str = "(i={}, j={}, k={})".format(i, j, k)

	return pni.getUnknownSymbol() + "[{}]{}".format(index, position_str)
'''


def getTermStructureId( node ):
	if node.__class__ == Coefficient:
		return "C"
	elif node.__class__ == Unknown:
		#raise ValueError("where is pni defined !?!?!?")
		#print(pni)
		index = getUnknownIndex(pni, node.getExpr().getArguments())
		return node.getSymbol() + "[{}]".format(index)
	elif node.__class__ == Derivation:
		return "d{}({})".format(node.getVariable().getSymbol(), getTermStructureId(node.child))
	elif node.__class__ == Product:
		return "({})*({})".format(getTermStructureId(node.a), getTermStructureId(node.b))
	else:
		raise ValueError("getTermStructureId: unknown term structure node type")



def factorize(equ, pni):
	if equ.__class__ != meh.Addition:
		raise ValueError("factorize: can only work with additions")

	terms = equ.getOperands()
	term_structures = {}
	for term in terms:
		term_structure = TermStructure(term, pni.getUnknownSymbol())
		key = getTermStructureId(term_structure.getRoot())
		if not key in term_structures:
			term_structures[key] = term_structure
		else:
			# merge new term structure into existing one
			term_structures[key].merge(term_structure)

	factorized_terms = []
	for key, value in term_structures.items():
		factorized_terms.append(value.getExpr())

	return meh.Addition(factorized_terms)


class DiscretizedFunction(meh.Function):
	def __init__(self, original_function, voxel):
		# this class represents a function, which depended on continuous spatial
		# variable x and has been discretized into staggerd grid location.
		# parameters:
		# original_function: the original function (allows access symbol and arguments later)
		# voxel: the voxel at which this function will have to be evaluated according to x
		args = []
		lpos = []
		numArgs = original_function.numArguments()
		for i in range(numArgs):
			org_arg = original_function.getArgument(i)
			if isinstance(org_arg, meh.Variable) and org_arg.getSymbol() == "\\vec{x}":
				# we replace the positional argument with its discrete version
				# these are only placeholders and will be replaced by the voxel
				#args.append(meh.var("i"))
				#args.append(meh.var("j"))
				#args.append(meh.var("k"))
				args.append(meh.num(voxel[0]))
				args.append(meh.num(voxel[1]))
				args.append(meh.num(voxel[2]))

				# we have the problem that the function names are latex code,
				# which means that for function names containing subscripts
				# (such as \sigma_t), we cant easily add subscript arguments,
				# as this would break latex (\sigma_t_a is not allowed)
				# as a hack, we make ijk superscripts when there is a subscript
				# present
				# TODO:find a nicer solution
				if "_" not in original_function.getSymbol():
					lpos.append(-1)
					lpos.append(-1)
					lpos.append(-1)
				else:
					lpos.append(1)
					lpos.append(1)
					lpos.append(1)
			else:
				args.append(org_arg)
				lpos.append(original_function.getLatexArgumentPosition(i))

		super().__init__( original_function.getSymbol(), args)
		for i in range(len(args)):
			self.setLatexArgumentPosition(i, lpos[i])

		self.voxel = voxel
		self.original_function = original_function
	def deep_copy(self):
		return DiscretizedFunction( self.original_function, self.voxel )
	#def toHierarchy( self, indent = 0 ):
	#	return self.indent_string(indent) + "DiscretizedFunction {}\n".format(self.getSymbol())
	#def toLatex(self):
	#	l = self.get_l().toLatex()
	#	m = self.get_m().toLatex()
	#	i = "i"
	#	j = "j"
	#	k = "k"
	#	return "{}^{{{}, {}}}_{{ {}, {}, {} }}".format(self.getSymbol(), l, m, i, j, k )



class EmptyObject(object):
	pass


def apply_spatial_discretization_recursive( expr, info, level = 0 ):
	istr = meh.indent_string(level)
	if isinstance(expr, meh.Function) and expr.depends_on(meh.var('\\vec{x}')):
		#print("!!!!")
		#return expr
		# TODO: generalize this to any function which depends on x

		# we want to express the function at info.location, but due to
		# staggering, we might have it live at a different location.
		# Depending on the relation between location of evaluation and location
		# of where the function is defined, we might have to do some interpolation
		# between different values of that function at different locations.
		# This is what the following code is about.

		# eval_offset is the staggered grid location at which to evaluate the function
		eval_offset = info.location.getOffset()
		# function_offset is the staggered grid location at which the function is defined
		if expr.getSymbol() == info.pni.getUnknownSymbol():
			coeff_index = getUnknownIndex(info.pni, expr.getArguments())
			if coeff_index == None:
				# the coefficient vanishes (this happens in 2d case)
				return meh.num(0)				
			function_offset = info.pni.getOffset(coeff_index)
		else:
			# all other functions are always defined at cell centers
			function_offset = np.array( [1, 1, 1] , dtype=int)

		#print("eval_offset=")
		#print(eval_offset)
		#print("function_offset=")
		#print(function_offset)

		# now find in which dimension the staggered grid locations match
		same_axis = np.array([0 for i in range(3)])
		for i in range(3):
			if eval_offset[i] == function_offset[i]:
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
			u = DiscretizedFunction(expr, info.location.getVoxel())
			#u.interpolated = False
			# so we just return the expression as is
			return u
		else:
			#print(eval_offset)
			#print(coeff_index)
			#print(function_offset)
			#raise ValueError("interpolation not expected")
			# now we interpolate from the neighbous along the defined local axes
			offset_combinations = itertools.product(*[[-1, 1] for d in range(numAxes)])
			num_offset_combinations = 2**numAxes
			weight_expr = meh.frac(meh.num(1), meh.num(num_offset_combinations))
			terms = []
			for o in offset_combinations:
				offset = np.array([0, 0, 0])
				for i in range(numAxes):
					offset += o[i]*local_axes[i]
				u = DiscretizedFunction(expr, info.location.getShiftedLocation(offset).getVoxel())
				#u.interpolated = True
				terms.append( meh.mul(weight_expr, u) )
			# return linear interpolation
			return meh.Addition( terms )
	elif expr.__class__ == meh.Derivation:
		if expr.getVariable().getSymbol() == "x":
			dimension = 0
		elif expr.getVariable().getSymbol() == "y":
			dimension = 1
		elif expr.getVariable().getSymbol() == "z":
			dimension = 2
		else:
			raise ValueError("apply_spatial_discretization_recursive::derivation: unable to identify derivation variable")

		if info.pni.is2D() and dimension == 2:
			# we have a derivative in z, although we only work in 2d domain
			return meh.num(0)

		# stepsize determines the stepsize of the stencil in number of half-voxels
		stepsize = info.pni.getStencilHalfSteps()
		step = np.zeros(3, dtype=int)
		step[dimension] = stepsize

		# 
		location = info.location

		#central_difference_weight = meh.var("h_inv[{}]".format(dimension))
		central_difference_weight = meh.var("h^{{-1}}_{}".format("xyz"[dimension]))

		nested_expr = expr.getExpr()

		info.location = location.getShiftedLocation(-step)
		#info.vars["\\vec{x}"] = info.location
		a = apply_spatial_discretization_recursive( meh.mul(meh.neg(central_difference_weight), nested_expr.deep_copy()), info, level+1)

		info.location = location.getShiftedLocation(step)
		#info.vars["\\vec{x}"] = info.location
		b = apply_spatial_discretization_recursive( meh.mul(central_difference_weight, nested_expr.deep_copy()), info, level+1)

		info.location = location
		#info.vars["\\vec{x}"] = info.location

		result = meh.add(a, b)
		if info.debug == True:
			print("{}result={}".format(meh.indent_string(level), str(result)))
		return result
	else:
		# recurse
		for i in range(expr.numChildren()):
			child = apply_spatial_discretization_recursive(expr.getChildren(i), info, level+1)
			expr.setChildren(i, child)
		return expr



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
	# TODO: simplify location buisness
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



def generate_cpp_recursive(expr, info=None, level=0):
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
		return str( "-{}".format(generate_cpp_recursive(expr.getExpr(), info, level+1) ))
	elif expr.__class__ == meh.Variable:
		# this variable could not be resolved at stencil generation time
		# therefore we directly pass it along to the cpp code and
		# expect the user/client code to deal with it
		if expr.getSymbol() == "h^{-1}_x":
			return "h_inv[0]"
		elif expr.getSymbol() == "h^{-1}_y":
			return "h_inv[1]"
		elif expr.getSymbol() == "h^{-1}_z":
			return "h_inv[2]"

		return expr.getSymbol()
	elif isinstance(expr, meh.Function):
		numArgs = expr.numArguments()
		# turn all arguments into cpp
		arg_str = ""
		for i in range(numArgs):
			arg_str += generate_cpp_recursive(expr.getArgument(i), info, level+1)
			if i < numArgs-1:
				arg_str += ", "
		symbol = expr.getSymbol()
		if symbol == "\\sigma_t":
			return "ctx.evalExtinction({})[color_channel]".format(arg_str)
		elif symbol == "\\sigma_a":
			return "ctx.evalAbsorption({})[color_channel]".format(arg_str)
		elif symbol == "\\sigma_s":
			return "ctx.evalScattering({})[color_channel]".format(arg_str)
		elif symbol == "Q":
			return "ctx.evalEmission({})[color_channel]".format(arg_str)
		elif symbol == "f":
			return "ctx.evalPhase({})[color_channel]".format(arg_str)
		else:
			raise ValueError("unable to convert function {} to c++ code.".format(symbol))

		#return info.fun_to_cpp(expr, info, level)
		return "0"
	#elif expr.__class__ == GridLocation3D:
	#	loc = expr
	#	#return "TODO"
	#	#return "domain.voxelToWorld(vd+V3d({}, {}, {}))".format(expr.getVoxel()[0]+expr.getOffset()[0]*0.5, expr.getVoxel()[1]+expr.getOffset()[1]*0.5, expr.getVoxel()[2]+expr.getOffset()[2]*0.5 )
	#
	#	#return "domain.voxelToWorld(vd+V3d({}, {}, {}))".format(expr.getVoxel()[0]+expr.getOffset()[0]*0.5, expr.getVoxel()[1]+expr.getOffset()[1]*0.5, expr.getVoxel()[2]+expr.getOffset()[2]*0.5 )
	#
	#	#'''
	#	if loc.shift.isZero():
	#		grid_index = info.pni.getGridIndex(loc.getOffset())
	#		start_code = "+ctx.getVoxelSpaceOffsetFromGrid2({})".format(grid_index)
	#		if hasattr(loc.start, "to_cpp_str"):
	#			start_code = loc.start.to_cpp_str
	#		return "domain.voxelToWorld(vd{})".format(start_code)
	#	else:
	#		grid_index = info.pni.getGridIndex(loc.start.getOffset())
	#		shift = loc.shift.getVoxel()+loc.shift.getOffset()*0.5
	#		return "domain.voxelToWorld(vd+ctx.getVoxelSpaceOffsetFromGrid2({})+V3d({}, {}, {}))".format(grid_index, shift[0], shift[1], shift[2])		#'''	
	elif expr.__class__ == DiscretizedFunction:
		return "!!!!!!"
	elif expr.__class__ == meh.Addition:
		numOperands = expr.numOperands()
		cpp_ops = []
		max_str_le = 0
		for i in range(numOperands):
			cpp_op = generate_cpp_recursive(expr.getOperand(i), info, level+1)
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
			op_cpp = generate_cpp_recursive(op, info, level+1)
			result += op_cpp
			if i<numOperands-1:
				result+= "*"
		result += ")"
		return result
	elif expr.__class__ == meh.Power:
		base_cpp = generate_cpp_recursive(expr.getBase(), info)
		exponent_cpp = generate_cpp_recursive(expr.getExponent(), info)
		result = "std::pow({}, {})".format(base_cpp, exponent_cpp)
		return result
	elif expr.__class__ == meh.Quotient:
		num_cpp = generate_cpp_recursive(expr.getNumerator(), info)
		denom_cpp = generate_cpp_recursive(expr.getDenominator(), info)
		result = "({}/{})".format(num_cpp, denom_cpp)
		return result

	else:
		#pass
		raise ValueError("unable to handle expression of type {}".format(expr.__class__.__name__))
		#print("test {}".format(expr.__class__.__name__))



def generate_stencil_code( pne, pni, stencil_name ):
	numCoeffs = pni.getNumCoeffs()

	# the width of the stencil (the number of neighbouring voxels the stencil touches)
	# this values is determined later during construction the code for building A
	stencil_width = 0

	pne_discretized = []

	# discretize derivatives by applying spatial stencil ---
	for coeff_index in range(numCoeffs):

		#if coeff_index != 3:
		#	continue

		l,m = pni.getLMIndex(coeff_index)



		# get the pn-equation associated with the current index ---
		e = pne[coeff_index]

		#if coeff_index == 1:
		#	print("\n---------------------")
		#	print("l={} m={}".format(l,m))
		#	meh.print_expr(e)


		# factorize pn-equation if multiple terms are involved ---
		if e.__class__ == meh.Addition:
			e = factorize(e, pni)
			e = meh.apply_recursive(e, meh.CleanupSigns())

		# collapse constants
		e = meh.apply_recursive(e, meh.FoldConstants2())
		e = meh.apply_recursive(e, meh.CleanupSigns())



		#'''
		# discretize spatial variable into staggered grid ---
		info = EmptyObject()
		info.debug = False
		info.pni = pni
		# here we set the location at which to evaluate the given expression. During discretization, this
		# location may be modified (due to derivations etc.)
		info.location = GridLocation3D(StaggeredGridLocation(np.array([0,0,0]), pni.getOffset(coeff_index)), StaggeredGridLocation(), "check") 
		e = apply_spatial_discretization_recursive(e, info)
		# The discretization has replaced the unknown with a sum which we now expand
		e = meh.apply_recursive(e, meh.DistributiveLaw())
		e = meh.apply_recursive(e, meh.CleanupSigns())
		e = meh.apply_recursive(e, meh.DistributiveLaw())
		e = meh.apply_recursive(e, meh.CleanupSigns())
		e = meh.apply_recursive(e, meh.DistributiveLaw())
		e = meh.apply_recursive(e, meh.CleanupSigns())
		e = meh.apply_recursive(e, meh.DistributiveLaw())
		e = meh.apply_recursive(e, meh.CleanupSigns())

		e = meh.apply_recursive(e, meh.FoldConstants())

		#meh.print_expr(e)

		# store result ---
		pne_discretized.append(e)
		#pne[coeff_index] = e
		#'''

	# source file
	cpp = ""

	cpp += "// This file was generated by {} on {}\n\n".format(__file__, datetime.datetime.now().strftime("%A, %d. %B %Y %I:%M%p"))
	cpp += "#include <PNSystem.h>\n\n"

	#file.write("// truncation order is directly linked to the generated stencil\n")
	#file.write("int PNSystem::g_order = {};\n\n".format(order))
	

	arg_sys = "PNSystem::Stencil::Context& ctx"
	#arg_voxel = "const V2i& voxel"
	#prototype = "void {}({},\n\t\t\t\t\t{})\n".format(stencil_name, arg_sys, arg_voxel)
	prototype = "void {}({})\n".format(stencil_name, arg_sys)
	cpp += prototype
	cpp += "{\n" 
	cpp += "\tV3i vi = ctx.getVoxelCoord();\n"
	cpp += "\tV3d vd = vi.cast<double>();\n"
	cpp += "\tconst Domain& domain = ctx.getDomain();\n"
	cpp += "\tconst PNVolume& problem = ctx.getProblem();\n"
	cpp += "\tV3d h_inv( 1.0/({}*domain.getVoxelSize()[0]), 1.0/({}*domain.getVoxelSize()[1]), 1.0/({}*domain.getVoxelSize()[2]) );\n".format(pni.getStencilHalfSteps(), pni.getStencilHalfSteps(), pni.getStencilHalfSteps())
	cpp += "\tint color_channel = 0;\n"
	cpp += "\n"

	# now generate cpp code ==================================

	# here we iterate over all rows of the coefficient matrix within a single voxel
	for coeff_index in range(numCoeffs):
		e = pne_discretized[coeff_index]
		l,m = pni.getLMIndex(coeff_index)

		debug = False
		#if l==1 and m == -1:
		#	debug = True

		if debug == True:
			print("TEST --------------")
			meh.print_expr(e)


		#if coeff_index == 1:
		#	print("\n---------------------")
		#	print("l={} m={}".format(l,m))
		#	meh.print_expr(e)
		cpp += "\t// row={} l={} m={} --------------------------\n".format(coeff_index, l, m)

		unknown_coeffs = {}
		rhs_coeffs = None

		# iterate all individual terms of the current pn-equation ---
		ops = None
		if e.__class__ == meh.Addition:
			ops = e.getOperands()
		else:
			ops = [e]
		for op in ops:
			# analyse the term structure of the current term
			term_structure = TermStructure(op, pni.getUnknownSymbol())
			root = term_structure.getRoot()

			unknown = None
			coeff = None

			# if the term is just a single coefficient...
			if root.__class__ == Coefficient:
				# we know it belongs to the RHS as it does not contain an unknown
				coeff = root

				# if no rhs_coeff has been found yet...
				if not rhs_coeffs:
					# we initialize the rhs coefficient with this one
					rhs_coeffs = coeff
				else:
					# we merge the existing rhs coefficient with this one (basically adding them up)
					rhs_coeffs.merge(coeff)
			elif root.__class__ == Product:
				# the term is a product...we expect an unknown as part of the product
				if root.a.__class__ == Unknown:
					unknown = root.a
					coeff = root.b
				elif root.b.__class__ == Unknown:
					unknown = root.b
					coeff = root.a
				else:
					raise ValueError("unexpected")

				# we evaluate the voxel offset from the unknown. Note that at this point the information about
				# staggered grid location is not necessary anymore as that is defined by the unknown index
				args = unknown.getExpr().getArguments()
				index = getUnknownIndex(pni, args)
				# args[2] -> i
				i_expr = args[2]
				# args[3] -> j
				j_expr = args[3]
				# args[4] -> k
				k_expr = args[4]

				if not (i_expr.canEvaluate() and j_expr.canEvaluate() and k_expr.canEvaluate()):
					raise ValueError( "generate_stencil_code: unable to evaluate expressions for i, j or k" )

				i = i_expr.evaluate()
				j = j_expr.evaluate()
				k = k_expr.evaluate()

				vox = np.array([i, j, k])

				# we keep track of the maximum voxel offset in the stencil_width variable
				# this is later used to drive the number of layers of boundary voxels...
				if np.max(np.abs(vox)) > stencil_width:
					stencil_width = np.max(np.abs(vox))

				# at this point the unknown is defined by its index and its voxel offset
				unknown_id = (index, i, j, k)

				# the unknown_coeffs dictionary holds the coefficient for all the different unknowns
				# these coefficients are what gets evaluated into the coefficient matrix A or rhs b
				# since the same unknown may appear multiple times in the equation, we merge all of its
				# coefficients together. This is basically a factorization.
				if not unknown_id in unknown_coeffs:
					unknown_coeffs[unknown_id] = coeff
				else:
					# merge new term structure into existing one
					unknown_coeffs[unknown_id].merge(coeff)
			else:
				print(root.__class__.__name__)
				raise ValueError("unexpected")

		'''
		# now we have factorized all the unknowns. That is, for each unique unknown, we have a coefficient
		# expression.
		# what we now want to do is, to further factorize these unknown expressions. This is done so that
		# we may detect terms which completely vanishes, because the coefficient is a sum of terms which
		# cancel out each other. This happens alot and reduces the code significantly.
		if coeff_index != 1:
			continue

		for unknown_id, coeffs in unknown_coeffs.items():
			coeff_index_j = unknown_id[0]
			i = unknown_id[1]
			j = unknown_id[2]
			k = unknown_id[3]

			if coeff_index_j != 1:
				continue

			coeff_expr = coeffs.getExpr()
			#meh.print_expr(coeff_expr)
		#print("CHECK!")
		'''


		# now we generate the code for all the remaining unknowns and their coefficients
		for unknown_id, coeffs in unknown_coeffs.items():
			coeff_index_j = unknown_id[0]
			i = unknown_id[1]
			j = unknown_id[2]
			k = unknown_id[3]

			#if coeff_index_j != 1:
			#	continue


			coeff_expr = coeffs.getExpr()
			if coeff_expr.__class__ == meh.Addition:
				cpp += "\t{\n"
				#cpp += "\t\tstd::complex<double> c(0.0, 0.0);\n"
				cpp += "\t\tdouble c = 0.0;\n"
				ops = coeff_expr.getOperands()
				for op in ops:
					code = generate_cpp_recursive(op);
					cpp += "\t\tc+={};\n".format(code)
				#cpp += "\t\tctx.coeff_A( {}, vi+V3i({}, {}, {}), {} ) += c.real();\n".format(coeff_index, i, j, k, coeff_index_j)
				cpp += "\t\tctx.coeff_A( {}, vi+V3i({}, {}, {}), {} ) += c;\n".format(coeff_index, i, j, k, coeff_index_j)
				cpp += "\t}\n"
			else:
				code = generate_cpp_recursive(coeff_expr);
				#cpp += "\t\tc+={};\n".format(code)
				cpp += "\tctx.coeff_A( {}, vi+V3i({}, {}, {}), {} ) += {};\n".format(coeff_index, i, j, k, coeff_index_j, code)

		#'''
		# RHS
		if not rhs_coeffs is None:
			cpp += "\t{\n"
			#cpp += "\t\tstd::complex<double> c(0.0, 0.0);\n"
			cpp += "\t\tdouble c = 0.0;\n"
			if rhs_coeffs.getExpr().__class__ == meh.Addition:
				ops = rhs_coeffs.getExpr().getOperands()
				for op in ops:
					code = generate_cpp_recursive(op);
					cpp += "\t\tc+={};\n".format(code)
			else:
				code = generate_cpp_recursive(rhs_coeffs.getExpr());
				cpp += "\t\tc+={};\n".format(code)
			#cpp += "\t\tctx.coeff_b( {} ) += c.real();\n".format(coeff_index)
			cpp += "\t\tctx.coeff_b( {} ) += c;\n".format(coeff_index)
			cpp += "\t}\n"
		#'''




		#print("// row={} ---------------------".format(coeff_index))
		#code = generate_cpp_recursive(e)
		#print(code)

	cpp += "}\n"
	cpp += "V3i {}_get_offset(int coeff)\n".format(stencil_name)
	cpp += "{\n"
	cpp += "\tswitch(coeff)\n"
	cpp += "\t{\n"
	for i in range(pni.getNumCoeffs()):
		offset = pni.getOffset(i)
		cpp += "\t\tcase {}:return V3i({}, {}, {});break;\n".format(i, offset[0], offset[1], offset[2])
	cpp += "\t\tdefault:throw std::runtime_error(\"unexpected coefficient index\");break;\n"
	cpp += "\t};\n"
	cpp += "}\n"
	cpp += "REGISTER_STENCIL({}, {}, {}, {})\n".format(stencil_name, order, pni.getNumCoeffs(), stencil_width)

	return cpp






if __name__ == "__main__":

	'''
	staggered = True
	order = 5
	#pni = PNInfo3D(order, staggered)
	pni = PNInfo2D(order, staggered)

	terms = []
	# each function returns a list of expression representing the PN-equations up to given order
	terms.append(rte_terms.fopn_real.transport_term(order))
	terms.append(rte_terms.fopn_real.collision_term(order))
	terms.append(rte_terms.fopn_real.scattering_term(order))
	terms.append(rte_terms.fopn_real.source_term(order))

	pn_equations = []
	# we currently have a set pf pn-equations per RTE term
	# here we collapse the different equations for each term
	# into a single one for each SH coefficient index
	for i in range(pni.getNumCoeffs()):
		# in case of 2d, the coefficient indices are different
		l,m = pni.getLMIndex(i)
		index3d = l*(l+1)+m 

		gg = []
		for term in terms:
			gg.append( term[index3d] )
		if len(gg) == 1:
			pn_equations.append( gg[0] )
		else:
			pn_equations.append( meh.Addition(gg) )

	#staggered_id = {True:"sg", False:"cg"}
	#stencil_name = "stencil2_{}_p{}_{}".format("fopn", order, staggered_id[staggered] )
	#cpp = generate_stencil_code( pn_equations, pni, stencil_name )
	
	#print(cpp)

	#path = "c:/projects/epfl/epfl17/cpp/pnsolver/src"
	#filename = "{}/{}.cpp".format(path, stencil_name)
	#file = open(filename, "w")
	#file.write(cpp)
	#file.close()
	'''

	#dimension = ["2d", "3d"]
	dimension = ["3d"]
	#order = [1]
	#order = [2]
	#order = [3,5]
	order = [1,2,3,4,5]
	#staggered = [False, True]
	staggered = [True, False]

	staggered_id = {True:"sg", False:"cg"}
	test = itertools.product(order, dimension, staggered)
	for c in test:
		order = c[0]
		dim = c[1]
		is_staggered = c[2]

		if dim == "2d":
			pni = PNInfo2D(order, is_staggered)
		elif dim == "3d":
			pni = PNInfo3D(order, is_staggered)
		else:
			raise ValueError("invalid dimension value")

		## test
		#for t in range(9):
		#	a,b = pni.getLMIndex(t)
		#	print("index={} l={} m={}".format(t, a, b))

		# generate the mathematical expressions representing the PN equations ---
		terms = []
		# each function returns a list of expression representing the PN-equations up to given order
		terms.append(rte_terms.fopn_real.transport_term(order))
		terms.append(rte_terms.fopn_real.collision_term(order))
		terms.append(rte_terms.fopn_real.scattering_term(order))
		terms.append(rte_terms.fopn_real.source_term(order))

		pn_equations = []
		# we currently have a set pf pn-equations per RTE term
		# here we collapse the different equations for each term
		# into a single one for each SH coefficient index
		for i in range(pni.getNumCoeffs()):
			# in case of 2d, the coefficient indices are different
			l,m = pni.getLMIndex(i)
			index3d = l*(l+1)+m 

			gg = []
			for term in terms:
				gg.append( term[index3d] )
			if len(gg) == 1:
				pn_equations.append( gg[0] )
			else:
				pn_equations.append( meh.Addition(gg) )

		# generate cpp stencil code ----
		stencil_name = "stencil_p{}_{}_{}".format(order, dim, staggered_id[is_staggered] )
		cpp = generate_stencil_code( pn_equations, pni, stencil_name )

		# write cpp code to file ----
		path = "c:/projects/epfl/epfl17/cpp/pnsolver/src"
		filename = "{}/{}.cpp".format(path, stencil_name)
		print("writing stencil code to {}".format(filename))
		file = open(filename, "w")
		file.write(cpp)
		file.close()
