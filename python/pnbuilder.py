
import numpy as np
import util
import time
import shtools
import scipy.io
import cas
import itertools


class GridLocation2D(object):
	def __init__(self, domain, voxel_i , voxel_j, offset):
		self.domain = domain

		(voxel_offset_i, offset_i) = divmod( offset[0], 2 )
		(voxel_offset_j, offset_j) = divmod( offset[1], 2 )

		self.voxel_i = voxel_i + voxel_offset_i
		self.voxel_j = voxel_j + voxel_offset_j
		self.voxel = np.array([voxel_i + voxel_offset_i,voxel_j + voxel_offset_j])
		self.offset = np.array([offset_i,offset_j])
		#self.pWS = self.domain.bound_min + np.multiply(np.array([self.voxel_i+self.offset[0]*0.5, self.voxel_j+self.offset[1]*0.5]), self.domain.voxelsize)
		#self.pWS = self.domain.bound_min + np.multiply(self.voxel+self.offset*0.5, self.domain.voxelsize)
		x = (self.voxel[0]+self.offset[0]*0.5)*self.domain.voxelsize[0]
		y = (self.voxel[1]+self.offset[1]*0.5)*self.domain.voxelsize[1]
		self.pWS = np.array([x, y])
	def getPWS(self):
		# get world space position
		return self.pWS
	def getOffset(self):
		return self.offset
	def getVoxel(self):
		return self.voxel
	def getShiftedLocation(self, offset):
		return GridLocation2D(self.domain, self.voxel_i, self.voxel_j, self.offset+offset)
	def __str__(self):
		return "voxel={} {} offset={} {}".format(self.voxel[0], self.voxel[1], self.offset[0], self.offset[1])





class EvalInfo(object):
	def __init__(self):
		self.term_vanishes = False
		self.debug = False

class Unknown(object):
	def __init__(self, l, m, voxel, weight):
		self.l = l
		self.m = m
		self.voxel = voxel
		self.weight = weight



class UnknownInfo(object):
	def __init__(self, unknowns):
		self.unknowns = unknowns
	def __str__(self):
		result = "unknowns:\n"
		for u in self.unknowns:
			result += "\t l={} m={} voxel={} {} weight={}\n".format(u.l, u.m, u.voxel[0], u.voxel[1], u.weight)
		return result
	def __mul__(self, other):
		for u in self.unknowns:
			u.weight *= other
		return self
	def __rmul__(self, lhs):
		return self * lhs
	def __add__(self, other):
		return UnknownInfo( self.unknowns + other.unknowns )
	def __radd__(self, lhs):
		return self + lhs
	def __sub__(self, other):
		return self+ (-1.0)*other
	def __neg__(self):
		return (-1.0)*self



def eval_term_recursive( expr, info, level=0 ):
	istr = cas.indent_string(level)

	# in case of a simple number
	if expr.__class__ == cas.Number:
		if info.debug == True:
			print("{}eval_term_recursive::Number".format(istr))
		result = expr.getValue()
		return result
	elif expr.__class__ == cas.Negate:
		if info.debug == True:
			print("{}eval_term_recursive::Negate".format(istr))

		result = -eval_term_recursive(expr.getOperand(0), info, level+1)
		return result
	elif expr.__class__ == cas.Quotient:
		if info.debug == True:
			print("{}eval_term_recursive::Quotient".format(istr))

		result = eval_term_recursive(expr.getNumerator(), info, level+1)/eval_term_recursive(expr.getDenominator(), info, level+1)
		return result
	elif expr.__class__ == cas.Addition:
		if info.debug == True:
			print("{}eval_term_recursive::Addition".format(istr))

		numOperands = expr.numOperands()
		result = 0.0
		for i in range(numOperands):
			result += eval_term_recursive(expr.getOperand(i), info, level+1)
		return result
	#elif expr.__class__ == cas.Variable:
	elif isinstance(expr, cas.Variable):
		if info.debug == True:
			print("{}eval_term_recursive::Variable".format(istr))

		result = None
		# TODO: if we come across the unknown, then return a stencil point
		if expr.__class__ == cas.ImaginaryUnit:
			result = complex(0.0, 1.0)
		elif expr.getSymbol() in info.vars:
			result = info.vars[expr.getSymbol()]
		else:
			raise ValueError("unable to resolve variable {}".format(expr.getSymbol()))

		if info.debug == True:
			print("{}result={}".format(istr, str(result)))

		return result
	# in case of a function
	elif isinstance(expr, cas.Function):
		if info.debug == True:
			print("{}eval_term_recursive::Function".format(istr))
		numArgs = expr.numArguments()

		# evaluate all arguments
		args = []
		for i in range(numArgs):
			args.append( eval_term_recursive(expr.getArgument(i), info, level+1) )
		if expr.getSymbol() == info.unknown_symbol:
			# currently, we assume that:
			# args[0] -> l;args[1] -> m;args[2] -> x
			l = args[0]
			m = args[1]

			if l < 0 or l > info.builder.N:
				info.term_vanishes = True
				return 0.0

			coeff_index = info.builder.shIndex(l, m)
			# check bounds
			if coeff_index == None or coeff_index >=info.builder.numCoeffs:
				info.term_vanishes = True
				return 0.0

			numDimensions = 2

			# check if the location, at which to evaluate the unknown,
			# matches the actual grid location of the unknown
			# this is true for first order equation with non-anisotropic media
			location_offset = info.location.getOffset()
			unknown_offset = info.builder.unknown_info[coeff_index]['offset']
			if (location_offset == unknown_offset).all():
				# unknown location and eval location are the same spot
				# no interpolation needed
				u = Unknown(l, m, info.location.getVoxel(), 1.0)
				result = UnknownInfo([u])
			elif location_offset[0] == unknown_offset[0]:
				# in the previous if-clause, we checked for location and unknown to be exactly equal
				# now if their offset matches only in one dimension, then we can simply interpolate
				# between the two neighbouring datapoints in that dimension

				# TODO: generalize to 3D

				u0 = Unknown(l, m, info.location.getShiftedLocation(np.array([0, 1])).getVoxel(), 0.5)
				u1 = Unknown(l, m, info.location.getShiftedLocation(np.array([0, -1])).getVoxel(), 0.5)
				return UnknownInfo([u0,u1])
			elif location_offset[1] == unknown_offset[1]:
				u0 = Unknown(l, m, info.location.getShiftedLocation(np.array([1, 0])).getVoxel(), 0.5)
				u1 = Unknown(l, m, info.location.getShiftedLocation(np.array([-1, 0])).getVoxel(), 0.5)
				result = UnknownInfo([u0,u1])
			#elif (location_offset[0]+1)%2 == unknown_offset[0] and (location_offset[1]+1)%2 == unknown_offset[1]:
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
					u = Unknown(l, m, info.location.getShiftedLocation(np.array(o)).getVoxel(), weight)
					unknowns.append(u)
				result = UnknownInfo(unknowns)

		elif expr.getSymbol() in info.functions:
			# evaluate function
			result = info.functions[expr.getSymbol()]( *args )
		elif not expr.body2 is None:
			result = expr.body2(*args)
		else:
			raise ValueError("function {} not defined for evaluation".format(expr.getSymbol()))

		if info.debug == True:
			print("{}result={}".format(cas.indent_string(level), str(result)))


		return result
	# derivation...
	elif expr.__class__ == cas.Derivation:
		if info.debug == True:
			print("{}eval_term_recursive::Derivation".format(cas.indent_string(level)))
		if expr.getVariable().getSymbol() == "x":
			# NB: grid has x-dimension along columns (therefore dim=1)
			dimension = 0
		elif expr.getVariable().getSymbol() == "y":
			dimension = 1
		elif expr.getVariable().getSymbol() == "z":
			dimension = 2
		else:
			raise ValueError("unable to identify derivation variable")

		if dimension >= info.voxelsize.shape[0]:
			# we have a derivative in z although we only work in 2d domain
			info.term_vanishes = True
			return 0.0

		# stepsize determines the stepsize of the stencil in number of half-voxels
		stepsize = info.builder.stencil_half_steps

		step = np.zeros(info.voxelsize.shape[0], dtype=int)
		step[dimension] = stepsize

		if info.debug == True:
		#if True:
			print("{}dimension={}".format(cas.indent_string(level), dimension))
			#print("step=")
			#print(step)
			print("{}info.location={}".format(cas.indent_string(level), info.location))



		location = info.location

		central_difference_weight = 1.0/(stepsize*info.voxelsize[dimension])

		nested_expr = expr.getExpr()
		#print("{} class={}".format(cas.indent_string(level+1), nested_expr.__class__))

		# idea behind this is to evaluate the child expression for the different positions
		# of the discretization stencils
		info.location = location.getShiftedLocation(-step)
		info.vars["\\vec{x}"] = info.location.getPWS()
		if info.debug == True:
			print("{}-step={} {}".format(cas.indent_string(level), info.location.voxel[0], info.location.voxel[1]))

		#info.prefix += "d" + expr.getVariable().getSymbol()
		a = eval_term_recursive(nested_expr, info, level+1)

		info.location = location.getShiftedLocation(step)
		info.vars["\\vec{x}"] = info.location.getPWS()
		if info.debug == True:
			print("{}+step={} {}".format(cas.indent_string(level), info.location.voxel[0], info.location.voxel[1]))
		b = eval_term_recursive(nested_expr, info, level+1)

		info.location = location
		info.vars["\\vec{x}"] = info.location.getPWS()
		return central_difference_weight*(b - a)
	elif expr.__class__ == cas.Multiplication:
		if info.debug == True:
			print("{}eval_term_recursive::Multiplication".format(cas.indent_string(level)))
		numOperands = expr.numOperands()
		result = 1.0
		for i in range(numOperands):
			result = result * eval_term_recursive(expr.getOperand(i), info, level+1)
		return result
	elif expr.__class__ == cas.Power:
		if info.debug == True:
			print("{}eval_term_recursive::Power".format(cas.indent_string(level)))
		return eval_term_recursive(expr.getBase(), info)**eval_term_recursive(expr.getExponent(), info, level+1)
	else:
		raise ValueError("unable to handle expression of type {}".format(expr.__class__.__name__))



class PNBuilder(object):
	'''
	'''


	def __init__(self, order, domain):
		'''The init function takes approximation order and domain (including discretization info)'''

		self.N = order
		self.domain = domain

		self.index_to_lm = [] # this will map equation index to l,m indices
		self.lm_to_index = {} # this dict will map l,m indices to the equation index
		# NB: we dont use l(l+1)+m because in 2d, we skipp odd (l+m) entries, which changes the sequence.
		# iterate all sh coefficients for given truncation order
		for l in range(0, self.N+1):
			for m in range(-l, l+1):
				# in 2d, we only need to solve for moments where l+m is even
				if (l+m) % 2 == 0:
					self.index_to_lm.append( (l, m) )
					self.lm_to_index[(l,m)] = len(self.index_to_lm)-1
		self.numCoeffs = len(self.index_to_lm)
		#print("numCoeffs={}".format(self.numCoeffs))

		self.build_S()
		#self.build_M()

		# X and Y are scalar fields which store the x and y components of the center for each voxel
		offset = (0.5,0.5)
		self.X = np.array([np.arange(self.domain.res_x) for y in range(self.domain.res_y)]).T*self.domain.h_x + offset[0]*self.domain.h_x
		self.Y = np.array([np.arange(self.domain.res_y) for x in range(self.domain.res_x)])*self.domain.h_y + offset[1]*self.domain.h_y


		# extract the coefficients to the SH coefficients Llm
		self.terms = []

		self.unknown_info = [ {} for i in range(self.numCoeffs)]

		self.stencil_half_steps = 1



	def shIndex(self, l, m):
		key = (l,m)
		#NB: we dont use l(l+1)+m because in 2d, we skipp odd (l+m) entries, which changes the sequence.
		if key in self.lm_to_index:
			return self.lm_to_index[key]
		return None

	def build_S(self):
		'''builds the S matrix, which converts from complex-valued to real valued coefficients'''
		# build S matrix ( we iterate over l, m to make sure that the order is correct)

		self.S = np.zeros((self.numCoeffs, self.numCoeffs),dtype=complex)
		count = 0
		for l in range(0, self.N+1):
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

	def place_unknown( self, coeff_index, grid_id ):
		self.unknown_info[coeff_index]['grid_id'] = grid_id
		self.unknown_info[coeff_index]['offset'] = np.array( [grid_id[0], grid_id[1]] , dtype=int)

	def set_stencil_half_steps(self, stencil_half_steps):
		self.stencil_half_steps = stencil_half_steps

	def add_terms(self, expr):
		if expr.__class__ == cas.Multiplication or expr.__class__ == cas.Negate or isinstance(expr, cas.Function):
			self.terms.append(expr)
		elif expr.__class__ == cas.Addition:
			numTerms = expr.numOperands()
			for i in range(numTerms):
				self.terms.append(expr.getOperand(i))
		else:
			ValueError("expected expression to be addition or multiplication")



	#def assemble_global_matrix(self):
	#	'''This function assembles the global matrix which expresses the PN equations for all voxels
	#	and all SH coefficients in one _big_ matrix.'''

	def get_global_index( self, voxel_i, voxel_j, coeff ):
		'''Returns the equation index for the given h coeffient at the given voxel.'''
		voxel = voxel_j*self.domain.res_x + voxel_i
		return voxel*self.numCoeffs + coeff


	def build_global( self, functions ):
		'''The run method takes all problem specific inputs, such as absorption-
		and scattering coefficient fields and source field and assembles a global
		matrix which is used to solve for all sh coefficients of all voxels.'''
		print("building global systen Ax=b...")

		numVoxels = self.domain.res_x*self.domain.res_y

		# coefficient matrix A and rhs b of our global problem
		A_complex = np.zeros( (numVoxels*self.numCoeffs, numVoxels*self.numCoeffs), dtype = complex )
		b_complex = np.zeros( (numVoxels*self.numCoeffs), dtype=complex )
		A_real = np.zeros( (numVoxels*self.numCoeffs, numVoxels*self.numCoeffs) )
		b_real = np.zeros( (numVoxels*self.numCoeffs) )

		# Now we build the global coefficient matrix A and rhs b by iterating over all elements (voxels)
		# and within each voxel, we iterate over all sh coefficients. Each row within the global system
		# is associated with a specific sh coefficient of a specific voxel.
		# Then, for each row, we evaluate all terms and accumulate the results into A and b.
		for voxel_x in range(self.domain.res_x):
			print("voxel_x={}".format(voxel_x))

			#if voxel_x > 0:
			#	continue

			#if voxel_x > 10:
			#	break

			for voxel_y in range(self.domain.res_y):

				#if voxel_y > 0:
				#	continue

				# temp
				#global_i_min = self.get_global_index(voxel_x, voxel_y, 0)
				#global_i_max = self.get_global_index(voxel_x, voxel_y, self.numCoeffs)
				#if 6389 < global_i_min or 6389 > global_i_max:
				#	continue


				# here we iterate the local matrix vector product between some matrix and u
				# M represents the coupling between coefficients at the same location
				# u contains the sh coefficients for the current voxel
				# the spatial derivative will cause dependencies on u's from other voxels
				# which we can express easily because we have a global system
				for local_i in range(self.numCoeffs):

					#if local_i > 0:
					#	continue

					# find the equation index within our global system
					global_i = self.get_global_index(voxel_x, voxel_y, local_i)

					#if global_i != 6389:
					#	continue

					# get (grid)location of the unknown which is associated with the current row
					location = GridLocation2D(self.domain, voxel_x, voxel_y, self.unknown_info[local_i]['offset'])
					pWS = np.array([self.X[voxel_x, voxel_y], self.Y[voxel_x, voxel_y]])
					pWS2 = location.getPWS()
					#print("voxel={} {} center={} {}".format(voxel_i, voxel_j, pWS[0], pWS[1]))
					#print("pWS={} {} check={} {}".format(pWS[0], pWS[1], pWS2[0], pWS2[1]))
					#continue


					l = self.index_to_lm[local_i][0]#-> use this for 3d: l,m = shtools.lmIndex(i)
					m = self.index_to_lm[local_i][1]

					# now we have a specific row in of our Ax=b system
					# what we now do is to simply iterate over all lhs terms, evaluate them and accumulate
					# coefficients into A
					# the tricky part is to take derivatives into account, which is done by using stencils
					# in a generalized way
					#for term in self.lhs_terms:
					term_index = 0
					#print("numTerms={}".format(len(self.terms)))
					for term in self.terms:

						

						#print("term={}".format(term.toLatex()))
						#print(cas.hierarchy(term))

						info = EvalInfo()
						info.unknown_symbol = "L"

						info.vars = {}
						info.vars["\\vec{x}"] = pWS
						info.vars["l'"] = l
						info.vars["m'"] = m
						info.coeff_equ = local_i
						info.prefix = ""

						info.functions = functions
						info.voxelsize = self.domain.voxelsize
						# location at which to evaluate the current term
						# this is driven by the unknown which is associated with the current row
						info.location = location
						info.builder = self
						#info.debug = True

						#if global_i == 6447 and term_index == 1:
						#	info.debug = True
						#	print("term debug {} ------------------------------".format(term_index))

						if info.debug == True:
							cas.print_expr(term)
							print("location={}".format(info.location))
						result = eval_term_recursive( term, info )

						if info.term_vanishes:
							if info.debug == True:
								print("term vanishes")
							# term vanishes
							# happens for example if it involves z-derivative in 2d
							pass
						elif result.__class__ == UnknownInfo:
							if info.debug == True:
								print("term_index={} unknown info".format(term_index))
								print("#unknowns={}".format(len(result.unknowns)))

							# this is a lhs term
							# weights are going into A
							for u in result.unknowns:
								if u.l < 0 or u.l > self.N or u.weight == 0.0:
									if info.debug == True:
										print("skipping unknown u.l={} u.weight={} {}".format(u.l, np.real(u.weight), np.imag(u.weight)))
									continue

								local_j = self.shIndex( u.l, u.m )
								#print("coefficient {} depends on coefficient {} with weight={}".format(local_i, local_j, u.weight) )

								# TODO: think about how to handle boundaries
								# TODO: probably we will have that weight=0.0 for out of bound locations
								if u.voxel[0] < 0 or u.voxel[0] >= self.domain.res_x or u.voxel[1] < 0 or u.voxel[1] >= self.domain.res_y:
									continue

								# voxel and shcoeff-index define the final column within the current row
								# of the global matrix A
								global_j = self.get_global_index(u.voxel[0], u.voxel[1], local_j)


								#if global_i == 6447 and global_j == 6448:
								#	print("voxel_i={} voxel_j={} term_index={} u.l={} u.m={}".format(voxel_x, voxel_y, term_index, u.l, u.m))
								#	print("u.weight={}".format(u.weight))

								A_complex[global_i, global_j] += u.weight
						else:
							if info.debug == True:
								print("b term")

							#print("no unknown!")
							#if np.isnan(result):
							#	print("term index: {}".format(ii))
							#	print("l={} m={}".format(l,m))
							#	raise ValueError()
							# this is a rhs term (because it has no unknowns and evaluated to a number)
							#print("result={} {}".format(np.real(result), np.imag(result)))
							#this goes straight into b
							b_complex[global_i] += result
						term_index += 1

				# now transform all blocks of the current block-row into real variables
				# TODO: this can be optimized by analysing which blocks are zero
				#'''
				for voxel_x2 in range(self.domain.res_x):
					for voxel_y2 in range(self.domain.res_y):
						block_i = self.get_global_index(voxel_x, voxel_y, 0)
						block_j = self.get_global_index(voxel_x2, voxel_y2, 0)
						M_complex = A_complex[block_i:block_i + self.numCoeffs, block_j:block_j + self.numCoeffs]

						A_real[block_i:block_i + self.numCoeffs, block_j:block_j + self.numCoeffs] = np.real(self.S.dot(M_complex.dot(self.S_inv)))
						b_real[block_i:block_i + self.numCoeffs] = np.real(self.S.dot(b_complex[block_i:block_i + self.numCoeffs]))
				#'''
				#A_real = np.real(A_complex)
				#b_real = np.real(b_complex)






		# todo: transform this into a system of real variables
		return (A_real,b_real)








def test_sigma_t( pWS ):
	#return 1.123
	return pWS[0]

if __name__ == "__main__":
	x = cas.tensor("\\vec{x}", rank=1, dimension=3)
	x.setComponent(0, cas.var("x"))
	x.setComponent(1, cas.var("y"))
	x.setComponent(2, cas.var("z"))
	x.collapsed = True

	L_coeffs = cas.SHCoefficient( "L", cas.var("l'"), cas.var("m'"), x )

	#expr = cas.num(10)
	#expr = cas.var("\\vec{x}")
	#expr = x
	#expr = cas.fun("\\sigma_t", x)
	#dx_sigma_t = cas.deriv( cas.fun("\\sigma_t", x), x.getComponent(0) )
	#expr = dx_sigma_t
	#expr = cas.mul(cas.num(2.0), dx_sigma_t)
	#expr = L_coeffs
	#expr = cas.mul( cas.num(2.0), L_coeffs )
	dx_L = cas.deriv(cas.deriv( L_coeffs, x.getComponent(1) ), x.getComponent(0))
	expr = cas.mul(cas.num(2.0), dx_L)
	#expr = cas.neg(10)




	coord = np.array([0,0], dtype=int)
	voxelsize = np.array([1.0, 1.0])
	pWS = np.zeros(2)
	#stencil = util.stencil2d(coord, pWS, voxelsize)

	unknown_symbol = "L"
	vars = {}
	vars["\\vec{x}"] = pWS
	vars["l'"] = 0
	vars["m'"] = 0
	funs = { "\\sigma_t":test_sigma_t }
	eval_term( expr, unknown_symbol, vars, funs, coord, voxelsize )