
import numpy as np
import util
import time
import shtools
import scipy.io
import cas



class EvalInfo(object):
	pass

class Unknown(object):
	pass


class UnknownInfo(object):
	def __init__(self, unknowns):
		self.unknowns = unknowns
	def __str__(self):
		result = "unknowns:\n"
		for u in self.unknowns:
			result += "\t l={} m={} coord={} {} weight={}\n".format(u.l, u.m, u.coord[0], u.coord[1], u.weight)
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

	# in case of a simple number
	if expr.__class__ == cas.Number:
		#print("{}eval_term_recursive::Number".format(cas.indent_string(level)))
		return expr.getValue()
	elif expr.__class__ == cas.Negate:
		#print("{}eval_term_recursive::Negate".format(cas.indent_string(level)))
		return -eval_term_recursive(expr.getOperand(0), info, level+1)
	elif expr.__class__ == cas.Quotient:
		#print("{}eval_term_recursive::Quotient".format(cas.indent_string(level)))
		return eval_term_recursive(expr.getNumerator(), info, level+1)/eval_term_recursive(expr.getDenominator(), info, level+1)
	elif expr.__class__ == cas.Addition:
		#print("{}eval_term_recursive::Addition".format(cas.indent_string(level)))
		numOperands = expr.numOperands()
		sum = 0.0
		for i in range(numOperands):
			sum += eval_term_recursive(expr.getOperand(i), info, level+1)
		return sum
	#elif expr.__class__ == cas.Variable:
	elif isinstance(expr, cas.Variable):
		#print("{}eval_term_recursive::Variable".format(cas.indent_string(level)))
		# TODO: if we come across the unknown, then return a stencil point
		if expr.__class__ == cas.ImaginaryUnit:
			return complex(0.0, 1.0)
		elif expr.getSymbol() in info.vars:
			return info.vars[expr.getSymbol()]
		else:
			raise ValueError("unable to resolve variable {}".format(expr.getSymbol()))
	# in case of a function
	elif isinstance(expr, cas.Function):
		#print("{}eval_term_recursive::Function".format(cas.indent_string(level)))
		numArgs = expr.numArguments()
		# evaluate all arguments
		args = []
		for i in range(numArgs):
			args.append( eval_term_recursive(expr.getArgument(i), info, level+1) )
		if expr.getSymbol() == info.unknown_symbol:
			u = Unknown()
			# currently, we assume that:
			# args[0] -> l;args[1] -> m;args[2] -> x
			u.l = args[0]
			u.m = args[1]
			u.x = args[2]
			u.coord = info.coord
			u.weight = 1.0
			return UnknownInfo([u])
		elif expr.getSymbol() in info.functions:
			# evaluate function
			return info.functions[expr.getSymbol()]( *args )
		elif not expr.body2 is None:
			return expr.body2(*args)
		raise ValueError("function {} not defined for evaluation".format(expr.getSymbol()))
	# derivation...
	elif expr.__class__ == cas.Derivation:
		#print("{}eval_term_recursive::Derivation".format(cas.indent_string(level)))
		step = np.zeros(info.voxelsize.shape[0], dtype=int)
		pWS = info.vars['\\vec{x}']
		coord = info.coord
		if expr.getVariable().getSymbol() == "x":
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

		step[dimension] = 1
		stepWS = step*info.voxelsize
		central_difference_weight = 1.0/(2.0*info.voxelsize[dimension])

		nested_expr = expr.getExpr()
		#print("{} class={}".format(cas.indent_string(level+1), nested_expr.__class__))

		# idea behind this is to evaluate the child expression for the different positions
		# of the discretization stencils
		info.vars['\\vec{x}'] = pWS - stepWS
		info.coord = coord - step
		a = eval_term_recursive(nested_expr, info, level+1)

		info.vars['\\vec{x}'] = pWS + stepWS
		info.coord = coord + step
		b = eval_term_recursive(nested_expr, info, level+1)

		info.vars['\\vec{x}'] = pWS
		return central_difference_weight*(b - a)
	elif expr.__class__ == cas.Multiplication:
		#print("{}eval_term_recursive::Multiplication".format(cas.indent_string(level)))
		numOperands = expr.numOperands()
		result = 1.0
		for i in range(numOperands):
			result = result * eval_term_recursive(expr.getOperand(i), info, level+1)
		return result
	elif expr.__class__ == cas.Power:
		#print("{}eval_term_recursive::Power".format(cas.indent_string(level)))
		return eval_term_recursive(expr.getBase(), info)**eval_term_recursive(expr.getExponent(), info, level+1)
	else:
		raise ValueError("unable to handle expression of type {}".format(expr.__class__.__name__))


def eval_term(expr, unknown_symbol, vars, funs, coord, voxelsize):
	info = EvalInfo()
	info.unknown_symbol = unknown_symbol
	info.vars = vars
	info.functions = funs
	info.voxelsize = voxelsize
	info.coord = coord
	info.term_vanishes = False
	result = eval_term_recursive(expr, info)

	if info.term_vanishes == True:
		return None

	return result
	



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


	def build_global( self, sigma_a = None, sigma_s = None, phase_shcoeffs = None, source_shcoeffs = None ):
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

			if voxel_x > 10:
				break

			for voxel_y in range(self.domain.res_y):

				# get location of voxel center
				pWS = np.array([self.X[voxel_x, voxel_y], self.Y[voxel_x, voxel_y]])
				#print("voxel={} {} center={} {}".format(voxel_i, voxel_j, pWS[0], pWS[1]))

				# here we iterate the local matrix vector product between some matrix and u
				# M represents the coupling between coefficients at the same location
				# u contains the sh coefficients for the current voxel
				# the spatial derivative will cause dependencies on u's from other voxels
				# which we can express easily because we have a global system
				for local_i in range(self.numCoeffs):

					# find the equation index within our global system
					global_i = self.get_global_index(voxel_x, voxel_y, local_i)


					l = self.index_to_lm[local_i][0]#-> use this for 3d: l,m = shtools.lmIndex(i)
					m = self.index_to_lm[local_i][1]

					# now we have a specific row in of our Ax=b system
					# what we now do is to simply iterate over all lhs terms, evaluate them and accumulate
					# coefficients into A
					# the tricky part is to take derivatives into account, which is done by using stencils
					# in a generalized way
					#for term in self.lhs_terms:
					ii = 0
					#print("numTerms={}".format(len(self.terms)))
					for term in self.terms:
						#print("term={}".format(term.toLatex()))
						#print(cas.hierarchy(term))
						vars = {}
						vars["\\vec{x}"] = pWS
						vars["l'"] = l
						vars["m'"] = m
						functions = {}
						functions["\\sigma_t"] = lambda pWS: sigma_a(pWS) + sigma_s(pWS)
						functions["\\sigma_a"] = sigma_a
						functions["\\sigma_s"] = sigma_s
						functions["f_p"] = phase_shcoeffs
						functions["q"] = source_shcoeffs
						coord = np.array([voxel_x, voxel_y], dtype=int)
						result = eval_term( term, "L", vars, functions, coord, self.domain.voxelsize )
						if result is None:
							# term vanishes
							# happens for example if it involves a derivative in z in 2d
							pass
						elif result.__class__ == UnknownInfo:
							# this is a lhs term
							# weights are going into A
							for u in result.unknowns:
								if u.l < 0 or u.l > self.N:
									continue

								local_j = self.shIndex( u.l, u.m )
								# check bounds
								if local_j == None or local_j >=self.numCoeffs:
									continue

								# TODO: think about how to handle boundaries
								if u.coord[0] < 0 or u.coord[0] >= self.domain.res_x or u.coord[1] < 0 or u.coord[1] >= self.domain.res_y:
									continue

								# position and shcoeff-index define the final column within the current row
								# of the global matrix A
								global_j = self.get_global_index(u.coord[0], u.coord[1], local_j)

								A_complex[global_i, global_j] += u.weight
						else:
							#print("no unknown!")
							#if np.isnan(result):
							#	print("term index: {}".format(ii))
							#	print("l={} m={}".format(l,m))
							#	raise ValueError()
							# this is a rhs term (because it has no unknowns and evaluated to a number)
							#this goes straight into b
							b_complex[global_i] += result
						ii += 1

				# now transform all blocks of the current block-row into real variables
				# TODO: this can be optimized by analysing which blocks are zero
				for voxel_x2 in range(self.domain.res_x):
					for voxel_y2 in range(self.domain.res_y):
						block_i = self.get_global_index(voxel_x, voxel_y, 0)
						block_j = self.get_global_index(voxel_x2, voxel_y2, 0)
						M_complex = A_complex[block_i:block_i + self.numCoeffs, block_j:block_j + self.numCoeffs]

						A_real[block_i:block_i + self.numCoeffs, block_j:block_j + self.numCoeffs] = np.real(self.S.dot(M_complex.dot(self.S_inv)))
						b_real[block_i:block_i + self.numCoeffs] = np.real(self.S.dot(b_complex[block_i:block_i + self.numCoeffs]))





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