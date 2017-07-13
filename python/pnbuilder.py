
import numpy as np
import util
import time
import shtools
import scipy.io
import cas




'''
def find_lm_offset(Llm):
	# find the l,m arguments
	# this is tricky as we need to inspect the expressions
	# alternatively we could use the latex code as hash, but
	# I somehow didn't like it really
	l = Llm.getArgument(0)
	l_offset = 0
	m = Llm.getArgument(1)
	m_offset = 0

	if l.__class__ == cas.Variable:
		pass
	elif l.__class__ == cas.Addition and l.getOperand(1).__class__ == cas.Number:
		l_offset = l.getOperand(1).getValue()
	elif l.__class__ == cas.Addition and l.getOperand(1).__class__ == cas.Negate and l.getOperand(1).getOperand(0).__class__ == cas.Number:
		l_offset = -l.getOperand(1).getOperand(0).getValue()
	else:
		raise ValueError


	if m.__class__ == cas.Variable:
		pass
	elif m.__class__ == cas.Addition and m.getOperand(1).__class__ == cas.Number:
		m_offset = m.getOperand(1).getValue()
	elif m.__class__ == cas.Addition and m.getOperand(1).__class__ == cas.Negate and m.getOperand(1).getOperand(0).__class__ == cas.Number:
		m_offset = -m.getOperand(1).getOperand(0).getValue()
	else:
		raise ValueError

	return (l_offset, m_offset)

class Unknown(object):
	def __init__(self, l_offset, m_offset):
		self.l_offset = l_offset
		self.m_offset = m_offset
	def getId(self):
		return "u"
	def __str__(self):
		return "L l_offset={} m_offset={}\n".format(self.l_offset, self.m_offset)
	def __eq__(self, other):
		if other.__class__ == Unknown and self.l_offset == other.l_offset and self.m_offset == other.m_offset:
			return True
		return False
	def __ne__(self, other):
			return not self == other

class Coefficients(object):
	def __init__(self, *args):
		self.coeffs = list(args)
	def getId(self):
		return "c"
	def append_coefficient(self, expr):
		self.coeffs.append(expr)
	def __str__(self):
		result = "coefficients:\n"
		for c in self.coeffs:
			result += "\t" + str(c) + "\n"
		return result
	def __ne__(self, other):
		return not self == other
	def __eq__(self, other):
		if other.__class__ == Coefficients:
			return True


class Derivation(object):
	def __init__(self, var):
		self.var = var
	def getId(self):
		return "d{}".format(self.var)
	def __str__(self):
		return "Derivation in {}\n".format(self.var)
	def __ne__(self, other):
		return not self == other
	def __eq__(self, other):
		if other.__class__ == Derivation and self.var == other.var:
			return True


class CoefficientChain(object):
	def __init__(self, l_offset, m_offset):
		self.chain = [Unknown(l_offset, m_offset)]
	def numElements(self):
		return len(self.chain)
	def getElements(self):
		return self.chain
	def append_coefficient(self, expr):
		if self.chain[-1].__class__ != Coefficients:
			self.chain.append(Coefficients())
		self.chain[-1].append_coefficient(expr)
	def derive(self, var):
		self.chain.append(Derivation(var))

	def getId(self):
		result = ""
		for c in reversed(self.chain):
			result += c.getId()
		return result
	def __str__(self):
		result = ""
		for e in self.chain:
			result += str(e)
		return result
	def __getitem__(self, index):
		return self.chain[index]

class ChainElementData(object):
	def __init__(self, id):
		self.coefficients = {}
		self.id = id
	def getId(self):
		return self.id
	def append_coefficients(self, l_offset, m_offset, coeffs):
		key = (l_offset, m_offset)
		if not key in self.coefficients:
			self.coefficients[key] = []
		self.coefficients[key].append(coeffs)




def extract_coefficients( parent, unknown_symbol = "L" ):
	coeff_chain = None
	others = []
	for i in range(parent.numChildren()):
		child = parent.getChildren(i)
		result = extract_coefficients(child, unknown_symbol)

		if result == None:
			others.append(child)
		else:
			coeff_chain = result

	# now handle the coefficients to the unknown ---

	# no coefficient found yet in any of our childs?
	if coeff_chain == None:
		#print("{} {}".format(unknown_symbol,parent.getSymbol() ))
		# check if parent is the unkown we are looking for
		if parent.__class__ == cas.SHCoefficient and parent.getSymbol() == unknown_symbol:
			(l_offset, m_offset) = find_lm_offset(parent)
			coeff_chain = CoefficientChain( l_offset, m_offset )
			return coeff_chain
		else:
			# non of our childs found Llm nor is it the parent
			# so we are looking in the wrong branch of the expression tree
			return None
	else:
		# one of our children contained the unkown
		# we now analyse its siblings in the expression tree
		# to find the coefficients
		if parent.__class__ == cas.Negate:
			coeff_chain.append_coefficient(cas.Number(-1))
			return coeff_chain
		elif parent.__class__ == cas.Multiplication:
			for op in others:
				coeff_chain.append_coefficient(op)
			return coeff_chain
		elif parent.__class__ == cas.Derivation:
			coeff_chain.derive(parent.getVariable().getSymbol())
			return coeff_chain
'''



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



def eval_term_recursive( expr, info ):

	# in case of a simple number
	if expr.__class__ == cas.Number:
		return expr.getValue()
	elif expr.__class__ == cas.Negate:
		return -eval_term_recursive(expr.getOperand(0), info)
	#elif expr.__class__ == cas.Variable:
	elif isinstance(expr, cas.Variable):
		# TODO: if we come across the unknown, then return a stencil point
		if expr.getSymbol() in info.vars:
			return info.vars[expr.getSymbol()]
		else:
			raise ValueError("unable to resolve variable {}".format(expr.getSymbol()))
	# in case of a function
	elif isinstance(expr, cas.Function):
		numArgs = expr.numArguments()
		# evaluate all arguments
		args = []
		for i in range(numArgs):
			args.append( eval_term_recursive(expr.getArgument(i), info) )
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
		raise ValueError("function {} not defined for evaluation".format(expr.getSymbol()))
	# derivation...
	elif expr.__class__ == cas.Derivation:
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
		step[dimension] = 1
		stepWS = step*info.voxelsize
		central_difference_weight = 1.0/(2.0*info.voxelsize[dimension])

		nested_expr = expr.getExpr()

		# idea behind this is to evaluate the child expression for the different positions
		# of the discretization stencils
		info.vars['\\vec{x}'] = pWS - stepWS
		info.coord = coord - step
		a = eval_term_recursive(nested_expr, info)

		info.vars['\\vec{x}'] = pWS + stepWS
		info.coord = coord + step
		b = eval_term_recursive(nested_expr, info)

		info.vars['\\vec{x}'] = pWS
		return central_difference_weight*(b - a)
	elif expr.__class__ == cas.Multiplication:
		numOperands = expr.numOperands()
		result = 1.0
		for i in range(numOperands):
			result = result * eval_term_recursive(expr.getOperand(i), info)
		return result
	else:
		raise ValueError("unable to handle expression of type {}".format(expr.__class__.__name__))


def eval_term(expr, unknown_symbol, vars, funs, coord, voxelsize):
	info = EvalInfo()
	info.unknown_symbol = unknown_symbol
	info.vars = vars
	info.functions = funs
	info.voxelsize = voxelsize
	info.coord = coord
	return eval_term_recursive(expr, info)
	



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
		self.lhs_term_table = {}
		self.lhs_terms = []
		self.rhs_terms = []
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


	def add_lhs(self, expr):
		unknown_symbol = "L"


		terms = []
		if expr.__class__ == cas.Multiplication or expr.__class__ == cas.Negate:
			terms.append(expr)
		elif expr.__class__ == cas.Addition:
			numTerms = expr.numOperands()
			for i in range(numTerms):
				terms.append(expr.getOperand(i))
		else:
			ValueError("expected expression to be addition or multiplication")


		# we iterate over all operands of the outer sum
		for m in terms:
			# extract coefficients and their structure
			term_chain = extract_coefficients(m, unknown_symbol = unknown_symbol)

			term_elements = reversed(term_chain.getElements())
			for e in term_elements:
				if e.getId() == "u":
					# todo: mark as lhs term
					pass
				elif e.getId() == "c":
					# finalize all coefficients
					if len(e.coeffs) == 0:
						raise ValueError("expected at least one coefficient expression")
					elif len(e.coeffs) == 1:
						e.coeff_expr = e.coeffs[0]
					else:
						e.coeff_expr = cas.Multiplication(e.coeffs)

			self.lhs_terms.append(term_chain)

	def add_rhs(self, expr):

		terms = []
		if expr.__class__ == cas.Multiplication or expr.__class__ == cas.Negate:
			terms.append(expr)
		elif expr.__class__ == cas.Addition:
			numTerms = expr.numOperands()
			for i in range(numTerms):
				terms.append(expr.getOperand(i))
		else:
			ValueError("expected expression to be addition or multiplication")


		# we iterate over all operands of the outer sum
		for m in terms:
			self.rhs_terms.append(m)


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
		A = np.zeros( (numVoxels*self.numCoeffs, numVoxels*self.numCoeffs) )
		b = np.zeros( (numVoxels*self.numCoeffs) )

		# Now we build the global coefficient matrix A and rhs b by iterating over all elements (voxels)
		# and within each voxel, we iterate over all sh coefficients. Each row within the global system
		# is associated with a specific sh coefficient of a specific voxel.
		# Then, for each row, we evaluate all terms and accumulate the results into A and b.
		for voxel_x in range(self.domain.res_x):
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
					for term in self.terms:
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

						if result.__class__ == UnknownInfo:
							# this is a lhs term
							# weights are going into A
							for u in result.unknowns:
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

								A[global_i, global_j] += u.weight
						else:
							# this is a rhs term (because it has no unknowns and evaluated to a number)
							#this goes straight into b
							b[global_i] += result



						'''
						# find which coefficient the current term is dealing with
						# for that we know that the first term element contains information
						# about the unknown, such as the l,m offset
						local_j = self.shIndex( l+term[0].l_offset, m+term[0].m_offset )

						# check bounds
						if local_j == None or local_j >=self.numCoeffs:
							continue


						stencil = util.stencil2d(voxel_x, voxel_y)

						term_elements = reversed(term.getElements())
						for e in term_elements:
							if e.getId() == "u":
								# we reached the unknown and are done
								# now accumulate the coefficients
								for p in stencil.getPoints():
									# position and shcoeff-index define the final column within the current row
									# of the global matrix A

									# TODO: think about how to handle boundaries
									if p.coord[0] < 0 or p.coord[0] >= self.domain.res_x or p.coord[1] < 0 or p.coord[1] >= self.domain.res_y:
										continue

									global_j = self.get_global_index(p.coord[0], p.coord[1], local_j)

									A[global_i, global_j] += p.weight
							elif e.getId() == "c":
								# a coefficient expression which we want to evalute to a weight for our
								# stencil points
								for p in stencil.getPoints():
									# get world position of current stencil point
									pWS = np.array([self.X[p.coord[0], p.coord[1]], self.Y[p.coord[0], p.coord[1]]])
									# here we evaluate the expression to arrive at a concrete number
									variables = {"l'":l, "m'":m, "x":pWS[0], "y":pWS[1], "z":0.0, "\\vec{x}":pWS}
									functions = {"\\sigma_t":lambda pWS: sigma_a(pWS) + sigma_s(pWS), "\\sigma_s":sigma_s, "f_p":phase_shcoeffs}
									p.weight *= cas.eval(e.coeff_expr, variables, functions)
							elif e.getId() == "dx":
								# discretize derivative in x
								stencil.dx(domain.h_x)
							elif e.getId() == "dy":
								# discretize derivative in y
								stencil.dy(domain.h_y)
							elif e.getId() == "dz":
								# discretize derivative in z
								stencil.dz(domain.h_z)
							else:
								raise ValueError("unknown term element")
						'''
					## evaluate rhs terms -------
					#for term in self.rhs_terms:
					#	pass

		return (A,b)








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