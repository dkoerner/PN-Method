import numpy as np
import shtools
import util
import lspn
import cas
import pnbuilder
import scipy.io
import problems



def term0(location, theta, phi, problem):
	omega = shtools.sphericalDirection(theta, phi)
	L = problem["L"]

	result = 0.0
	
	#'''
	for i in range(2):
		for j in range(2):
			c = 0.0
			if i ==0 and j == 0:
				c = L.dxdx(location, omega)
			elif i ==0 and j == 1:
				c = L.dxdy(location, omega)
			elif i ==1 and j == 0:
				c = L.dydx(location, omega)
			elif i ==1 and j == 1:
				c = L.dydy(location, omega)
			else:
				raise ValueError("sdadasdasdsd")
			result += -omega[i]*omega[j]*c
	#'''

	#result += omega[0]*L.dx(x, omega)
	#result += L.coeff_functions[shtools.shIndex(1,1)].dx(location)

	return result


def term1( location, theta, phi, problem ):
	sigma_t = problem["sigma_t_instance"]
	L = problem["L"]
	omega = shtools.sphericalDirection(theta, phi)
	return -(omega[0]*sigma_t.dx(x) + omega[1]*sigma_t.dy(x) + omega[2]*sigma_t.dz(x))*L(x, omega)

def term2( location, theta, phi, problem ):
	sigma_t = problem["sigma_t_instance"](x)
	L = problem["L"]
	omega = shtools.sphericalDirection(theta, phi)
	return sigma_t*sigma_t*L(x, omega)
	#return sigma_t*sigma_t


def term3( location, theta, phi, problem ):
	# NB: assuming constant phase function....
	omega = shtools.sphericalDirection(theta, phi)
	sigma_t = problem["sigma_t_instance"]
	sigma_s = problem["sigma_s_instance"]
	#f = problem["f_instance"]
	L00 = problem["L"].coeff_functions[0]

	result = 0.0
	result += omega[0]*(sigma_s.dx(x)*L00(x) + sigma_s(x)*L00.dx(x))
	result += omega[1]*(sigma_s.dy(x)*L00(x) + sigma_s(x)*L00.dy(x))
	result += omega[2]*(sigma_s.dz(x)*L00(x) + sigma_s(x)*L00.dz(x))
	result *= np.sqrt(4.0*np.pi)/(4.0*np.pi)

	return result



	#int_L = L.coeff_functions[0](x)*np.sqrt(4*np.pi)
	#return -sigma_t*sigma_s*(1.0/(4.0*np.pi))*int_L

def term4( location, theta, phi, problem ):
	# NB: assuming constant phase function....
	omega = shtools.sphericalDirection(theta, phi)
	sigma_t = problem["sigma_t_instance"](x)
	sigma_s = problem["sigma_s_instance"](x)
	#f = problem["f_instance"]
	L = problem["L"]
	int_L = L.integral_over_solid_angle(x)
	return -sigma_t*sigma_s*(1.0/(4.0*np.pi))*int_L

def term5( location, theta, phi, problem ):
	omega = shtools.sphericalDirection(theta, phi)
	q = problem["q_instance"]
	return -(omega[0]*q.dx(x, omega) + omega[1]*q.dy(x, omega) + omega[2]*q.dz(x, omega))

def term6( location, theta, phi, problem ):
	sigma_t = problem["sigma_t_instance"](x)
	q = problem["q_instance"]
	omega = shtools.sphericalDirection(theta, phi)
	return sigma_t*q(x, omega)



def test_term( label, problem, voxel, term, term_expr ):

	#print("\n\ntesting term {}----------------------".format(label))
	staggered = problem["staggered"]
	order = problem["order"]
	domain = problem["domain"]
	pnb = pnbuilder.PNBuilder(order, domain)

	if staggered == True:
		# staggered grid (and 3 point stencil)
		pnb.place_unknown( 0, (1,1) )
		pnb.place_unknown( 1, (1,0) )
		pnb.place_unknown( 2, (0,1) )
		pnb.set_stencil_half_steps(1)

	pnb.add_terms(term_expr)
	(A,b) = pnb.build_global( problem )


	if not "x_complex" in problem:
		# now construct solution vector from known L at the
		# respective coefficient locations
		# this step is actually redundant, because L is already defined in terms of SH coefficients
		# however, it serves as a sanity check...
		x_complex = np.zeros(domain.numVoxels*pnb.numCoeffs, dtype=complex)
		for voxel_i in range(domain.res_x):
			for voxel_j in range(domain.res_y):
				for i in range(pnb.numCoeffs):
					print("voxel_i={} voxel_j={} coeff={}".format(voxel_i, voxel_j, i))
					location = pnb.get_unknown_location( voxel_i, voxel_j, i )
					L = problem["L"]
					coeffs = shtools.project_sh(lambda theta, phi: L(location, shtools.sphericalDirection(theta, phi)), order)
					# NB: we take into account, that for 2d, pnb will have different index and lm ordering
					(l,m) = pnb.lmIndex(i)
					x_complex[pnb.get_global_index(voxel_i, voxel_j, i)] = coeffs[shtools.shIndex(l,m)]

		data = {}
		data['x_complex'] = x_complex
		filename = "C:/projects/epfl/epfl17/python/sopn/data_debug_term.mat"
		scipy.io.savemat(filename, data)

		problem["x_complex"] = x_complex
	else:
		x_complex = problem["x_complex"]

	b_complex_pnb = pnb.A_complex.dot(x_complex)

	# groundtruth ----------------------------
	#print("\n\ngroundtruth ----------------------")

	# here we project the term under investigation into SH. This is done at the respective locations
	# of all the unknowns
	b_complex_gt = np.zeros(domain.numVoxels*pnb.numCoeffs, dtype=complex)
	for i in range(pnb.numCoeffs):
		voxel_i = voxel[0]
		voxel_j = voxel[1]
		(l,m) = pnb.lmIndex(i)
		location = pnb.get_unknown_location( voxel_i, voxel_j, i )
		global_i = pnb.get_global_index(voxel_i, voxel_j, i)
		# NB: we take into account, that for 2d, pnb will have different index and lm ordering
		coeff = shtools.project_sh_coeff(lambda theta, phi: term(location, theta, phi, problem), l, m)
		b_complex_gt[global_i] = coeff

	for i in range(pnb.numCoeffs):
		voxel_i = voxel[0]
		voxel_j = voxel[1]
		index = pnb.get_global_index(voxel_i, voxel_j, i)
		print( "{} {}".format(np.real(b_complex_gt[index]), np.real(b_complex_pnb[index])) )
		print( "{} {}".format(np.imag(b_complex_gt[index]), np.imag(b_complex_pnb[index])) )



def test_rhs_term( label, problem, voxel, term, term_expr ):
	print("\n\ntesting term {}----------------------".format(label))
	staggered = problem["staggered"]
	order = problem["order"]
	domain = problem["domain"]
	pnb = pnbuilder.PNBuilder(order, domain)

	if staggered == True:
		# staggered grid (and 3 point stencil)
		pnb.place_unknown( 0, (1,1) )
		pnb.place_unknown( 1, (1,0) )
		pnb.place_unknown( 2, (0,1) )
		pnb.set_stencil_half_steps(1)


	pnb.add_terms(term_expr)
	(A,b) = pnb.build_global( problem )

	b_complex_pnb = pnb.b_complex

	# groundtruth ----------------------------
	#print("\n\ngroundtruth ----------------------")

	# here we project the term under investigation into SH. This is done at the respective locations
	# of all the unknowns
	b_complex_gt = np.zeros(domain.numVoxels*pnb.numCoeffs, dtype=complex)
	for i in range(pnb.numCoeffs):
		voxel_i = voxel[0]
		voxel_j = voxel[1]
		(l,m) = pnb.lmIndex(i)
		location = pnb.get_unknown_location( voxel_i, voxel_j, i )
		global_i = pnb.get_global_index(voxel_i, voxel_j, i)
		coeffs = shtools.project_sh(lambda theta, phi: term(location, theta, phi, problem), order)
		# NB: we take into account, that for 2d, pnb will have different index and lm ordering
		b_complex_gt[global_i] = coeffs[shtools.shIndex(l,m)]

	
	#print("x=")
	#print(x_complex)


	#print("Ax=")
	#print(b_complex)

	for i in range(pnb.numCoeffs):
		voxel_i = voxel[0]
		voxel_j = voxel[1]
		index = pnb.get_global_index(voxel_i, voxel_j, i)
		print( "{} {}".format(np.real(b_complex_gt[index]), np.real(b_complex_pnb[index])) )
		print( "{} {}".format(np.imag(b_complex_gt[index]), np.imag(b_complex_pnb[index])) )


def print_expr(expr, switch = True):
	#pass
	if switch == True:
		print("\n----------------------------\n")
		print("$$\n" + cas.latex(expr) + "\n$$")



def debug_terms_simple():
	problem = problems.debug()

	#problem = problems.checkerboard2d_rasterized()
	#filename = "C:/projects/epfl/epfl17/python/sopn/dat_checkerboard2d_rasterized.mat"
	#data = {}
	#data["sigma_t"] = problem["\\sigma_t"].voxels
	#scipy.io.savemat(filename, data)


	#'''
	filename = "C:/projects/epfl/epfl17/python/sopn/data_debug_term_staggered.mat"
	x_complex_shape = scipy.io.loadmat(filename)["x_complex"].shape
	problem["x_complex"] = scipy.io.loadmat(filename)["x_complex"].reshape(x_complex_shape[0]*x_complex_shape[1])

	voxel = np.array([2,2], dtype=int)

	test_term( "term0", problem, voxel, term0, lspn.lspn_sotransport_term() )
	test_term( "term1", problem, voxel, term1, lspn.lspn_extinction_directional_derivative_term() )
	test_term( "term2", problem, voxel, term2, lspn.lspn_squared_extinction_term() )
	test_term( "term3", problem, voxel, term3, lspn.lspn_directional_derivative_scattering_term() )
	test_term( "term4", problem, voxel, term4, lspn.lspn_extinction_scattering_term() )
	test_rhs_term( "term5", problem, voxel, term5, lspn.lspn_directional_derivative_source_term() )
	test_rhs_term( "term6", problem, voxel, term6, lspn.lspn_extinction_source_term() )
	#'''	


class Test(object):
	def __init__( self, location, problem, term ):
		self.location = location
		self.problem = problem
		self.term = term
	def __call__( self, angle ):
		return self.term(self.location, angle[0], angle[1], self.problem)
		#return 1.0


def debug_A_term(debug_voxel_x, debug_voxel_y):

	terms = []
	terms.append((0, lspn.lspn_sotransport_term()))
	#terms.append((1, lspn_extinction_directional_derivative_term()))
	#terms.append((2, lspn_squared_extinction_term()))
	#terms.append((3, lspn_directional_derivative_scattering_term2()))
	#terms.append((4, lspn_extinction_scattering_term()))
	#terms.append((5, lspn_directional_derivative_source_term()))
	#terms.append((6, lspn_extinction_source_term()))

	for (term_index, term) in terms:
		problem = problems.blurred(problems.checkerboard2d(), 10.0)

		problem["debug_voxel_x"] = debug_voxel_x
		problem["debug_voxel_y"] = debug_voxel_y

		pnb = pnbuilder.PNBuilder(problem["order"], problem["domain"])

		# staggered grid (and 3 point stencil)
		if problem["staggered"] == True:
			pnb.place_unknown( 0, (1,1) )
			pnb.place_unknown( 1, (1,0) )
			pnb.place_unknown( 2, (0,1) )
			pnb.set_stencil_half_steps(1)

		pnb.add_terms(term)

		problem["id"] += "_term{}".format(term_index)

		(A,b) = pnb.build_global( problem )

		#data = {}
		#data['A'] = A
		#data['b'] = b.reshape((pnb.domain.numVoxels*pnb.numCoeffs, 1))
		#scipy.io.savemat("C:/projects/epfl/epfl17/python/sopn/debug_terms/debug_A_term{}.mat".format(term_index), data)
		return A





if __name__ == "__main__":
	#lspn.lspn_sotransport_term(True)
	#exit(1)

	#debug_terms_simple()

	id = "checkerboard_blur10.0_term1"
	data = lspn.load_pn_solution("C:/projects/epfl/epfl17/python/sopn/solution_{}.mat".format(id))
	pnb = data["pnb"]
	x_real = data["x_real"]

	problem = {}

	# L -----
	x_complex = pnb.to_complex(x_real)
	coefficient_grids = [problems.Constant(0.0) for i in range(shtools.numSHCoeffs(pnb.N))]
	for i in range(pnb.numCoeffs):
		(l,m) = pnb.lmIndex(i)
		offset = pnb.get_unknown_offset(i)
		grid = problems.CoefficientGrid(pnb.domain, pnb.numCoeffs, i, offset, x_complex )
		coefficient_grids[shtools.shIndex(l,m)] = grid
	problem["L"] = problems.SHEXP( pnb.N, coefficient_grids )

	# sigma_t -----
	problem["sigma_t_instance"] = problems.VoxelGrid(pnb.domain, data["sigma_t"])
	problem["sigma_t_instance"].build_derivatives()

	problem["sigma_s_instance"] = problems.VoxelGrid(pnb.domain, data["sigma_s"])
	problem["sigma_s_instance"].build_derivatives()



	terms = []
	terms.append((0, term0))
	#terms.append((1, term1))
	#terms.append((2, term2))
	#terms.append((3, term3))
	#terms.append((4, term4))
	##terms.append(term5)
	##terms.append(term6)


	debug_voxel_x = 40
	debug_voxel_y = 35
	voxels_min = np.array([debug_voxel_x, debug_voxel_y])
	voxels_max = np.array([debug_voxel_x+1, debug_voxel_y+1])

	#voxels_min = np.array([20, 20])
	#voxels_max = np.array([50, 50])


	data = {}
	
	#'''
	for (term_index, term) in terms:

		# get the matrix A for the current term
		A_term = scipy.io.loadmat("C:/projects/epfl/epfl17/python/sopn/debug_terms/debug_A_term{}.mat".format(term_index))["A"]

		#A_term = debug_A_term(debug_voxel_x, debug_voxel_y)
		x_term_real = A_term.dot(x_real)


		x_term_complex_gt = np.zeros(pnb.domain.numVoxels*pnb.numCoeffs, dtype=complex)
		for voxel_i in range(voxels_min[0], voxels_max[0]):
			print("voxel_i={}".format(voxel_i))
			for voxel_j in range(voxels_min[1], voxels_max[1]):
				# here we project the term under investigation into SH. This is done at the respective locations
				# of all the unknowns
				#for i in range(pnb.numCoeffs):
				for i in range(1):
					(l,m) = pnb.lmIndex(i)
					location = pnb.get_unknown_location(voxel_i, voxel_j, i)
					global_i = pnb.get_global_index(voxel_i, voxel_j, i)
					# NB: we take into account, that for 2d, pnb will have different index and lm ordering
					#print("---------------")
					#coeff = shtools.project_sh_coeff(lambda theta, phi: term(location, theta, phi, problem), l, m)
					#coeff = shtools.project_sh_coeff_x(lambda angle: term(pWS, angle[0], angle[1], problem), l, m)
					coeff = shtools.project_sh_coeff_x(Test(location, problem, term), l, m)
					#print(coeff)
					#print("done")
					x_term_complex_gt[global_i] = coeff
		x_term_real_gt = pnb.to_real(x_term_complex_gt)
		debug_i = pnb.get_global_index(debug_voxel_x,debug_voxel_y,0)
		print(x_term_real[debug_i])
		print(x_term_real_gt[debug_i])
		#print(pnb.to_complex(x_term_real)[debug_i])
		#print(x_term_complex_gt[debug_i])

		#data['x_term{}'.format(term_index)] = x_term_real_gt.reshape((pnb.domain.numVoxels*pnb.numCoeffs, 1))
		#scipy.io.savemat("C:/projects/epfl/epfl17/python/sopn/debug_terms/data_term{}.mat".format(term_index), data)
	#'''

	'''
	num = 1000
	sigma_t = np.zeros((num,num))
	for i in range(num):
		for j in range(num):
			x = i/float(num-1)*7.0
			y = j/float(num-1)*7.0
			pWS = np.array([x, y])
			sigma_t[i, j] = problem["sigma_t_instance"](pWS)
	data = {"sigma_t":sigma_t}
	scipy.io.savemat("C:/projects/epfl/epfl17/python/sopn/debug_terms/interp_test.mat", data)
	'''

