import numpy as np
import shtools
import util
import lspn
import cas
import pnbuilder
import scipy.io




class Constant(object):
    def __init__(self, value):
        self.value = value
    def __call__(self, x):
        return self.value
    def dx(self, x):
        return 0.0
    def dy(self, x):
        return 0.0
    def dz(self, x):
        return 0.0

class Gradient(object):
    def __init__(self, normal):
        norm=np.linalg.norm(normal)
        if norm==0: 
            self.normal = normal
        else:
            self.normal = normal/norm
    def __call__(self, x):
        return np.dot(x, self.normal)
    def dx(self, x):
        return self.normal[0]
    def dy(self, x):
        return self.normal[1]
    def dz(self, x):
        return 0.0
        #return self.normal[2]

    
class RBF(object):
	def __init__(self, stddev, amplitude=1.0 ):
		self.amplitude = amplitude
		self.stddev = stddev
		self.variance = stddev*stddev
		self.normalization = 1.0/(np.power(stddev, 3.0)*np.power(2.0*np.pi, 1.5))
	def __call__(self, x):
		#return self.normalization*np.exp(-(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])/self.variance)
		return self.normalization*np.exp(-(x[0]*x[0]+x[1]*x[1])/self.variance)
	def dx(self, x):
		return -self(x)*2.0/(self.variance)*x[0]
	def dxdx(self, x):
		c = 2.0/self.variance
		return self(x)*(c*c*x[0]*x[0] - c)
	def dxdy(self, x):
		c = 2.0/self.variance
		return self(x)*c*c*x[0]*x[1]
	def dy(self, x):
		return -self(x)*2.0/(self.variance)*x[1]
	def dydx(self, x):
		return self.dxdy(x)
	def dydy(self, x):
		c = 2.0/self.variance
		return self(x)*(c*c*x[1]*x[1] - c)
	def dz(self, x):
		return 0.0


class SHEXP(object):
    def __init__(self, order, coeff_functions):
        self.order = order
        self.coeff_functions = coeff_functions
    def __call__(self, x, omega):
        (theta, phi) = shtools.sphericalCoordinates(omega)
        coeffs = [f(x) for f in self.coeff_functions]
        return shtools.sh_sum(self.order, coeffs, theta, phi)
    def dx(self, x, omega):
    	(theta, phi) = shtools.sphericalCoordinates(omega)
    	coeffs_dx = [f.dx(x) for f in self.coeff_functions]
    	return shtools.sh_sum(self.order, coeffs_dx, theta, phi)
    def dxdx(self, x, omega):
    	(theta, phi) = shtools.sphericalCoordinates(omega)
    	coeffs_dxdx = [f.dxdx(x) for f in self.coeff_functions]
    	return shtools.sh_sum(self.order, coeffs_dxdx, theta, phi)
    def dxdy(self, x, omega):
    	(theta, phi) = shtools.sphericalCoordinates(omega)
    	coeffs_dxdy = [f.dxdy(x) for f in self.coeff_functions]
    	return shtools.sh_sum(self.order, coeffs_dxdy, theta, phi)
    def dydx(self, x, omega):
    	return self.dxdy(x, omega)
    def dy(self, x, omega):
    	(theta, phi) = shtools.sphericalCoordinates(omega)
    	coeffs_dy = [f.dy(x) for f in self.coeff_functions]
    	return shtools.sh_sum(self.order, coeffs_dy, theta, phi)
    def dydy(self, x, omega):
    	(theta, phi) = shtools.sphericalCoordinates(omega)
    	coeffs_dydy = [f.dydy(x) for f in self.coeff_functions]
    	return shtools.sh_sum(self.order, coeffs_dydy, theta, phi)
    def dz(self, x, omega):
    	(theta, phi) = shtools.sphericalCoordinates(omega)
    	coeffs_dz = [f.dz(x) for f in self.coeff_functions]
    	return shtools.sh_sum(self.order, coeffs_dz, theta, phi)
    def dzdz(self, x, omega):
    	(theta, phi) = shtools.sphericalCoordinates(omega)
    	coeffs_dzdz = [f.dzdz(x) for f in self.coeff_functions]
    	return shtools.sh_sum(self.order, coeffs_dzdz, theta, phi)

class Derivative(object):
    def __init__(self, fun, dimension = 0, stepsize = 0.01):
        self.fun = fun
        self.step = np.zeros(3)
        self.step[dimension] = stepsize
        self.stepsize = stepsize
    def __call__(self, x):
        return (self.fun( x+self.step )-self.fun( x-self.step ))/(2.0*self.stepsize)

def problem_debug( order = 1 ):

	sigma_t = Gradient(np.array([1.0, 1.0]))
	sigma_s = Gradient(np.array([1.0, 1.0]))

	functions = {}
	functions["\\sigma_t"] = lambda pWS: sigma_t(pWS)
	functions["\\sigma_s"] = lambda pWS: sigma_s(pWS)
	functions["sigma_t_instance"] = sigma_t
	functions["sigma_s_instance"] = sigma_s
	q_instance = SHEXP(order, [Constant(1.0), Constant(0.0), Constant(0.0), Constant(0.0)])
	def q_fun( l, m, pWS ):
		shindex = shtools.shIndex(l,m)
		if shindex is None or l > order or m < -l or m > l:
			return 0.0
		return q_instance.coeff_functions[shindex](pWS)
	functions["q"] = q_fun
	functions["q_instance"] = q_instance


	f_instance = SHEXP(order, [Constant(1.0), Constant(0.0), Constant(0.0), Constant(0.0)])
	def f_fun( l, m, pWS ):
		shindex = shtools.shIndex(l,m)
		if shindex is None or l > order or m < -l or m > l:
			return 0.0
		return f_instance.coeff_functions[shindex](pWS)
	functions["f_p"] = f_fun
	functions["f_instance"] = f_instance


	#functions["L"] = SHEXP(order, [Constant(1.0), Constant(0.2), Constant(0.3), Constant(0.4)])
	functions["L"] = SHEXP(order, [RBF(1.0), RBF(0.2), RBF(0.3), RBF(0.4)])
	return functions


def term1_factorized( x, l, m, problem ):
    #return shtools.c_lm(l-1, m-1)*0.5*sigma_t.dx(x)*L.coeff_functions[shtools.shIndex(l-1, m-1)](x)
    sigma_t = problem["sigma_t_instance"]
    L = problem["L"]


    result = 0.0
    if not pnb.shIndex(l-1, m-1) is None:
	    result += shtools.c_lm(l-1, m-1)*0.5*sigma_t.dx(x)*L.coeff_functions[shtools.shIndex(l-1, m-1)](x)

    if not pnb.shIndex(l+1, m-1) is None:
    	result += -shtools.d_lm(l+1, m-1)*0.5*sigma_t.dx(x)*L.coeff_functions[shtools.shIndex(l+1, m-1)](x)

    if not pnb.shIndex(l-1, m+1) is None:
    	result += -shtools.e_lm(l-1, m+1)*0.5*sigma_t.dx(x)*L.coeff_functions[shtools.shIndex(l-1, m+1)](x)

    if not pnb.shIndex(l+1, m+1) is None:
    	result += shtools.f_lm(l+1, m+1)*0.5*sigma_t.dx(x)*L.coeff_functions[shtools.shIndex(l+1, m+1)](x)

    if not pnb.shIndex(l-1, m-1) is None:
    	result += -shtools.c_lm(l-1, m-1)*np.complex(0.0, 1.0)*0.5*sigma_t.dy(x)*L.coeff_functions[shtools.shIndex(l-1, m-1)](x)

    if not pnb.shIndex(l+1, m-1) is None:
    	result += shtools.d_lm(l+1, m-1)*np.complex(0.0, 1.0)*0.5*sigma_t.dy(x)*L.coeff_functions[shtools.shIndex(l+1, m-1)](x)

    if not pnb.shIndex(l-1, m+1) is None:
    	result += -shtools.e_lm(l-1, m+1)*np.complex(0.0, 1.0)*0.5*sigma_t.dy(x)*L.coeff_functions[shtools.shIndex(l-1, m+1)](x)

    if not pnb.shIndex(l+1, m+1) is None:
    	result += shtools.f_lm(l+1, m+1)*np.complex(0.0, 1.0)*0.5*sigma_t.dy(x)*L.coeff_functions[shtools.shIndex(l+1, m+1)](x)

    if not pnb.shIndex(l-1, m) is None:
    	result += -shtools.a_lm(l-1, m)*sigma_t.dz(x)*L.coeff_functions[shtools.shIndex(l-1, m)](x)

    if not pnb.shIndex(l+1, m) is None:
    	result += -shtools.b_lm(l+1, m)*sigma_t.dz(x)*L.coeff_functions[shtools.shIndex(l+1, m)](x)

    return result

def term0(x, theta, phi, l, m, problem):
	omega = shtools.sphericalDirection(theta, phi)
	L = problem["L"]

	result = 0.0

	shindex = shtools.shIndex(l,m)
	if shindex is None or l<0 or m< -l or m > l:
		return 0.0

	Llm = L.coeff_functions[shindex]
	
	for i in range(2):
		for j in range(2):
			c = 0.0
			if i ==0 and j == 0:
				c = L.dxdx(x, omega)
			elif i ==0 and j == 1:
				c = L.dxdy(x, omega)
			elif i ==1 and j == 0:
				c = L.dydx(x, omega)
			elif i ==1 and j == 1:
				c = L.dydy(x, omega)
			else:
				raise ValueError("sdadasdasdsd")
			result += -omega[i]*omega[j]*c
	return result


def term1( x, theta, phi, l, m, problem ):
	sigma_t = problem["sigma_t_instance"]
	L = problem["L"]
	omega = shtools.sphericalDirection(theta, phi)
	return -(omega[0]*sigma_t.dx(x) + omega[1]*sigma_t.dy(x) + omega[2]*sigma_t.dz(x))*L(x, omega)

def term2( x, theta, phi, l, m, problem ):
	sigma_t = problem["sigma_t_instance"](x)
	L = problem["L"]
	omega = shtools.sphericalDirection(theta, phi)
	return sigma_t*sigma_t*L(x, omega)


def term3( x, theta, phi, l, m, problem ):
	# NB: assuming constant phase function....
	omega = shtools.sphericalDirection(theta, phi)
	sigma_t = problem["sigma_t_instance"]
	sigma_s = problem["sigma_s_instance"]
	f = problem["f_instance"]
	L00 = problem["L"].coeff_functions[0]

	result = 0.0
	result += omega[0]*(sigma_s.dx(x)*L00(x) + sigma_s(x)*L00.dx(x))
	result += omega[1]*(sigma_s.dy(x)*L00(x) + sigma_s(x)*L00.dy(x))
	result += omega[2]*(sigma_s.dz(x)*L00(x) + sigma_s(x)*L00.dz(x))
	result *= np.sqrt(4.0*np.pi)/(4.0*np.pi)

	return result



	#int_L = L.coeff_functions[0](x)*np.sqrt(4*np.pi)
	#return -sigma_t*sigma_s*(1.0/(4.0*np.pi))*int_L

def term4( x, theta, phi, l, m, problem ):
	# NB: assuming constant phase function....
	omega = shtools.sphericalDirection(theta, phi)
	sigma_t = problem["sigma_t_instance"](x)
	sigma_s = problem["sigma_s_instance"](x)
	f = problem["f_instance"]
	L = problem["L"]
	int_L = L.coeff_functions[0](x)*np.sqrt(4*np.pi)
	return -sigma_t*sigma_s*(1.0/(4.0*np.pi))*int_L

def term5( x, theta, phi, l, m, problem ):
	omega = shtools.sphericalDirection(theta, phi)
	q = problem["q_instance"]
	return -(omega[0]*q.dx(x, omega) + omega[1]*q.dy(x, omega) + omega[2]*q.dz(x, omega))

def term6( x, theta, phi, l, m, problem ):
	sigma_t = problem["sigma_t_instance"](x)
	q = problem["q_instance"]
	omega = shtools.sphericalDirection(theta, phi)
	return sigma_t*q(x, omega)


def test_term( order, domain, problem, term, term_expr ):
	voxel = np.array([2,2], dtype=int)
	pnb = pnbuilder.PNBuilder(order, domain)
	pnb.add_terms(term_expr)
	(A,b) = pnb.build_global( problem )

	# now construct solution vector from known L at the
	# respective coefficient locations
	# this step is actually redundant, because L is already defined in terms of SH coefficients
	# however, it serves as a sanity check...
	'''
	x_complex = np.zeros(domain.numVoxels*pnb.numCoeffs, dtype=complex)
	for voxel_i in range(domain.res_x):
		for voxel_j in range(domain.res_y):
			for i in range(pnb.numCoeffs):
				print("voxel_i={} voxel_j={} coeff={}".format(voxel_i, voxel_j, i))
				pWS = pnb.get_unknown_location( voxel_i, voxel_j, i ).getPWS()
				L = problem["L"]
				coeffs = shtools.project_sh(lambda theta, phi: L(pWS, shtools.sphericalDirection(theta, phi)), order)
				# NB: we take into account, that for 2d, pnb will have different index and lm ordering
				(l,m) = pnb.lmIndex(i)
				x_complex[pnb.get_global_index(voxel_i, voxel_j, i)] = coeffs[shtools.shIndex(l,m)]
	'''
	filename = "C:/projects/epfl/epfl17/python/sopn/data_debug_term.mat"
	x_complex = scipy.io.loadmat(filename)["x_complex"].reshape(domain.numVoxels*pnb.numCoeffs)

	#data = {}
	#data['x_complex'] = x_complex
	#data['A_complex'] = pnb.A_complex
	#scipy.io.savemat(filename, data)

	b_complex_pnb = pnb.A_complex.dot(x_complex)
	print("\n\n----------------------")

	#print("A=")
	#print(pnb.A_complex)
	#print("x=")
	#print(x_complex)
	#print("Ax=")
	#print()


	# groundtruth ----------------------------
	#print("\n\ngroundtruth ----------------------")

	# here we project the term under investigation into SH. This is done at the respective locations
	# of all the unknowns
	b_complex_gt = np.zeros(domain.numVoxels*pnb.numCoeffs, dtype=complex)
	for i in range(pnb.numCoeffs):
		voxel_i = voxel[0]
		voxel_j = voxel[1]
		(l,m) = pnb.lmIndex(i)
		pWS = pnb.get_unknown_location( voxel_i, voxel_j, i ).getPWS()
		global_i = pnb.get_global_index(voxel_i, voxel_j, i)

		debug = False
		if i == 1:
			#print(coeffs)
			debug = True

		coeffs = shtools.project_sh(lambda theta, phi: term(pWS, theta, phi, l, m, problem), order, debug)
		# NB: we take into account, that for 2d, pnb will have different index and lm ordering
		b_complex_gt[global_i] = coeffs[shtools.shIndex(l,m)]

	for i in range(pnb.numCoeffs):
		voxel_i = voxel[0]
		voxel_j = voxel[1]
		index = pnb.get_global_index(voxel_i, voxel_j, i)
		print( "{} {}".format(np.real(b_complex_gt[index]), np.real(b_complex_pnb[index])) )
		print( "{} {}".format(np.imag(b_complex_gt[index]), np.imag(b_complex_pnb[index])) )



def test_rhs_term( order, domain, problem, term, term_expr ):
	voxel = np.array([2,2], dtype=int)
	pnb = pnbuilder.PNBuilder(order, domain)
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
		pWS = pnb.get_unknown_location( voxel_i, voxel_j, i ).getPWS()
		global_i = pnb.get_global_index(voxel_i, voxel_j, i)
		coeffs = shtools.project_sh(lambda theta, phi: term(pWS, theta, phi, l, m, problem), order)
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




if __name__ == "__main__":
	order = 1
	domain = util.Domain2D(0.1, 5)

	problem = problem_debug()


	#test_term( order, domain, problem, term0, lspn.lspn_sotransport_term() )
	#test_term( order, domain, problem, term1, lspn.lspn_extinction_directional_derivative_term() )
	#test_term( order, domain, problem, term2, lspn.lspn_squared_extinction_term() )
	#test_term( order, domain, problem, term3, lspn.lspn_directional_derivative_scattering_term() )
	#test_term( order, domain, problem, term4, lspn.lspn_extinction_scattering_term() )
	#test_rhs_term( order, domain, problem, term5, lspn.lspn_directional_derivative_source_term() )
	#test_rhs_term( order, domain, problem, term6, lspn.lspn_extinction_source_term() )



