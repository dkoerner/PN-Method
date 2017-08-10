
import numpy as np
import util
import shtools
import cas
import pnbuilder
from scipy.ndimage.filters import gaussian_filter
import scipy.io







class Constant(object):
    def __init__(self, value):
        self.value = value
    def __call__(self, x):
        return self.value
    def dx(self, x):
        return 0.0
    def dxdx(self, x):
        return 0.0
    def dxdy(self, x):
        return 0.0
    def dy(self, x):
        return 0.0
    def dydy(self, x):
        return 0.0
    def dydx(self, x):
        return 0.0
    def dz(self, x):
        return 0.0

class Sum(object):
	def __init__(self, a, b):
		self.a = a
		self.b = b
	def __call__(self, x):
		return self.a(x) + self.b(x)
	def dx(self, x):
		return self.a.dx(x) + self.b.dx(x)
	def dy(self, x):
		return self.a.dy(x) + self.b.dy(x)
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
    def integral_over_solid_angle(self, x):
    	return self.coeff_functions[0](x)*np.sqrt(4*np.pi)

class VoxelGrid(object):
	def __init__( self, domain, voxels, offset = np.array([0.5, 0.5]) ):
		self.domain = domain
		self.voxels = voxels
		self.offset = offset



	def sample( self, voxel ):
		return self.voxels[voxel[0], voxel[1]]
		
	def eval( self, pVS ):
		voxel = np.floor(pVS).astype(int)

		# clamp to domain boundaries
		if voxel[0] < 0:
			voxel[0] = 0
		elif voxel[0] >= self.domain.res_x:
			voxel[0] = self.domain.res_x-1
		if voxel[1] < 0:
			voxel[1] = 0
		elif voxel[1] >= self.domain.res_y:
			voxel[1] = self.domain.res_y-1

		return self.sample(voxel)

		
	def __call__(self, pWS):
		pVS = self.domain.worldToVoxel(pWS)
		# take offset into account
		pVS -= self.offset
		return self.eval(pVS)

	def create_derivative_grid(self, dimension = 0):
		step = np.zeros(2, dtype=int)
		step[dimension] = 1
		voxels = np.zeros( (self.domain.res_x, self.domain.res_y), dtype=self.voxels.dtype )
		for voxel_i in range(self.domain.res_x):
			for voxel_j in range(self.domain.res_y):
				c = np.array([voxel_i, voxel_j])
				p = c + step
				m = c - step
				if not self.in_bounds(p):
					# forward differences
					voxels[voxel_i, voxel_j] = (self.sample(c) - self.sample(m))/self.domain.voxelsize[dimension]
				elif not self.in_bounds(m):
					# backward differences
					voxels[voxel_i, voxel_j] = (self.sample(p) - self.sample(c))/self.domain.voxelsize[dimension]
				else:
					# finite differences
					voxels[voxel_i, voxel_j] = (self.sample(p) - self.sample(m))/(2.0*self.domain.voxelsize[dimension])
		return VoxelGrid(self.domain, voxels, self.offset)

	def in_bounds(self, voxel):
		if voxel[0] < 0 or voxel[0]>=self.domain.res_x or voxel[1] < 0 or voxel[1]>=self.domain.res_y:
			return False
		return True


	def blur( self, stddev ):
		self.voxels = gaussian_filter(self.voxels, stddev)

	def build_derivatives(self):
		# precompute dx, dy, dz grids so that we just have to look it up later
		self.dx = self.create_derivative_grid(0)
		self.dxdx = self.dx.create_derivative_grid(0)
		self.dxdy = self.dx.create_derivative_grid(1)
		self.dy = self.create_derivative_grid(1)
		self.dydy = self.dy.create_derivative_grid(1)
		self.dydx = self.dxdy
		self.dz = Constant(0.0)


class CoefficientGrid(VoxelGrid):
	def __init__(self, domain, stride, index, offset, x):
		self.domain = domain
		self.stride = stride

		voxels = np.zeros( (domain.res_x, domain.res_y), dtype=complex )
		for voxel_i in range(domain.res_y):
			for voxel_j in range(domain.res_x):
				voxels[voxel_i, voxel_j] = x[self.get_global_index(voxel_i, voxel_j, index)]

		super().__init__(domain, voxels, offset)

		self.build_derivatives()



	def get_global_index( self, voxel_i, voxel_j, coeff_index ):
		v = voxel_j*self.domain.res_x + voxel_i
		return v*self.stride + coeff_index

	'''
	def dx(self, x):
		return self.dx_grid(x)
	def dxdx(self, x):
		return self.dxdx_grid(x)
	def dxdy(self, x):
		return self.dxdy_grid(x)
	def dy(self, x):
		return self.dy_grid(x)
	def dydx(self, x):
		return self.dxdy_grid(x)
	def dydy(self, x):
		return self.dydy_grid(x)
	def dz(self, x):
		#return self.dz_grid(x)
		return 0.0
	'''


class Derivative(object):
    def __init__(self, fun, dimension = 0, stepsize = 0.01):
        self.fun = fun
        self.step = np.zeros(3)
        self.step[dimension] = stepsize
        self.stepsize = stepsize
    def __call__(self, x):
        return (self.fun( x+self.step )-self.fun( x-self.step ))/(2.0*self.stepsize)












def blob2d():
	# define problem ---
	def source( pWS ):
		x = pWS[0]
		y = pWS[1]
		if x > 3.0 and x < 4.0 and y > 3.0 and y < 4.0:
			return 3.0
		return 0.0

	def sigma_a( pWS ):
		x = pWS[0]
		y = pWS[1]
		if x > 1.0 and x < 6.0 and y > 1.0 and y < 6.0:
			return 0.0
		return 3.0

	def sigma_s( pWS ):
		x = pWS[0]
		y = pWS[1]
		if x > 1.0 and x < 6.0 and y > 1.0 and y < 6.0:
			return 3.0
		return 0.0


	def phase_shcoeffs( l, m, pWS ):
		if l == 0:
			return 1.0
		return 0.0

	def source_shcoeffs( l, m, pWS ):
		if l==0:
			x = pWS[0]
			y = pWS[1]
			if x > 3.0 and x < 4.0 and y > 3.0 and y < 4.0:
				return 1.0
			return 0.0
		return 0.0

	problem = {}
	problem["id"] = "blob2d"
	problem["\\sigma_t"] = lambda pWS: sigma_a(pWS) + sigma_s(pWS)
	problem["\\sigma_a"] = sigma_a
	problem["\\sigma_s"] = sigma_s
	problem["f_p"] = phase_shcoeffs
	problem["q"] = source_shcoeffs

	problem["order"] = 1
	problem["domain"] = util.Domain2D(7.0, 70)
	problem["staggered"] = True


	return problem



def checkerboard2d():
	# define problem ---
	def source( pWS ):
		x = pWS[0]
		y = pWS[1]
		if x > 3.0 and x < 4.0 and y > 3.0 and y < 4.0:
			return 1.0
		return 0.0

	def sigma_a( pWS ):
		x = pWS[0]
		y = pWS[1]
		cx = np.ceil(x)
		cy = np.ceil(y)
		g = 0
		if np.ceil((x+y)/2.0)*2.0 == (cx+cy) and cx > 1.0 and cx < 7.0 and cy > 1.0 and cy-2.0*np.abs(cx-4.0) < 4:
			g = 1
		return (1.0-g)*0 + g*10

	def sigma_s( pWS ):
		x = pWS[0]
		y = pWS[1]
		cx = np.ceil(x)
		cy = np.ceil(y)
		g = 0
		if np.ceil((x+y)/2.0)*2.0 == (cx+cy) and cx > 1.0 and cx < 7.0 and cy > 1.0 and cy-2.0*np.abs(cx-4.0) < 4:
			g = 1
		return (1.0-g)*1 + g*0


	def phase_shcoeffs( l, m, pWS ):
		if l == 0:
			return 1.0
		return 0.0

	def source_shcoeffs( l, m, pWS ):
		if l==0:
			x = pWS[0]
			y = pWS[1]
			if x > 3.0 and x < 4.0 and y > 3.0 and y < 4.0:
				return 1.0
			return 0.0
		return 0.0

	problem = {}
	problem["id"] = "checkerboard"
	problem["\\sigma_t"] = lambda pWS: sigma_a(pWS) + sigma_s(pWS)
	problem["\\sigma_a"] = sigma_a
	problem["\\sigma_s"] = sigma_s
	problem["f_p"] = phase_shcoeffs
	problem["q"] = source_shcoeffs

	problem["order"] = 1
	problem["domain"] = util.Domain2D(7.0, 70)
	problem["staggered"] = True


	return problem


def blurred( problem, stddev = 1.0 ):

	domain = problem["domain"]

	problem["id"] = "{}_blur{}".format(problem["id"], stddev)

	sigma_t_grid = domain.rasterize(problem["\\sigma_t"])
	albedo_grid = domain.rasterize(problem["\\sigma_s"])/sigma_t_grid

	problem["\\sigma_a"] = VoxelGrid(domain, domain.rasterize(problem["\\sigma_a"]))
	problem["\\sigma_a"].blur(stddev)
	problem["\\sigma_s"] = VoxelGrid(domain, domain.rasterize(problem["\\sigma_s"]))
	problem["\\sigma_s"].blur(stddev)
	problem["\\sigma_t"] = VoxelGrid(domain, problem["\\sigma_a"].voxels+problem["\\sigma_s"].voxels)

	return problem




def diffusion():
	def point_source(pWS):
		if np.linalg.norm(pWS - np.array([3.5, 3.5])) < 0.2:
			return 1.0
		return 0.0
	#pnb = pnbuilder.PNBuilder(order, domain)
	#pnb.add_terms(diffusion_terms())

	problem = {}
	problem["id"] = "diffusion"
	problem["\\sigma_t"] = lambda pWS: 10.0
	problem["q"] = point_source

	problem["order"] = 0
	problem["domain"] = util.Domain2D(7.0, 70)

	return problem









def debug():
	order = 1

	#sigma_t = Gradient(np.array([1.0, 1.0]))
	sigma_a = Constant(1.0)
	sigma_s = Gradient(np.array([1.0, 1.0]))
	sigma_t = Sum( sigma_a, sigma_s )


	problem = {}
	problem["id"] = "debug_term"
	problem["\\sigma_t"] = lambda pWS: sigma_t(pWS)
	problem["\\sigma_s"] = lambda pWS: sigma_s(pWS)
	problem["sigma_t_instance"] = sigma_t
	problem["sigma_s_instance"] = sigma_s
	q_instance = SHEXP(order, [Constant(1.0), Constant(0.0), Constant(0.0), Constant(0.0)])
	def q_fun( l, m, pWS ):
		shindex = shtools.shIndex(l,m)
		if shindex is None or l > order or m < -l or m > l:
			return 0.0
		return q_instance.coeff_functions[shindex](pWS)
	problem["q"] = q_fun
	problem["q_instance"] = q_instance


	f_instance = SHEXP(order, [Constant(1.0), Constant(0.0), Constant(0.0), Constant(0.0)])
	def f_fun( l, m, pWS ):
		shindex = shtools.shIndex(l,m)
		if shindex is None or l > order or m < -l or m > l:
			return 0.0
		return f_instance.coeff_functions[shindex](pWS)
	problem["f_p"] = f_fun
	problem["f_instance"] = f_instance


	#functions["L"] = SHEXP(order, [Constant(1.0), Constant(0.2), Constant(0.3), Constant(0.4)])
	problem["L"] = SHEXP(order, [RBF(1.0), RBF(0.2), RBF(0.3), RBF(0.4)])


	problem["staggered"] = True
	problem["order"] = order
	problem["domain"] = util.Domain2D(0.1, 5)

	return problem