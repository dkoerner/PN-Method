# This module contains different classes for working with field variables.
# These are used to describe the RTE parameters and also for working with
# PN solutions.

import numpy as np
import util
from scipy.ndimage.filters import gaussian_filter
import itertools
import pnbuilder




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
    def __call__(self, location, omega):
        (theta, phi) = util.sphericalCoordinates(omega)
        coeffs = [f(location) for f in self.coeff_functions]
        return util.sh_sum(self.order, coeffs, theta, phi)
    def dx(self, location, omega):
    	(theta, phi) = util.sphericalCoordinates(omega)
    	coeffs_dx = [f.dx(location) for f in self.coeff_functions]
    	return util.sh_sum(self.order, coeffs_dx, theta, phi)
    def dxdx(self, location, omega):
    	(theta, phi) = util.sphericalCoordinates(omega)
    	coeffs_dxdx = [f.dxdx(location) for f in self.coeff_functions]
    	return util.sh_sum(self.order, coeffs_dxdx, theta, phi)
    def dxdy(self, location, omega):
    	(theta, phi) = util.sphericalCoordinates(omega)
    	coeffs_dxdy = [f.dxdy(location) for f in self.coeff_functions]
    	return util.sh_sum(self.order, coeffs_dxdy, theta, phi)
    def dydx(self, location, omega):
    	return self.dxdy(location, omega)
    def dy(self, location, omega):
    	(theta, phi) = util.sphericalCoordinates(omega)
    	coeffs_dy = [f.dy(location) for f in self.coeff_functions]
    	return util.sh_sum(self.order, coeffs_dy, theta, phi)
    def dydy(self, location, omega):
    	(theta, phi) = util.sphericalCoordinates(omega)
    	coeffs_dydy = [f.dydy(location) for f in self.coeff_functions]
    	return util.sh_sum(self.order, coeffs_dydy, theta, phi)
    def dz(self, location, omega):
    	(theta, phi) = util.sphericalCoordinates(omega)
    	coeffs_dz = [f.dz(location) for f in self.coeff_functions]
    	return util.sh_sum(self.order, coeffs_dz, theta, phi)
    def dzdz(self, location, omega):
    	(theta, phi) = util.sphericalCoordinates(omega)
    	coeffs_dzdz = [f.dzdz(location) for f in self.coeff_functions]
    	return util.sh_sum(self.order, coeffs_dzdz, theta, phi)
    def integral_over_solid_angle(self, location):
    	return self.coeff_functions[0](location)*np.sqrt(4*np.pi)

class VoxelGrid(object):
	def __init__( self, domain, voxels, offset = np.array([0.5, 0.5]) ):
		self.domain = domain
		self.voxels = voxels
		self.offset = offset



	def sample( self, voxel_i, voxel_j ):
		return self.voxels[voxel_i, voxel_j]
		

	def eval( self, pVS ):
		v0 = np.floor(pVS).astype(int)
		v1 = v0+np.array([1, 1])
		t = pVS - v0

		# clamp to domain boundaries
		if v0[0] < 0:
			v0[0] = 0
		elif v0[0] >= self.domain.res_x:
			v0[0] = self.domain.res_x-1
		if v0[1] < 0:
			v0[1] = 0
		elif v0[1] >= self.domain.res_y:
			v0[1] = self.domain.res_y-1
		if v1[0] < 0:
			v1[0] = 0
		elif v1[0] >= self.domain.res_x:
			v1[0] = self.domain.res_x-1
		if v1[1] < 0:
			v1[1] = 0
		elif v1[1] >= self.domain.res_y:
			v1[1] = self.domain.res_y-1

		#return self.sample(v0[0], v0[1])

		result = 0.0
		result += self.sample(v0[0], v0[1])*(1.0-t[0])*(1.0-t[1])
		result += self.sample(v1[0], v0[1])*t[0]*(1.0-t[1])
		result += self.sample(v0[0], v1[1])*(1.0-t[0])*t[1]
		result += self.sample(v1[0], v1[1])*t[0]*t[1]
		return result


		
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
					voxels[voxel_i, voxel_j] = (self.sample(c[0], c[1]) - self.sample(m[0], m[1]))/self.domain.voxelsize[dimension]
				elif not self.in_bounds(m):
					# backward differences
					voxels[voxel_i, voxel_j] = (self.sample(p[0], p[1]) - self.sample(c[0], c[1]))/self.domain.voxelsize[dimension]
				else:
					# finite differences
					voxels[voxel_i, voxel_j] = (self.sample(p[0], p[1]) - self.sample(m[0], m[1]))/(2.0*self.domain.voxelsize[dimension])
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



class CoefficientGrid(object):
	def __init__(self, domain, stride, index, offset, x):
		self.domain = domain
		self.stride = stride
		self.offset = offset

		voxels = np.zeros( (domain.res_x, domain.res_y), dtype=complex )
		for voxel_i in range(domain.res_y):
			for voxel_j in range(domain.res_x):
				voxels[voxel_i, voxel_j] = x[self.get_global_index(voxel_i, voxel_j, index)]
		self.voxels = voxels

		#super().__init__(domain, voxels, offset)

		#self.build_derivatives()


	def get_global_index( self, voxel_i, voxel_j, coeff_index ):
		v = voxel_j*self.domain.res_x + voxel_i
		return v*self.stride + coeff_index

	def sample( self, voxel ):
		return self.voxels[voxel[0], voxel[1]]

	def __call__(self, location):
		numDimensions = 2
		# check if the location, at which to evaluate the unknown,
		# matches the actual grid location of the unknown
		# this is true for first order equation with non-anisotropic media
		location_offset = location.getOffset()
		unknown_offset = self.offset
		if (location_offset == unknown_offset).all():
			# unknown location and eval location are the same spot
			# no interpolation needed
			return self.sample(location.getVoxel())
		elif location_offset[0] == unknown_offset[0]:
			# in the previous if-clause, we checked for location and unknown to be exactly equal
			# now if their offset matches only in one dimension, then we can simply interpolate
			# between the two neighbouring datapoints in that dimension

			# TODO: generalize to 3D

			l = self.sample(location.getShiftedLocation(np.array([0, 1])).getVoxel())
			r = self.sample(location.getShiftedLocation(np.array([0, -1])).getVoxel())
			return 0.5*(l+r)
		elif location_offset[1] == unknown_offset[1]:
			l = self.sample(location.getShiftedLocation(np.array([1, 0])).getVoxel())
			r = self.sample(location.getShiftedLocation(np.array([-1, 0])).getVoxel())
			return 0.5*(l+r)
		#elif (location_offset[0]+1)%2 == unknown_offset[0] and (location_offset[1]+1)%2 == unknown_offset[1]:
		else:
			# the offsets of evaluation and unknown do not match in any dimension, therefore 
			# we can conclude that the location of the unknown is on the diagonals to the eval location

			# the unknown is located diagonal to the position at which we want to evaluate it
			# therefore we will do an equal weighted sum of all for surrounding unknowns
			offset_combinations = itertools.product(*[[-1, 1] for d in range(numDimensions)])
			num_offset_combinations = 2**numDimensions

			result = 0.0
			for o in offset_combinations:
				result += self.sample(location.getShiftedLocation(np.array(o)).getVoxel())
			return result/num_offset_combinations




	def dx(self, location):
		step = np.array([1,0], dtype=int)
		a = self(location.getShiftedLocation(-step))
		b = self(location.getShiftedLocation(step))
		return (b-a)/self.domain.voxelsize[0]
	def dxdx(self, location):
		step = np.array([1,0], dtype=int)
		a = self.dx(location.getShiftedLocation(-step))
		b = self.dx(location.getShiftedLocation(step))
		return (b-a)/self.domain.voxelsize[0]
	def dxdy(self, location):
		step = np.array([0,1], dtype=int)
		a = self.dx(location.getShiftedLocation(-step))
		b = self.dx(location.getShiftedLocation(step))
		return (b-a)/self.domain.voxelsize[1]
	def dy(self, location):
		step = np.array([0,1], dtype=int)
		a = self(location.getShiftedLocation(-step))
		b = self(location.getShiftedLocation(step))
		return (b-a)/self.domain.voxelsize[1]
	def dydy(self, location):
		step = np.array([0,1], dtype=int)
		a = self.dy(location.getShiftedLocation(-step))
		b = self.dy(location.getShiftedLocation(step))
		return (b-a)/self.domain.voxelsize[1]
	def dydx(self, location):
		return self.dxdy(location)
	def dz(self, location):
		return 0.0
	def test(self):
		for s in range(10):
			for t in range(10):
				print(self.voxels[s,t])

def from_pn_solution( x_complex, pnb ):
	'''Produces a radiance field from given PN solution.'''

	# initialize all coefficients with a zero field
	coefficient_fields = [Constant(0.0) for i in range(util.num_sh_coeffs(pnb.order()))]

	# now for each coefficient
	for i in range(pnb.num_coeffs()):
		(l,m) = pnb.lm_index(i)
		offset = pnb.unknown_offset(i)
		grid = CoefficientGrid(pnb.domain, pnb.num_coeffs(), i, offset, x_complex )

		'''
		res = pnb.domain_cpp.resolution()
		voxels = np.zeros( (res[0], res[1]), dtype=complex )
		for voxel_i in range(res[1]):
			for voxel_j in range(res[0]):
				voxels[voxel_i, voxel_j] = x_complex[pnb.global_index(voxel_i, voxel_j, i)]
		grid = pnbuilder.CoefficientGrid( voxels, pnb.domain_cpp, offset )
		'''

		coefficient_fields[util.sh_index(l,m)] = grid

	#coefficient_fields[0].test()
	return SHEXP( pnb.order(), coefficient_fields )
