# This module contains different classes for working with field variables.
# These are used to describe the RTE parameters and also for working with
# PN solutions.

import numpy as np
import util
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
    #def __init__(self, order, coeff_functions):
    #    self.order = order
    #    self.coeff_functions = coeff_functions
    def __init__(self, order):
    	self.order = order
    	self.coeff_functions = [Constant(0.0) for i in range(util.num_sh_coeffs(order))]
    def set_coeff(self, index, field):
    	self.coeff_functions[index] = field
    def __call__(self, pWS, omega):
        (theta, phi) = util.sphericalCoordinates(omega)
        coeffs = [f(pWS) for f in self.coeff_functions]
        return util.sh_sum(self.order, coeffs, theta, phi)
    def dx(self, pWS, omega):
    	(theta, phi) = util.sphericalCoordinates(omega)
    	coeffs_dx = [f.dx(pWS) for f in self.coeff_functions]
    	return util.sh_sum(self.order, coeffs_dx, theta, phi)
    def dxdx(self, pWS, omega):
    	(theta, phi) = util.sphericalCoordinates(omega)
    	coeffs_dxdx = [f.dxdx(pWS) for f in self.coeff_functions]
    	return util.sh_sum(self.order, coeffs_dxdx, theta, phi)
    def dxdy(self, pWS, omega):
    	(theta, phi) = util.sphericalCoordinates(omega)
    	coeffs_dxdy = [f.dxdy(pWS) for f in self.coeff_functions]
    	return util.sh_sum(self.order, coeffs_dxdy, theta, phi)
    def dydx(self, pWS, omega):
    	return self.dxdy(pWS, omega)
    def dy(self, pWS, omega):
    	(theta, phi) = util.sphericalCoordinates(omega)
    	coeffs_dy = [f.dy(pWS) for f in self.coeff_functions]
    	return util.sh_sum(self.order, coeffs_dy, theta, phi)
    def dydy(self, pWS, omega):
    	(theta, phi) = util.sphericalCoordinates(omega)
    	coeffs_dydy = [f.dydy(pWS) for f in self.coeff_functions]
    	return util.sh_sum(self.order, coeffs_dydy, theta, phi)
    def dz(self, pWS, omega):
    	(theta, phi) = util.sphericalCoordinates(omega)
    	coeffs_dz = [f.dz(pWS) for f in self.coeff_functions]
    	return util.sh_sum(self.order, coeffs_dz, theta, phi)
    def dzdz(self, pWS, omega):
    	(theta, phi) = util.sphericalCoordinates(omega)
    	coeffs_dzdz = [f.dzdz(pWS) for f in self.coeff_functions]
    	return util.sh_sum(self.order, coeffs_dzdz, theta, phi)
    def integral_over_solid_angle(self, pWS):
    	return self.coeff_functions[0](pWS)*np.sqrt(4*np.pi)



class VoxelGrid(object):
	def __init__( self, voxels, domain, offset = np.array([0.5, 0.5]) ):
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
		res = self.domain.resolution()
		if v0[0] < 0:
			v0[0] = 0
		elif v0[0] >= res[0]:
			v0[0] = res[0]-1
		if v0[1] < 0:
			v0[1] = 0
		elif v0[1] >= res[1]:
			v0[1] = res[1]-1
		if v1[0] < 0:
			v1[0] = 0
		elif v1[0] >= res[0]:
			v1[0] = res[0]-1
		if v1[1] < 0:
			v1[1] = 0
		elif v1[1] >= res[1]:
			v1[1] = res[1]-1

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
	def dx(self, pWS):
		step = np.array([0.5,0.0])*self.domain.voxelSize()
		a = self(pWS-step)
		b = self(pWS+step)
		return (b-a)/self.domain.voxelSize()[0]
	def dxdx(self, pWS):
		step = np.array([0.5,0.0])*self.domain.voxelSize()
		a = self.dx(pWS-step)
		b = self.dx(pWS+step)
		return (b-a)/self.domain.voxelSize()[0]
	def dxdy(self, pWS):
		step = np.array([0.0,0.5])*self.domain.voxelSize()
		a = self.dx(pWS-step)
		b = self.dx(pWS+step)
		return (b-a)/self.domain.voxelSize()[1]
	def dy(self, pWS):
		step = np.array([0.0,0.5])*self.domain.voxelSize()
		a = self(pWS-step)
		b = self(pWS+step)
		return (b-a)/self.domain.voxelSize()[1]
	def dydy(self, pWS):
		step = np.array([0.0,0.5])*self.domain.voxelSize()
		a = self.dy(pWS-step)
		b = self.dy(pWS+step)
		return (b-a)/self.domain.voxelSize()[1]
	def dydx(self, pWS):
		return self.dxdy(pWS)
	def dz(self, pWS):
		return 0.0
	#def in_bounds(self, voxel):
	#	if voxel[0] < 0 or voxel[0]>=self.domain.res_x or voxel[1] < 0 or voxel[1]>=self.domain.res_y:
	#		return False
	#	return True




def from_pn_solution( x_complex, pnb ):
	'''Produces a radiance field from given PN solution.'''


	# now for each coefficient
	shexp = SHEXP(pnb.order())
	#shexp = pnbuilder.SHEXP(pnb.order())
	for i in range(pnb.num_coeffs()):
		(l,m) = pnb.lm_index(i)
		offset = pnb.unknown_offset(i)

		res = pnb.domain_cpp.resolution()
		voxels = np.zeros( (res[0], res[1]), dtype=complex )
		for voxel_i in range(res[1]):
			for voxel_j in range(res[0]):
				voxels[voxel_i, voxel_j] = x_complex[pnb.global_index(voxel_i, voxel_j, i)]
		grid = VoxelGrid( voxels, pnb.domain_cpp, offset*0.5 )
		#grid = pnbuilder.VoxelGrid( voxels, pnb.domain_cpp, offset*0.5 )

		shexp.set_coeff(i,grid)

	return shexp
