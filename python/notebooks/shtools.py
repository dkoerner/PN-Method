import numpy as np
import math
from scipy.special import sph_harm, lpmv
from sympy.physics.quantum.cg import CG
import matplotlib.pyplot as plt
from multiprocessing import Pool



# utility functions =========================================
def numSHCoeffs(order):
    return (order + 1) * (order + 1)

def shIndex( l, m):
	if l<0 or m < -l or m > l:
		return None
	return l * (l + 1) + m

lm_index_table = None
lm_index_table_order = 7
def lmIndex(i):
	global lm_index_table

	if lm_index_table == None:
		global lm_index_table_order
		lm_index_table = [0 for i in range(numSHCoeffs(lm_index_table_order+1))]
		for l in range(lm_index_table_order+1):
			for m in range(-l, l+1):
				lm_index_table[shIndex(l, m)] = (l,m)


	return lm_index_table[i]


def sh_sum( order, coeffs, theta, phi ):
    result = 0.0
    for l in range(order+1):
        for m in range(-l, l+1):
            result+=coeffs[shIndex(l,m)]*sph_harm(m, l, phi, theta)
    return result

def project_sh( fun, order, debug=False ):
    coeffs = np.zeros(numSHCoeffs(order), dtype=complex)
    for l in range(order+1):
        for m in range(-l, l+1):
            # NB: we use the complex conjugate of the spherical harmonics function for projection
            # see https://en.wikipedia.org/wiki/Spherical_harmonics#Spherical_harmonics_expansion
            coeffs[shIndex(l,m)] = integrate_sphere( lambda theta, phi: fun(theta, phi)*np.conj(sph_harm(m, l, phi, theta)), debug )
    return coeffs            

def project_sh_coeff( fun, l, m ):
	return integrate_sphere( lambda theta, phi: fun(theta, phi)*np.conj(sph_harm(m, l, phi, theta)) )


def test(x):
	print(x)
	return x
project_sh_coeff_x_data = None
def project_sh_coeff_x( fun, l, m ):
	global project_sh_coeff_x_data
	if project_sh_coeff_x_data is None:
		min_theta = 0
		max_theta = np.pi
		min_phi = 0
		max_phi = 2.0*np.pi

		resolution_theta=128 # resolution along theta angle
		resolution_phi=256 # resolution along phi angle
		dtheta = (max_theta-min_theta)/float(resolution_theta)
		dphi = (max_phi-min_phi)/float(resolution_phi)

		pixel_area = dtheta*dphi

		theta_phi = np.zeros( (resolution_theta, resolution_phi, 2) )

		ins = []
		for t in range(resolution_theta):
			for p in range(resolution_phi):
				phi = dphi*p
				theta = dtheta*t
				theta_phi[t, p, 0] = theta
				theta_phi[t, p, 1] = phi
				ins.append((theta, phi))

		project_sh_coeff_x_data = {}
		project_sh_coeff_x_data["theta_phi"] = theta_phi
		project_sh_coeff_x_data["pixel_area"] = pixel_area
		project_sh_coeff_x_data["pool"] = Pool(5)
		project_sh_coeff_x_data["ins"] = ins
		

	theta_phi = project_sh_coeff_x_data["theta_phi"]

	if not (l,m) in project_sh_coeff_x_data:
		pixel_area = project_sh_coeff_x_data["pixel_area"]
		factor = np.zeros( (theta_phi.shape[0],theta_phi.shape[1]) , dtype=complex )
		for t in range(theta_phi.shape[0]):
			for p in range(theta_phi.shape[1]):
				theta = theta_phi[t, p, 0]
				phi = theta_phi[t, p, 1]
				factor[t, p] = np.conj(sph_harm(m, l, phi, theta))*pixel_area*np.sin(theta)
		project_sh_coeff_x_data[(l,m)] = factor.reshape(theta_phi.shape[0]*theta_phi.shape[1])

	factor = project_sh_coeff_x_data[(l,m)]


	res = project_sh_coeff_x_data["pool"].map( fun, project_sh_coeff_x_data["ins"] )
	result = np.sum(np.array(res)*factor)
	'''
	# integrate sphere
	for t in range(theta_phi.shape[0]):
		for p in range(theta_phi.shape[1]):
			theta = theta_phi[t, p, 0]
			phi = theta_phi[t, p, 1]
			result+= fun(theta, phi)*factor[t,p]
	'''

	return result


# condon-shortley phase
def csp( m ):
	if m % 2 == 0:
		return 1.0
	return -1.0


# these functions are the coefficients for the recursive relation of the sh basis function
# (see p. 4 in the starmap paper)
def a_lm( l, m ):
	base = (l-m+1)*(l+m+1)/((2*l+3)*(2*l+1))
	if base < 0.0:
		return np.nan
	return np.sqrt(base)
def b_lm( l, m ):
	base = (l-m)*(l+m)/((2*l+1)*(2*l-1))
	if base < 0.0:
		return np.nan
	return np.sqrt(base)
def c_lm( l, m ):
	base = (l+m+1)*(l+m+2)/((2*l+3)*(2*l+1))
	if base < 0.0:
		return np.nan
	return np.sqrt(base)
def d_lm( l, m ):
	base = (l-m)*(l-m-1)/((2*l+1)*(2*l-1))
	if base < 0.0:
		return np.nan
	return np.sqrt(base)
def e_lm( l, m ):
	base = (l-m+1)*(l-m+2)/((2*l+3)*(2*l+1))
	if base < 0.0:
		return np.nan
	return np.sqrt(base)
def f_lm( l, m ):
	base = (l+m)*(l+m-1)/((2*l+1)*(2*l-1))
	if base < 0.0:
		return np.nan
	return np.sqrt(base)

    
def integrate_sphere( fun, debug = False ):
    result = 0.0

    min_theta = 0
    max_theta = np.pi
    min_phi = 0
    max_phi = 2.0*np.pi

    resolution_theta=128 # resolution along theta angle
    resolution_phi=256 # resolution along phi angle
    dtheta = (max_theta-min_theta)/float(resolution_theta)
    dphi = (max_phi-min_phi)/float(resolution_phi)

    pixel_area = dtheta*dphi

    for t in range(resolution_theta):
        for p in range(resolution_phi):
            phi = dphi*p
            theta = dtheta*t
            result+=fun(theta, phi)*pixel_area*np.sin(theta)

    return result

def plot_spherical_function( fun, vmin=-0.5, vmax=0.5 ):
    # plot sh functions
    theta = np.arange(0.0, 1.0, 0.01)*np.pi
    phi = np.arange(0.0, 1.0, 0.01)*2.0*np.pi

    f_img = np.zeros((theta.shape[0], phi.shape[0]))
    for j in range(phi.shape[0]):
        for i in range(theta.shape[0]):
            f_img[i, j] = fun(theta[i], phi[j])
    plt.imshow(f_img, interpolation="nearest", cmap='jet', vmin=vmin, vmax=vmax)
    


def sphericalCoordinates( d ):
	theta = np.arccos(d[2])
	# np.atan2 returns wrong values
	phi = math.atan2(d[1], d[0])
	if phi < 0.0:
		phi += 2.0*np.pi
	return (theta, phi)

def sphericalDirection( theta, phi ):
	sinTheta = np.sin(theta)
	sinPhi = np.sin(phi)
	cosTheta = np.cos(theta)
	cosPhi = np.cos(phi)
	return np.array([sinTheta * cosPhi, sinTheta * sinPhi, cosTheta])


def coordinateSystem( n ):
	'''returns a matrix from given normal'''
	if np.abs(n[0]) > np.abs(n[1]):
		invLen = 1.0/np.sqrt(n[0]*n[0] + n[2]*n[2])
		c = np.array([n[2]*invLen, 0.0, -n[0]*invLen])
	else:
		invLen = 1.0/np.sqrt(n[1]*n[1] + n[2]*n[2])
		c = np.array([0.0, n[2]*invLen, -n[1]*invLen])
	b = np.cross(c, n)
	return np.array( [[b[i], c[i], n[i]] for i in range(3)] )





def DoubleFactorial(x):
    s = 1.0;
    n = x;
    while n > 1.0:
        s *= n
        n -= 2.0
    return s


def P_lm(l, m, x):
    # Compute Pmm(x) = (-1)^m(2m - 1)!!(1 - x^2)^(m/2), where !! is the double
    # factorial.
    pmm = 1.0
    # P00 is defined as 1.0, do don't evaluate Pmm unless we know m > 0
    if m > 0:
        if m % 2 == 0:
            sign = 1
        else:
            sign = -1
        pmm = sign * DoubleFactorial(2 * m - 1) * np.power(1 - x * x, m / 2.0)

    if l == m:
        #Pml is the same as Pmm so there's no lifting to higher bands needed
        return pmm
    # Compute Pmm+1(x) = x(2m + 1)Pmm(x)
    pmm1 = x * (2 * m + 1) * pmm
    if  l == m + 1:
        # Pml is the same as Pmm+1 so we are done as well
        return pmm1;
    # Use the last two computed bands to lift up to the next band until l is
    # reached, using the recurrence relationship:
    # Pml(x) = (x(2l - 1)Pml-1 - (l + m - 1)Pml-2) / (l - m)
    #for (int n = m + 2; n <= l; n++):
    for n in range(m+2, l+1):
        pmn = (x * (2 * n - 1) * pmm1 - (n + m - 1) * pmm) / (n - m)
        pmm = pmm1
        pmm1 = pmn
    # Pmm1 at the end of the above loop is equal to Pml
    return pmm1

def Y_lm( l, m, theta, phi ):
    if m < 0:
        if m % 2 == 0:
            sign = 1.0
        else:
            sign = -1.0
        return np.conj(sign*Y_lm(l, -m, theta, phi))
    return np.sqrt( (2.0*l+1)/(4.0*np.pi)*np.math.factorial(l-m)/np.math.factorial(l+m) )*np.complex(np.cos(m*phi), np.sin(m*phi))*P_lm(l, m, np.cos(theta))



class PhaseFunctionSH(object):
	def __init__(self, filename):
		file = open(filename, "r")
		self.order = int(file.readline().rstrip())
		self.numCoeffs = (self.order + 1) * (self.order + 1)
		self.P = np.zeros((self.numCoeffs, self.numCoeffs), dtype=complex)
		for i in range(self.numCoeffs):
			for j in range(self.numCoeffs):
				t = [float(v) for v in file.readline().rstrip().split()]
				self.P[i, j] = np.complex(t[0], t[1])
		file.close()
	def __call__(self, wi, wo):
		(theta_i, phi_i) = sphericalCoordinates(wi)
		(theta_o, phi_o) = sphericalCoordinates(wo)
		result = np.complex(0.0, 0.0)
		for l in range(self.order+1):
			for m in range(-l, l+1):
				for l2 in range(self.order+1):
					for m2 in range(-l2, l2+1):
						i = shIndex(l, m)
						j = shIndex(l2, m2)
						result += self.P[i, j]*Y_lm( l, m, theta_i, phi_i )*Y_lm( l2, m2, theta_o, phi_o )
		return result
        
class SGGX3D(object):
	def __init__(self, m, projectedAreas):
		S_11 = projectedAreas[0]*projectedAreas[0];
		S_22 = projectedAreas[1]*projectedAreas[1];
		S_33 = projectedAreas[2]*projectedAreas[2];
		self.S = m.dot(np.diag([S_11, S_22, S_33])).dot(m.T)
		self.S_xx = self.S[0,0]
		self.S_xy = self.S[0,1]
		self.S_xz = self.S[0,2]
		self.S_yy = self.S[1,1]
		self.S_yz = self.S[1,2]
		self.S_zz = self.S[2,2]

        
	def D(self, wm):
		'''returns distribution of normals (Eq.12 in the paper)'''
		detS = self.S_xx*self.S_yy*self.S_zz - self.S_xx*self.S_yz*self.S_yz - self.S_yy*self.S_xz*self.S_xz - self.S_zz*self.S_xy*self.S_xy + 2.0*self.S_xy*self.S_xz*self.S_yz
		den = wm[0]*wm[0]*(self.S_yy*self.S_zz - self.S_yz*self.S_yz) + wm[1]*wm[1]*(self.S_xx*self.S_zz - self.S_xz*self.S_xz) + wm[2]*wm[2]*(self.S_xx*self.S_yy - self.S_xy*self.S_xy)+2.0*(wm[0]*wm[1]*(self.S_xz*self.S_yz - self.S_zz*self.S_xy) + wm[0]*wm[2]*(self.S_xy*self.S_yz - self.S_yy*self.S_xz) + wm[1]*wm[2]*(self.S_xy*self.S_xz - self.S_xx*self.S_yz))
		D = np.power(np.abs(detS), 1.5) / (np.pi*den*den)
		return D;

	def __call__( self, d ):
		sigma_squared = d[0]*d[0]*self.S_xx + d[1]*d[1]*self.S_yy + d[2]*d[2]*self.S_zz + 2.0 * (d[0]*d[1]*self.S_xy + d[0]*d[2]*self.S_xz + d[1]*d[2]*self.S_yz)
		# conditional to avoid numerical errors
		if sigma_squared > 0.0:
			return np.sqrt(sigma_squared)
		return 0.0
        
        
class PhaseFunctionSGGX(object):
	def __init__(self, sggx):
		self.sggx = sggx
	def __call__(self, wi, wo):
		wh = wi+wo
		wh = wh/np.linalg.norm(wh)
		return 0.25*self.sggx.D(wh)/self.sggx(wi)






# here we build a table of precomputed Clebsch-Gordan coefficients
CG_table_order = 7
CG_table = None
def init_shtools( max_order=7, build_CG_table = False ):
	if build_CG_table == True:
		global CG_table_order
		global CG_table
		CG_table_order = max_order
		print("building CG table with order={}".format(CG_table_order))
		CG_table = np.zeros( (numSHCoeffs(CG_table_order), numSHCoeffs(CG_table_order), numSHCoeffs(CG_table_order)) )

		for l in range(CG_table_order+1):
			for m in range(-l, l+1):
				for l1 in range(CG_table_order+1):
					for m1 in range(-l1, l1+1):
						for l2 in range(CG_table_order+1):
							for m2 in range(-l2, l2+1):
								CG_table[shIndex(l, m), shIndex(l1, m1), shIndex(l2, m2)] = CG(l, m, l1, m1, l2, m2).doit().evalf()


















