import numpy as np
from scipy.special import sph_harm, lpmv
from sympy.physics.quantum.cg import CG
import matplotlib.pyplot as plt



# utility functions =========================================
def numSHCoeffs(order):
    return (order + 1) * (order + 1)

def shIndex( l, m):
    return l * (l + 1) + m;

def sh_sum( order, coeffs, theta, phi ):
    result = 0.0
    for l in range(order+1):
        for m in range(-l, l+1):
            result+=coeffs[shIndex(l,m)]*sph_harm(m, l, phi, theta)
    return result

def project_sh( fun, order ):
    coeffs = np.zeros(numSHCoeffs(order), dtype=complex)
    for l in range(order+1):
        for m in range(-l, l+1):
            # NB: we use the complex conjugate of the spherical harmonics function for projection
            # see https://en.wikipedia.org/wiki/Spherical_harmonics#Spherical_harmonics_expansion
            coeffs[shIndex(l,m)] = integrate_sphere( lambda theta, phi: fun(theta, phi)*np.conj(sph_harm(m, l, phi, theta)) )
    return coeffs            

def integrate_sphere( fun ):
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
    
# here we build a table of precomputed Clebsch-Gordan coefficients
CG_table_order = 7
print("building CG table with order={}".format(CG_table_order))
CG_table = np.zeros( (numSHCoeffs(CG_table_order), numSHCoeffs(CG_table_order), numSHCoeffs(CG_table_order)) )

for l in range(CG_table_order+1):
    for m in range(-l, l+1):
        for l1 in range(CG_table_order+1):
            for m1 in range(-l1, l1+1):
                for l2 in range(CG_table_order+1):
                    for m2 in range(-l2, l2+1):
                        CG_table[shIndex(l, m), shIndex(l1, m1), shIndex(l2, m2)] = CG(l, m, l1, m1, l2, m2).doit().evalf()











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




def sphericalCoordinates( d ):
	theta = np.acos(d[2])
	phi = np.atan2(d[1], d[0])
	if phi < 0.0:
		phi += 2.0*np.pi
	return (theta, phi)










def shIndex( l, m):
    return l * (l + 1) + m;

class PhaseFunctionSH(object):
	def __init__(self, order, P):
		self.P = P
		self.order = order
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
