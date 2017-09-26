# This module contains various helper and utility functions such as 
# related to spherical harmonics, I/O etc.

import numpy as np
import scipy.io
import math

from multiprocessing import Pool
from scipy.special import sph_harm, lpmv


# Spherical Harmonics related =====================================================

def num_sh_coeffs(order):
    return (order + 1) * (order + 1)

def sh_index( l, m):
	if l<0 or m < -l or m > l:
		return None
	return l * (l + 1) + m

def sh_sum( order, coeffs, theta, phi ):
    result = 0.0
    for l in range(order+1):
        for m in range(-l, l+1):
            result+=coeffs[sh_index(l,m)]*sph_harm(m, l, phi, theta)
    return result

#'''
project_sh_coeff_data = None
def project_sh_coeff( fun, l, m ):
	global project_sh_coeff_data
	if project_sh_coeff_data is None:
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

		project_sh_coeff_data = {}
		project_sh_coeff_data["theta_phi"] = theta_phi
		project_sh_coeff_data["pixel_area"] = pixel_area
		#project_sh_coeff_data["pool"] = Pool(5)
		project_sh_coeff_data["ins"] = ins
		

	theta_phi = project_sh_coeff_data["theta_phi"]

	if not (l,m) in project_sh_coeff_data:
		pixel_area = project_sh_coeff_data["pixel_area"]
		factor = np.zeros( (theta_phi.shape[0],theta_phi.shape[1]) , dtype=complex )
		for t in range(theta_phi.shape[0]):
			for p in range(theta_phi.shape[1]):
				theta = theta_phi[t, p, 0]
				phi = theta_phi[t, p, 1]
				factor[t, p] = np.conj(sph_harm(m, l, phi, theta))*pixel_area*np.sin(theta)
		project_sh_coeff_data[(l,m)] = factor.reshape(theta_phi.shape[0]*theta_phi.shape[1])

	factor = project_sh_coeff_data[(l,m)]


	#res = project_sh_coeff_data["pool"].map( fun, project_sh_coeff_data["ins"] )
	res = [fun(theta, phi) for (theta, phi) in project_sh_coeff_data["ins"]]
	


	result = np.sum(np.array(res)*factor)

	return result
#'''


'''
def project_sh_coeff( fun, l, m ):
	return integrate_sphere( lambda theta, phi: fun(theta, phi)*np.conj(sph_harm(m, l, phi, theta)) )

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
'''
# misc ==============================================================
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


def compare_matrices( A0, A1, id0, id1 ):
    if A0.shape != A1.shape:
        raise ValueError("unmatched shape {}={} {}={}".format(id0, str(A0.shape), id1, str(A1.shape)))

    if A0.dtype != A1.dtype:
        raise ValueError("unmatched dtype")

    diff = A0-A1
    abs_real_diff = np.abs(np.real(diff))
    abs_imag_diff = np.abs(np.imag(diff))

    max_index = np.unravel_index(np.argmax(abs_real_diff), abs_real_diff.shape)
    mm = abs_real_diff[max_index]
    print("real: max={} max element=({}, {})".format(mm,max_index[0], max_index[1]))
    max_index = np.unravel_index(np.argmax(abs_imag_diff), abs_imag_diff.shape)
    mm = abs_imag_diff[max_index]
    print("imag: max={} max element=({}, {})".format(mm,max_index[0], max_index[1]))

# I/O related =====================================================
def rasterize( fun, domain, offset = np.array([0.5, 0.5]), dtype=float ):
	res = domain.resolution()
	shape = (res[0], res[1])
	voxels = np.zeros(shape, dtype=dtype)
	for i in range(0, res[0]):
		for j in range(0, res[1]):
			pVS = np.array([i, j]) + offset
			pWS = domain.voxelToWorld(pVS)
			voxels[i, j] = fun(pWS)
	return voxels

def write_pn_system(filename, sys, problem, x=None):
	domain = problem["domain"]
	A = sys.get_A_real()
	b = sys.get_b_real()
	#A = sys.get_boundary_A_real()
	#b = sys.get_boundary_b_real()

	info = {}
	info["order"] = sys.getOrder()
	info["numCoeffs"] = sys.getNumCoefficients()
	info["resolution"] = domain.resolution()

	data = {}

	data["id"] = problem["id"]

	data["info"] = info
	if not A is None:
		data['A'] = A
	if not b is None:
		data['b'] = b
	if not x is None:
		data['x'] = x

	data['sigma_s'] = rasterize(lambda pWS: np.real(problem["sigma_s"](pWS)), domain)
	data['sigma_a'] = rasterize(lambda pWS: np.real(problem["sigma_a"](pWS)), domain)
	data['sigma_t'] = rasterize(lambda pWS: np.real(problem["sigma_t"](pWS)), domain)
	data['q00']     = rasterize(lambda pWS: np.real(problem['q'][0](pWS)), domain)

	print("writing PN system to {}".format(filename))
	scipy.io.savemat(filename, data)


def load_pn_system( filename ):
    #print("loading PN solution from {}".format(filename))
    data = scipy.io.loadmat(filename)
   
    result = {}
    if "sigma_t" in data:
        result["sigma_t"] = data["sigma_t"]
    if "sigma_a" in data:
        result["sigma_a"] = data["sigma_a"]
    if "sigma_s" in data:
        result["sigma_s"] = data["sigma_s"]
    if "q00" in data:
        result["q00"] = data["q00"]
    if "x" in data:
        result["x"] = data["x"]
    if "b" in data:
        result["b"] = data["b"]
    if "A" in data:
        result["A"] = data["A"]
    if "info" in data:
        info = data["info"]
        result["order"] = info["order"][0][0][0][0]
        result["numCoeffs"] = info["numCoeffs"][0][0][0][0]
        result["resolution"] = np.array(info["resolution"][0][0][0])
    else:
        result["order"] = 1
        result["numCoeffs"] = 3
        result["resolution"] = np.array([70, 70])
        
    #print("\torder={}  numCoeffs={}  resolution={} {}".format(result["order"], result["numCoeffs"], result["resolution"][0], result["resolution"][1]))
        
    return result



def extract_coefficient_field( x, res, numCoeffs, coeff = 0 ):
	# returns 2d array containing the value for a specific coefficient
	# out of the solution vector
	res_x = res[0]
	res_y = res[1]

	u0 = np.zeros( (res_x, res_y), dtype=x.dtype )
	for voxel_i in range(res_x):
		for voxel_j in range(res_y):
			i = (voxel_j*res_x + voxel_i)*numCoeffs + coeff
			u0[voxel_i, voxel_j] = x[i, 0]
	return u0