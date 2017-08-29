# This module contains various helper and utility functions such as 
# related to spherical harmonics, I/O etc.

import numpy as np
import pnbuilder
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

def write_pn_system(pnb, problem, A, b):
	data = {}
	if not A is None:
		data['A'] = A
	if not b is None:
		data['b'] = b.reshape((pnb.domain.numVoxels*pnb.numCoeffs, 1))
	data['pnb_info'] = pnb.get_info()

	#filename = "C:/projects/epfl/epfl17/python/sopn/data_{}.mat".format(problem["id"])
	#scipy.io.savemat(filename, data)

	data['sigma_s'] = rasterize(problem["\\sigma_s"], pnb.domain)
	data['sigma_a'] = rasterize(problem["\\sigma_a"], pnb.domain)
	data['sigma_t'] = rasterize(problem["\\sigma_t"], pnb.domain)
	data['q']       = rasterize(lambda pWS: problem['q'](0,0,pWS), pnb.domain)

	filename = "C:/projects/epfl/epfl17/python/sopn/system_{}.mat".format(problem["id"])
	print("writing PN system to {}".format(filename))
	scipy.io.savemat(filename, data)



def load_pn_solution( filename ):
	print("loading PN solution from {}".format(filename))
	data = scipy.io.loadmat(filename)
	pnb_info = {}

	pnb_info["order"] = data["pnb_info"]["order"][0][0][0][0]
	pnb_info["numCoeffs"] = data["pnb_info"]["numCoeffs"][0][0][0][0]
	pnb_info["domain_size"] = float(data["pnb_info"]["domain_size"][0][0][0][0])
	pnb_info["domain_res"] = data["pnb_info"]["domain_res"][0][0][0][0]
	pnb_info["coeff_offsets"] = data["pnb_info"]["coeff_offsets"][0][0]

	pnb = pnbuilder.from_info(pnb_info)


	result = {}
	if "x" in data:
		x_real = data["x"].reshape((pnb.domain.num_voxels()*pnb.num_coeffs()))
		result["x_real"] = x_real
	if "b" in data:
		b_real = data["b"].reshape((pnb.domain.num_voxels()*pnb.num_coeffs()))
		result["b_real"] = b_real
	if "A" in data:
		A_real = data["A"]
		result["A_real"] = A_real
	
	
	result["pnb"] = pnb
	result["sigma_t"] = data["sigma_t"]
	result["sigma_a"] = data["sigma_a"]
	result["sigma_s"] = data["sigma_s"]
	result["q"] = data["q"]

	#return A_real, x_real, b_real, pnb
	return result

	'''
	# convert to complex variables
	x_complex = pnb.to_complex(x_real)

	# now construct field of SHexpansions...
	coeff_fields = []
	for i in range(numCoeffs):
		offset = coeff_offsets[i]*0.5
		coeff_fields.append( problems.CoefficientGrid(domain, numCoeffs, i, offset, x_complex) )
	return problems.SHEXP(order, coeff_fields), domain, x_real
	'''
