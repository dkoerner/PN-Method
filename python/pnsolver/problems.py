# This module contains a set of problems for testing and studying.
# A problem essentially is a dictionary, containing fields for the different RTE parameters.

import numpy as np
import pnsolver
import util

from scipy.ndimage.filters import gaussian_filter




# This is the problem from the starmap paper. An emitting square at the center, surrounded by
# some highly absorbing squares, embedded within a scattering medium.
def checkerboard():

	def sigma_a( pWS ):
		x = pWS[0]
		y = pWS[1]
		cx = np.ceil(x)
		cy = np.ceil(y)
		g = 0
		if np.ceil((x+y)/2.0)*2.0 == (cx+cy) and cx > 1.0 and cx < 7.0 and cy > 1.0 and cy-2.0*np.abs(cx-4.0) < 4:
			g = 1
		return (1.0-g)*0 + g*10
		#return 0.0
		#return 5.0

	def sigma_s( pWS ):
		x = pWS[0]
		y = pWS[1]
		cx = np.ceil(x)
		cy = np.ceil(y)
		g = 0
		if np.ceil((x+y)/2.0)*2.0 == (cx+cy) and cx > 1.0 and cx < 7.0 and cy > 1.0 and cy-2.0*np.abs(cx-4.0) < 4:
			g = 1
		return (1.0-g)*1 + g*0
		#return 0.0
		#return 10.0 - sigma_a(pWS) # constant sigma_t
		#return 5.0


	def phase_shcoeffs( l, m, pWS ):
		if l == 0 and m == 0:
			return 1.0
		return 0.0

	def source_shcoeffs( l, m, pWS ):
		if l==0 and m == 0:
			x = pWS[0]
			y = pWS[1]
			if x > 3.0 and x < 4.0 and y > 3.0 and y < 4.0:
				return 1.0
				#return 100.0
			return 0.0
		return 0.0

	problem = {}

	# here we set some general parameters 
	size = 7.0
	#res = 20
	#res = 512
	#res = 64
	res = 128
	#res = 71
	#res = 2
	#res = 200
	#res = 150
	#res = 200
	domain = pnsolver.Domain( np.array([size, size, 1]), np.array([res, res, 1]), np.array([0.0, 0.0, 0.0]))
	#voxel_area = domain.voxelSize()[0]*domain.voxelSize()[1]/0.01
	#voxel_area = 1.0
	#print("voxel_area={}".format(voxel_area))
	#print("1/voxel_area={}".format(1.0/voxel_area))

	problem["id"] = "checkerboard"
	problem["domain"] = domain

	# RTE parameters ------------
	offset = np.array([1.0, 1.0, 1.0])
	problem["sigma_t"] = pnsolver.VoxelGridField( util.rasterize(lambda pWS: sigma_a(pWS) + sigma_s(pWS), domain, dtype=complex), domain, offset*0.5 )
	problem["sigma_a"] = pnsolver.VoxelGridField( util.rasterize(sigma_a, domain, dtype=complex), domain, offset*0.5 )
	problem["sigma_s"] = pnsolver.VoxelGridField( util.rasterize(sigma_s, domain, dtype=complex), domain, offset*0.5 )
	problem["f_p"] = [pnsolver.Constant(1.0)]
	problem["q"] = [pnsolver.VoxelGridField( util.rasterize(lambda pWS: source_shcoeffs(0, 0, pWS), domain, dtype=complex), domain, offset*0.5 )]

	return problem


def homogeneous():

	def sigma_a( pWS ):
		return 2.0

	def sigma_s( pWS ):
		return 2.0

	def phase_shcoeffs( l, m, pWS ):
		if l == 0 and m == 0:
			return 1.0
		return 0.0

	def source_shcoeffs( l, m, pWS ):
		if l==0 and m == 0:
			x = pWS[0]
			y = pWS[1]
			if x > 3.0 and x < 4.0 and y > 3.0 and y < 4.0:
				return 1.0
			return 0.0
		return 0.0

	problem = {}

	# here we set some general parameters 
	size = 7.0
	#res = 20
	#res = 512
	res = 64
	#res = 71
	#res = 2
	#res = 200
	#res = 150
	#res = 200
	domain = pnsolver.Domain( np.array([size, size, 1]), np.array([res, res, 1]), np.array([0.0, 0.0, 0.0]))
	#voxel_area = domain.voxelSize()[0]*domain.voxelSize()[1]/0.01
	#voxel_area = 1.0
	#print("voxel_area={}".format(voxel_area))
	#print("1/voxel_area={}".format(1.0/voxel_area))

	problem["id"] = "checkerboard"
	problem["domain"] = domain

	# RTE parameters ------------
	offset = np.array([1.0, 1.0, 1.0])
	problem["sigma_t"] = pnsolver.VoxelGridField( util.rasterize(lambda pWS: sigma_a(pWS) + sigma_s(pWS), domain, dtype=complex), domain, offset*0.5 )
	problem["sigma_a"] = pnsolver.VoxelGridField( util.rasterize(sigma_a, domain, dtype=complex), domain, offset*0.5 )
	problem["sigma_s"] = pnsolver.VoxelGridField( util.rasterize(sigma_s, domain, dtype=complex), domain, offset*0.5 )
	problem["f_p"] = [pnsolver.Constant(1.0)]
	problem["q"] = [pnsolver.VoxelGridField( util.rasterize(lambda pWS: source_shcoeffs(0, 0, pWS), domain, dtype=complex), domain, offset*0.5 )]

	return problem

# This is the problem from the starmap paper. An emitting square at the center, surrounded by
# some highly absorbing squares, embedded within a scattering medium.
def checkerboard3d():

	def sigma_a( pWS ):
		x = pWS[0]
		y = pWS[1]
		z = pWS[2]
		cx = np.ceil(x)
		cy = np.ceil(y)
		cz = np.ceil(z)
		g = 0
		if( np.ceil((x+y)/2.0)*2.0 == (cx+cy) and
			np.ceil((z+y)/2.0)*2.0 == (cz+cy) and
		    cx > 1.0 and cx < 7.0 and
		    cz > 1.0 and cz < 7.0 and
		    cy > 1.0 and
		    cy-2.0*np.abs(cx-4.0) < 4 and
		    cy-2.0*np.abs(cz-4.0) < 4 ):
			g = 1
		return (1.0-g)*0 + g*10

	def sigma_s( pWS ):
		x = pWS[0]
		y = pWS[1]
		z = pWS[2]
		cx = np.ceil(x)
		cy = np.ceil(y)
		cz = np.ceil(z)
		g = 0
		if( np.ceil((x+y)/2.0)*2.0 == (cx+cy) and
			np.ceil((z+y)/2.0)*2.0 == (cz+cy) and
		    cx > 1.0 and cx < 7.0 and
		    cz > 1.0 and cz < 7.0 and
		    cy > 1.0 and
		    cy-2.0*np.abs(cx-4.0) < 4 and
		    cy-2.0*np.abs(cz-4.0) < 4 ):
			g = 1
		return (1.0-g)*1 + g*0


	def phase_shcoeffs( l, m, pWS ):
		if l == 0 and m == 0:
			return 1.0
		return 0.0

	def source_shcoeffs( l, m, pWS ):
		if l==0 and m == 0:
			x = pWS[0]
			y = pWS[1]
			z = pWS[2]
			if( x > 3.0 and x < 4.0 and
				y > 3.0 and y < 4.0 and
				z > 3.0 and z < 4.0 ):
				return 1.0
			return 0.0
		return 0.0

	#problem = {}

	# here we set some general parameters 
	size = 7.0
	#res = 20
	res = 64
	#res = 71
	#res = 2
	#res = 200
	#res = 150
	#res = 200
	domain = pnsolver.Domain( np.array([size, size, size]), np.array([res, res, res]), np.array([0.0, 0.0, 0.0]))
	#voxel_area = domain.voxelSize()[0]*domain.voxelSize()[1]/0.01
	#voxel_area = 1.0
	#print("voxel_area={}".format(voxel_area))
	#print("1/voxel_area={}".format(1.0/voxel_area))

	#problem["id"] = "checkerboard"
	#problem["domain"] = domain

	# RTE parameters ------------
	#offset = np.array([1.0, 1.0, 1.0])
	#problem["sigma_t"] = pnsolver.VoxelGridField( util.rasterize(lambda pWS: sigma_a(pWS) + sigma_s(pWS), domain, dtype=complex), domain, offset*0.5 )
	#problem["sigma_a"] = pnsolver.VoxelGridField( util.rasterize(sigma_a, domain, dtype=complex), domain, offset*0.5 )
	#problem["sigma_s"] = pnsolver.VoxelGridField( util.rasterize(sigma_s, domain, dtype=complex), domain, offset*0.5 )
	#problem["f_p"] = [pnsolver.Constant(1.0)]
	#problem["q"] = [pnsolver.VoxelGridField( util.rasterize(lambda pWS: source_shcoeffs(0, 0, pWS), domain, dtype=complex), domain, offset*0.5 )]

	t = pnsolver.PNVolume( domain )

	absorption_voxels = util.rasterize3d( domain, lambda pWS: sigma_a(pWS) )
	scattering_voxels = util.rasterize3d( domain, lambda pWS: sigma_s(pWS) )
	extinction_voxels = absorption_voxels+scattering_voxels

	extinction_field = pnsolver.VoxelGridField3d(extinction_voxels)
	albedo_field = pnsolver.VoxelGridField3d(scattering_voxels/extinction_voxels)

	t.setExtinctionAlbedo( extinction_field, albedo_field )
	t.setPhase( 0, 0, pnsolver.ConstantField3d(np.array([1.0, 1.0, 1.0])) )
	emission_voxels = util.rasterize3d( domain, lambda pWS: source_shcoeffs(0, 0, pWS) )
	t.setEmission( 0, 0, pnsolver.VoxelGridField3d(emission_voxels) )

	return t

def pointsource3d():

	# here we set some general parameters 
	size = 2.0
	#res = 20
	#res = 100
	#res = 100
	res = 80
	#res = 71
	#res = 2
	#res = 200
	#res = 150
	#res = 200
	domain = pnsolver.Domain( np.array([size, size, size]), np.array([res, res, res]), np.array([0.0, 0.0, 0.0]))

	t = pnsolver.PNVolume( domain )

	# RTE parameters ------------
	sigma_t = 8.0
	albedo = 0.9
	extinction_field = pnsolver.ConstantField3d(np.array([sigma_t, sigma_t, sigma_t]))
	albedo_field = pnsolver.ConstantField3d(np.array([albedo, albedo, albedo]))
	t.setExtinctionAlbedo( extinction_field, albedo_field )
	t.setPhase( 0, 0, pnsolver.ConstantField3d(np.array([1.0, 1.0, 1.0])) )

	# the emission field is a 3d voxelgrid where the center voxel is set to an emission value
	shape = (domain.getResolution()[0], domain.getResolution()[1], domain.getResolution()[2], 3)
	q_voxels = np.zeros(shape) # real valued 
	vs = domain.getVoxelSize()
	pointsource_center_WS = np.array([size*0.5, size*0.5, size*0.5])
	pointsource_center_VS = domain.worldToVoxel(pointsource_center_WS)
	center_voxel = np.array([int(pointsource_center_VS[0]), int(pointsource_center_VS[1]), int(pointsource_center_VS[2])])
	q_voxels[center_voxel[0], center_voxel[1], center_voxel[2], :] = 1.0/(vs[0]*vs[1]*vs[2])
	t.setEmission( 0, 0, pnsolver.VoxelGridField3d(q_voxels) )

	return t


def nebulae():
	size = 1.0
	res = 64
	domain = pnsolver.Domain( np.array([size, size, size]), np.array([res, res, res]), np.array([0.0, 0.0, 0.0]))
	t = pnsolver.PNVolume( domain )

	#sigma_t = 8.0
	#extinction_field = pnsolver.ConstantField3d(np.array([sigma_t, sigma_t, sigma_t]))
	#extinction_field = pnsolver.load_voxelgridfield3d("c:/projects/epfl/epfl17/python/pnsolver/results/nebulae/nebulae_extinction.vgrid.3d")
	#extinction_field = pnsolver.load_bgeo("c:/projects/epfl/epfl17/python/pnsolver/results/nebulae/nebulae64.bgeo")
	extinction_field = pnsolver.load_voxelgridfield3d("c:/projects/epfl/epfl17/python/pnsolver/results/nebulae/nebulae64_extinction.grid")
	
	albedo = 0.9
	albedo_field = pnsolver.ConstantField3d(np.array([albedo, albedo, albedo]))
	t.setExtinctionAlbedo( extinction_field, albedo_field )

	t.setPhase( 0, 0, pnsolver.ConstantField3d(np.array([1.0, 1.0, 1.0])) )

	# the emissionfield has been generated by computing the single scattered light
	emission_field = pnsolver.load_voxelgridfield3d("c:/projects/epfl/epfl17/python/pnsolver/results/nebulae/nebulae64_emission.grid")
	t.setEmission( 0, 0, emission_field )

	return t





def vacuum():

	def sigma_a( pWS ):
		return 0.0

	def sigma_s( pWS ):
		return 0.0
		#return 10.0 - sigma_a(pWS) # constant sigma_t


	def phase_shcoeffs( l, m, pWS ):
		if l == 0 and m == 0:
			return 1.0
		return 0.0

	def source_shcoeffs( l, m, pWS ):
		if l==0 and m == 0:
			x = pWS[0]
			y = pWS[1]
			if x > 3.0 and x < 4.0 and y > 3.0 and y < 4.0:
				return 1.0
				#return 100.0
			return 0.0
		return 0.0

	problem = {}

	# here we set some general parameters 
	size = 7.0
	res = 70
	domain = pnsolver.Domain( np.array([size, size]), np.array([res, res]), np.array([0.0, 0.0]))

	problem["id"] = "vacuum"
	problem["domain"] = domain

	# RTE parameters ------------
	offset = np.array([1.0, 1.0])
	problem["sigma_t"] = pnsolver.VoxelGrid( util.rasterize(lambda pWS: sigma_a(pWS) + sigma_s(pWS), domain, dtype=complex), domain, offset*0.5 )
	problem["sigma_a"] = pnsolver.VoxelGrid( util.rasterize(sigma_a, domain, dtype=complex), domain, offset*0.5 )
	problem["sigma_s"] = pnsolver.VoxelGrid( util.rasterize(sigma_s, domain, dtype=complex), domain, offset*0.5 )
	problem["f_p"] = [pnsolver.Constant(1.0)]
	problem["q"] = [pnsolver.VoxelGrid( util.rasterize(lambda pWS: source_shcoeffs(0, 0, pWS), domain, dtype=complex), domain, offset*0.5 )]

	return problem


# This basically takes an arbitrary problem and blurs it.
def blurred( problem, stddev = 1.0 ):

	domain = problem["domain"]
	offset = np.array([0.5, 0.5])

	problem["id"] = "{}_blur{}".format(problem["id"], stddev)

	albedo = lambda pWS : np.real(problem["sigma_s"](pWS))/np.real(problem["sigma_t"](pWS))
	albedo_voxels = util.rasterize(albedo, domain)
	sigma_t = lambda pWS : np.real(problem["sigma_t"](pWS))
	sigma_t_voxels = util.rasterize(sigma_t, domain)

	# blur albedo and sigma_t voxels
	albedo_blurred = pnsolver.VoxelGrid(gaussian_filter(albedo_voxels, stddev).astype(complex), domain, offset)
	sigma_t_blurred = pnsolver.VoxelGrid(gaussian_filter(sigma_t_voxels, stddev).astype(complex), domain, offset)

	# reconstruct sigma_s and sigma_a from blurred albedo and sigma_t
	sigma_s_blurred = lambda pWS: albedo_blurred(pWS)*sigma_t_blurred(pWS)
	sigma_a_blurred = lambda pWS: (1.0-albedo_blurred(pWS))*sigma_t_blurred(pWS)

	problem["sigma_a"] = pnsolver.VoxelGrid(util.rasterize(sigma_a_blurred, domain, dtype=complex), domain, offset)
	problem["sigma_s"] = pnsolver.VoxelGrid(util.rasterize(sigma_s_blurred, domain, dtype=complex), domain, offset)
	problem["sigma_t"] = pnsolver.VoxelGrid(util.rasterize(sigma_t_blurred, domain, dtype=complex), domain, offset)


	# fade to zero near border ----------------------------
	'''
	# fade out sigma_a and sigma_s near the border
	def fade(pWS):
		d0 = np.abs(pWS[0] - domain.bound_min[0])
		d1 = np.abs(pWS[0] - domain.bound_max[0])
		d2 = np.abs(pWS[1] - domain.bound_min[1])
		d3 = np.abs(pWS[1] - domain.bound_max[1])
		d = min( d0, d1, d2, d3 )

		left = 0.2
		right = 0.8
		if d < left:
			return 0.0
		elif d > right:
			return 1.0
		t = (d-left)/(right-left)
		return t
		
	fade_field = domain.rasterize(fade)
	problem["\\sigma_a"].voxels *= fade_field
	problem["\\sigma_s"].voxels *= fade_field
	'''

	#problem["\\sigma_t"] = field.VoxelGrid(sigma_a_voxels+sigma_s_voxels, domain)

	return problem

if __name__ == "__main__":
	#pass
	test = checkerboard3d()















