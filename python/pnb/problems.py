# This module contains a set of problems for testing and studying.
# A problem essentially is a dictionary, containing fields for the different RTE parameters.

import numpy as np
import field
import pnbuilder







# This is the problem from the starmap paper. An emitting square at the center, surrounded by
# some highly absorbing squares, embedded within a scattering medium.
def checkerboard2d():

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

	# here we set some general parameters 
	problem["id"] = "checkerboard"
	problem["order"] = 1
	problem["domain"] = pnbuilder.Domain2D(7.0, 70)
	problem["domain_cpp"] = pnbuilder.Domain( np.array([7.0, 7.0]), np.array([70, 70]), np.array([0.0, 0.0]))
	problem["staggered"] = True

	# here we set the RTE fields
	problem["\\sigma_t"] = lambda pWS: sigma_a(pWS) + sigma_s(pWS)
	problem["\\sigma_a"] = sigma_a
	problem["\\sigma_s"] = sigma_s
	problem["f_p"] = phase_shcoeffs
	problem["q"] = source_shcoeffs

	return problem


# This basically takes an arbitrary problem and blurs it.
def blurred( problem, stddev = 1.0 ):

	domain = problem["domain"]

	problem["id"] = "{}_blur{}".format(problem["id"], stddev)

	sigma_t_grid = domain.rasterize(problem["\\sigma_t"])
	albedo_grid = domain.rasterize(problem["\\sigma_s"])/sigma_t_grid

	problem["\\sigma_a"] = field.VoxelGrid(domain, domain.rasterize(problem["\\sigma_a"]))
	problem["\\sigma_a"].blur(stddev)
	problem["\\sigma_s"] = field.VoxelGrid(domain, domain.rasterize(problem["\\sigma_s"]))
	problem["\\sigma_s"].blur(stddev)


	# fade to zero near border ----------------------------
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


	problem["\\sigma_t"] = field.VoxelGrid(domain, problem["\\sigma_a"].voxels+problem["\\sigma_s"].voxels)

	return problem

if __name__ == "__main__":
	pass















