import numpy as np
import renderer
import util
import json



# just for plotting
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable




# openexr

def translation_matrix(direction):
	# Return matrix to translate by direction vector.
	M = np.identity(4)
	M[:3, 3] = direction[:3]
	return M

def rotation_matrix(angle, direction):
	# Return matrix to rotate about axis defined by point and direction.
	sina = np.sin(angle)
	cosa = np.cos(angle)
	#direction = unit_vector(direction[:3])
	# rotation matrix around unit vector
	R = np.diag([cosa, cosa, cosa])
	R += np.outer(direction, direction) * (1.0 - cosa)
	direction *= sina
	R += np.array([[ 0.0,         -direction[2],  direction[1]],
	                  [ direction[2], 0.0,          -direction[0]],
	                  [-direction[1], direction[0],  0.0]])
	M = np.identity(4)
	M[:3, :3] = R
	return M

def load_camera( filename, id ):
	scn_file = filename
	json_file = open(scn_file)
	dd = json.load(json_file)
	json_file.close()
	cam = dd["cameras"][id]
	resx = cam["resx"]
	resy = cam["resx"]

	# http://www.sidefx.com/docs/houdini/ref/cameralenses
	pixel_aspect = 1.0
	focal_length = cam["camera.fl"]
	aperture_x = cam["camera.horizontalFilmAperture"]
	aperture_y = (resy*aperture_x) / (resx*pixel_aspect)
	hfovy_deg = np.rad2deg(2.0*np.arctan( (aperture_y/2.0) / focal_length ))

	# assemble localtransform
	translate = translation_matrix( np.array([cam["transform.tx"], cam["transform.ty"], cam["transform.tz"]]) )
	rotate_x = rotation_matrix( np.deg2rad(cam["transform.rx"]), np.array([1.0, 0.0, 0.0]) )
	rotate_y = rotation_matrix( np.deg2rad(cam["transform.ry"]), np.array([0.0, 1.0, 0.0]) )
	rotate_z = rotation_matrix( np.deg2rad(cam["transform.rz"]), np.array([0.0, 0.0, 1.0]) )

	cam2world = translate.dot(rotate_z.dot(rotate_y.dot(rotate_x)))

	camera = renderer.create_perspective_camera(resx, resy, hfovy_deg)
	camera.setCameraToWorldTransform(cam2world)

	return camera


def create_scene_pointlight():
	volume = renderer.create_volume()
	volume.setBound( np.array([0.0, 0.0, 0.0]), np.array([2.0, 2.0, 2.0]) )
	sigma_t = 8.0
	sigma_t_field = renderer.create_constant_field3d( np.array([sigma_t, sigma_t, sigma_t]) )
	albedo = 0.9
	albedo_field = renderer.create_constant_field3d( np.array([albedo, albedo, albedo]) )
	volume.setExtinctionAlbedo( sigma_t_field, albedo_field )
	power = 4.0*np.pi
	light = renderer.create_point_light( np.array([1.0, 1.0, 1.0]), np.array([power, power, power]) )

	return volume, light




if __name__ == "__main__":

	




	# assemble scene to render
	camera = load_camera( "results/pointsource/result_pointsource.scn", "cam1" )
	integrator = renderer.create_simplept_integrator()
	volume, light = create_scene_pointlight()

	# get pnsolution
	pns_p5 = renderer.load_pnsolution( "results/pointsource/pointsource_p5.pns" )
	pns_p3 = renderer.load_pnsolution( "results/pointsource/pointsource_p3.pns" )
	pns_p1 = renderer.load_pnsolution( "results/pointsource/pointsource_p1.pns" )

	# render
	#renderer.render( volume, light, camera, integrator ).save("test.exr")

	# fluence
	size = 2.0
	r_list = np.linspace(1.0e-2, size*0.5, 100)
	#r_list = np.array([0.1])

	points = np.zeros( (r_list.shape[0], 3) )
	points[:, 0] = r_list+1.0
	points[:, 1] = 1.0
	points[:, 2] = 1.0

	seed = 123
	fluence = renderer.compute_fluence( volume, light, integrator, points, seed )
	fluence_pns_p5 = power*renderer.compute_fluence_pnsolution( pns_p5, points )
	fluence_pns_p3 = power*renderer.compute_fluence_pnsolution( pns_p3, points )
	fluence_pns_p1 = power*renderer.compute_fluence_pnsolution( pns_p1, points )

	#maybe the phase function normalization is the problem here...


	r_list2 = np.linspace(1.0e-2, size*0.5, 100)
	points2 = np.zeros( (r_list2.shape[0], 3) )
	points2[:, 0] = r_list2+1.0
	points2[:, 1] = 1.0
	points2[:, 2] = 1.0
	groundtruth_ms = np.array([power*util.grosjean(r, sigma_t, albedo, direct_light=False) for r in r_list2 ])
	groundtruth_ss = np.array([power*util.grosjean(r, sigma_t, albedo, multiple_scattered_light=False) for r in r_list2 ])

	solution_p5 = np.array([power*pns_p5.evalCoefficient(points2[i], 0) for i in range(points2.shape[0]) ])
	solution_p3 = np.array([power*pns_p3.evalCoefficient(points2[i], 0) for i in range(points2.shape[0]) ])
	solution_p1 = np.array([power*pns_p1.evalCoefficient(points2[i], 0) for i in range(points2.shape[0]) ])


	fig = plt.figure(figsize=(8,8));
	ax = fig.add_subplot(111)

	#plt.loglog( r_list2, groundtruth_ss, label="grosjean (single scattering term)" )
	#plt.loglog( r_list, fluence, label="MonteCarlo" )

	#plt.loglog( r_list2, groundtruth_ms, label="grosjean (multiple scattering term)" )
	#plt.loglog( r_list, fluence, label="MonteCarlo" )

	plt.loglog( r_list2, groundtruth_ss+groundtruth_ms, label="Grosjean" )
	plt.loglog( r_list, fluence, label="MonteCarlo" )
	plt.loglog( r_list2, solution_p5, label="P5" )
	plt.loglog( r_list, fluence_pns_p5, label="P5 check", linestyle=" ", marker="." )
	plt.loglog( r_list2, solution_p3, label="P3" )
	plt.loglog( r_list, fluence_pns_p3, label="P3 check", linestyle=" ", marker="." )
	plt.loglog( r_list2, solution_p1, label="P1" )
	plt.loglog( r_list, fluence_pns_p1, label="P1 check", linestyle=" ", marker="." )

	plt.legend(loc='best')

	plt.show()
	fig.savefig("foo.pdf", bbox_inches='tight')