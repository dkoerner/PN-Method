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

def load_camera( filename, id = "cam1" ):
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


def create_scene_pointlight( power, sigma_t, albedo ):
	volume = renderer.create_volume()
	volume.setBound( np.array([0.0, 0.0, 0.0]), np.array([2.0, 2.0, 2.0]) )
	sigma_t_field = renderer.create_constant_field3d( np.array([sigma_t, sigma_t, sigma_t]) )
	albedo_field = renderer.create_constant_field3d( np.array([albedo, albedo, albedo]) )
	volume.setExtinctionAlbedo( sigma_t_field, albedo_field )
	light = renderer.create_point_light( np.array([1.0, 1.0, 1.0]), np.array([power, power, power]) )
	return volume, light

def create_scene_nebulae():
	light_pos = np.array([-0.5, -0.5, -0.5])
	power = 4.0*np.pi
	albedo = 0.9

	radiance = 1.0
	light_dir = np.array([0.0, -1.0, 0.0])

	#offset = -0.5
	offset = 0.0

	volume = renderer.create_volume()
	volume.setBound( np.array([-0.5, -0.5, -0.5+offset]), np.array([0.5, 0.5, 0.5+offset]) )
	sigma_t_field = renderer.load_bgeo("nebulae200.bgeo")
	#sigma_t_field = renderer.load_bgeo("test_sphere.bgeo")
	albedo_field = renderer.create_constant_field3d( np.array([albedo, albedo, albedo]) )
	volume.setExtinctionAlbedo( sigma_t_field, albedo_field )
	#light = renderer.create_point_light( light_pos, np.array([power, power, power]) )
	light = renderer.create_directional_light( light_dir, np.array([radiance, radiance, radiance]) )
	return volume, light


def validate_pointsource_fluence():
	# assemble scene to render
	#camera = load_camera( "results/pointsource/result_pointsource.scn", "cam1" )
	integrator = renderer.create_simplept_integrator()
	power = 4.0*np.pi
	sigma_t = 8.0
	albedo = 0.9
	volume, light = create_scene_pointlight(power, sigma_t, albedo)

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


def render_groundtruth_radial_distribution( pWS, filename ):

	integrator = renderer.create_simplept_integrator()
	camera = renderer.create_sphere_camera( pWS, 128, 64 )

	# render groundtruth radial distribution of radiance at pWS into an image
	numSamples = 5000
	#numSamples = 1

	power = 4.0*np.pi
	sigma_t = 8.0
	albedo = 0.9

	volume, light = create_scene_pointlight(power, sigma_t, albedo)
	img_radiance = renderer.render( volume, light, camera, integrator, numSamples )
	img_radiance.save(filename)
	#exit(1)


if __name__ == "__main__":




	integrator = renderer.create_simplept_integrator()
	camera = load_camera("c:/projects/epfl/epfl17/python/pnsolver/results/nebulae/nebulae.scn")
	#camera = renderer.create_perspective_camera(512, 512, 45.0)
	volume, light = create_scene_nebulae()
	numSamples = 100
	img = renderer.render( volume, light, camera, integrator, numSamples )
	img.save("test.exr")

	#fig = plt.figure(figsize=(8,8));
	#ax = fig.add_subplot(111)
	#sigma_t = volume.getSlice(0.5)
	#img_view = plt.imshow(sigma_t, origin='lower',zorder=1, interpolation="nearest", cmap='gray', vmin=0.0, vmax=60.0)
	#plt.show()



	exit(1)


	#r_list = [0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9]
	#r_list = [0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9]
	#count = 2
	#for r in r_list:
	#	pWS = np.array([1.0+r, 1.0, 1.0])
	#	filename = "test.{}.exr".format(count)
	#	render_groundtruth_radial_distribution(pWS, filename)
	#	count = count+1
	#exit(1)

	
	power = 4.0*np.pi
	sigma_t = 8.0
	albedo = 0.9
	volume, light = create_scene_pointlight(power, sigma_t, albedo)


	#integrator = renderer.create_simplept_integrator()
	pns = renderer.load_pnsolution( "results/pointsource/pointsource_p5.pns" )
	integrator = renderer.create_pnispt_integrator(pns)


	res = 256
	points = np.zeros( (res, 3) )
	points[:, 0] = 1.7
	points[:, 1] = 1.0
	points[:, 2] = 1.0

	fluence = np.zeros( (res, res) )
	for i in range(res):
		print("computing fluence i={}".format(i))
		fluence[i,:] = renderer.compute_fluence( volume, light, integrator, points, 123+i*8 )


	renderer.save_exr("fluence2.exr", fluence)
	fig = plt.figure(figsize=(8,8));
	ax = fig.add_subplot(111)

	img_view = plt.imshow(fluence, origin='lower',zorder=1, interpolation="nearest", cmap='jet', vmin=np.min(fluence), vmax=np.max(fluence))

	#plt.plot( fluence )
	#plt.plot( fluence2, marker=".", linestyle=None )
	
	#plt.loglog( r_list2, groundtruth_ss+groundtruth_ms, label="Grosjean" )
	#plt.loglog( r_list, fluence, label="MonteCarlo" )
	#plt.loglog( r_list2, solution_p5, label="P5" )
	#plt.loglog( r_list, fluence_pns_p5, label="P5 check", linestyle=" ", marker="." )
	#plt.loglog( r_list2, solution_p3, label="P3" )
	#plt.loglog( r_list, fluence_pns_p3, label="P3 check", linestyle=" ", marker="." )
	#plt.loglog( r_list2, solution_p1, label="P1" )
	#plt.loglog( r_list, fluence_pns_p1, label="P1 check", linestyle=" ", marker="." )

	#plt.legend(loc='best')


	plt.show()


