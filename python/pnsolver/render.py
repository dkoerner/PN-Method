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




def create_scene_nebulae():
	albedo = 0.9

	#offset = -0.5
	offset = 0.0

	volume = renderer.create_volume()
	volume.setBound( np.array([-0.5, -0.5, -0.5+offset]), np.array([0.5, 0.5, 0.5+offset]) )


	sigma_t_field = renderer.load_bgeo("c:/projects/epfl/epfl17/python/pnsolver/results/nebulae/nebulae64.bgeo")
	#sigma_t_field = renderer.load_bgeo("c:/projects/epfl/epfl17/python/pnsolver/results/nebulae/nebulae200.bgeo")
	#sigma_t_field = renderer.load_bgeo("c:/projects/epfl/epfl17/python/pnsolver/results/nebulae/test.bgeo")
	#sigma_t_field.save("c:/projects/epfl/epfl17/python/pnsolver/results/nebulae/nebulae64_extinction.grid")

	#sigma_t_field = renderer.load_bgeo("test_sphere.bgeo")
	albedo_field = renderer.create_constant_field3d( np.array([albedo, albedo, albedo]) )
	volume.setExtinctionAlbedo( sigma_t_field, albedo_field )

	# point light
	#light_pos = np.array([-0.5, -0.5, -0.5])
	#power = 4.0*np.pi
	#light = renderer.create_point_light( light_pos, np.array([power, power, power]) )

	# directional light
	radiance = 1.0
	light_dir = np.array([0.0, -1.0, 0.0])
	light = renderer.create_directional_light( light_dir, np.array([radiance, radiance, radiance]) )

	return volume, light




def render_groundtruth_radial_distribution( pWS, filename ):

	integrator = renderer.create_simplept_integrator(True, -1)
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

	#renderer.test()
	#img = renderer.load_image("nebulae_p1_ms.exr").asMatrix()
	#renderer.save_exr("nebulae_p1_ms_check.exr", img*np.sqrt(4.0*np.pi)*np.sqrt(4.0*np.pi))
	#renderer.save_exr("nebulae_p1_ms_check.exr", img*4.0*np.pi)
	#renderer.save_exr("nebulae_p1_ms_check.exr", img*np.sqrt(4.0*np.pi))
	#*4.0*np.pi*np.sqrt(4.0*np.pi)
	#exit(1)

	'''
	bound_min = np.array([-0.5, -0.5, -0.5])
	bound_max = np.array([0.5, 0.5, 0.5])


	pns = renderer.load_pnsolution( "test.pns" )
	test = pns.getCoefficientField(0)
	print(np.min(test))
	print(np.max(test))
	gg = np.zeros((test.shape[0], test.shape[1], test.shape[2], 3))
	gg[:,:,:,0] = test
	gg[:,:,:,1] = test
	gg[:,:,:,2] = test
	renderer.save_bgeo("test.bgeo", renderer.VoxelGridField3d(gg), bound_min, bound_max)
	exit(1)
	'''


	'''
	sigma_t_field = renderer.load_bgeo("nebulae200.bgeo")
	renderer.save_bgeo("nebulae200_test.bgeo", sigma_t_field, bound_min, bound_max)
	exit(1)
	'''

	#'''

	# converting pn solution to bgeo file -----------------
	'''
	#pns = renderer.load_pnsolution( "c:/projects/epfl/epfl17/python/pnsolver/results/checkerboard/checkerboard_p1.pns" )
	pns = renderer.load_pnsolution( "c:/projects/epfl/epfl17/python/pnsolver/results/nebulae/nebulae_p1_3.pns" )
	#pns2 = renderer.load_pnsolution( "c:/projects/epfl/epfl17/python/pnsolver/results/checkerboard/checkerboard_p1_2.pns" )
	#pns = renderer.load_pnsolution( "c:/projects/epfl/epfl17/python/pnsolver/results/pointsource/pointsource_p1.pns" )
	img = pns.getCoefficientField(0)
	result = np.zeros( (img.shape[0], img.shape[1], img.shape[2], 3) )
	result[:,:,:,0] = img
	result[:,:,:,1] = img
	result[:,:,:,2] = img

	renderer.save_bgeo( "c:/projects/epfl/epfl17/python/pnsolver/results/nebulae/nebulae_p1_3.pns.bgeo", renderer.VoxelGridField3d(result), bound_min, bound_max )
	exit(1)
	'''


	# visualize slice of coefficient field from pn solution --------------
	'''
	fig = plt.figure(figsize=(8,8));

	ax = fig.add_subplot(121)
	#img_view = plt.imshow(img, origin='lower',zorder=1, interpolation="nearest", cmap='jet', norm=LogNorm(vmin=np.min(img), vmax=np.max(img)))
	img_view = plt.imshow(img, origin='lower',zorder=1, interpolation="nearest", cmap='jet', vmin=np.min(img), vmax=np.max(img))

	#ax = fig.add_subplot(122)
	#img_view = plt.imshow(img2, origin='lower',zorder=1, interpolation="nearest", cmap='jet', norm=LogNorm(vmin=np.min(img2), vmax=np.max(img2)))
	plt.show()

	exit(1)
	'''



	# compute unscattered fluence field for nebulae dataset and store to disk -----------
	'''
	volume, light = create_scene_nebulae()

	numSamples = 100

	# render unscattered light into emission field for PN solver ---
	unscattered_fluence_field = renderer.compute_unscattered_fluence( volume, light, numSamples, np.array([64, 64, 64]) )
	bound_min = np.array([-0.5, -0.5, -0.5])
	bound_max = np.array([0.5, 0.5, 0.5])
	#renderer.save_bgeo("c:/projects/epfl/epfl17/python/pnsolver/results/nebulae/nebulae64_emission.bgeo", unscattered_light_field, bound_min, bound_max)
	unscattered_fluence_field.save("c:/projects/epfl/epfl17/python/pnsolver/results/nebulae/nebulae64_unscattered_fluence.grid")
	#renderer.save_bgeo("c:/projects/epfl/epfl17/python/pnsolver/results/nebulae/nebulae200_emission.bgeo", unscattered_light_field, bound_min, bound_max)
	#unscattered_light_field.save("c:/projects/epfl/epfl17/python/pnsolver/results/nebulae/nebulae200_emission.grid")

	exit(1)
	'''

	# render pn-solution directly (instead of using it to boost standard MC) -----------------------------
	'''
	volume, light = create_scene_nebulae()
	camera = load_camera("c:/projects/epfl/epfl17/python/pnsolver/results/nebulae/nebulae.scn")
	#camera = renderer.create_perspective_camera(512, 512, 45.0)
	#translate = translation_matrix( np.array([0.0, 0.0, 1.9]) )
	#camera.setCameraToWorldTransform(translate)

	numSamples = 100
	#numSamples = 1
	
	#integrator_groundtruth_ss = renderer.create_simplept_integrator(True, 1)
	#img_ss = renderer.render( volume, light, camera, integrator_groundtruth_ss, numSamples )
	#img_ss.save("nebulae_groundtruth_ss.exr")
	#integrator_groundtruth_ms = renderer.create_simplept_integrator(False, -1)
	#img_ms = renderer.render( volume, light, camera, integrator_groundtruth_ms, numSamples )
	#img_ms.save("nebulae_groundtruth_ms.exr")


	pns = renderer.load_pnsolution( "c:/projects/epfl/epfl17/python/pnsolver/results/nebulae/nebulae_p1_2_ms.pns" )	
	#pns = renderer.load_pnsolution( "c:/projects/epfl/epfl17/python/pnsolver/results/nebulae/nebulae_p3_2_ms.pns" )	
	#pns = renderer.load_pnsolution( "c:/projects/epfl/epfl17/python/pnsolver/results/nebulae/nebulae_p5_2_ms.pns" )
	#pns = renderer.load_pnsolution( "c:/projects/epfl/epfl17/python/pnsolver/results/nebulae/nebulae_p5.pns" )	
	#pns = renderer.load_pnsolution( "c:/projects/epfl/epfl17/python/pnsolver/results/nebulae/nebulae_cda3.pns" )	

	#integrator_pn_ss = renderer.create_directpn_integrator(pns, True, False)
	#img_ss = renderer.render( volume, light, camera, integrator_pn_ss, numSamples )
	#img_ss.save("nebulae_p1_ss.exr")

	integrator_pn_ms = renderer.create_directpn_integrator(pns, False, True)
	img_ms = renderer.render( volume, light, camera, integrator_pn_ms, numSamples )
	img_ms.save("nebulae_p1_3_ms.exr")
	#img_ms.save("nebulae_p3_3_ms.exr")
	#img_ms.save("nebulae_p5_3_ms.exr")
	#img_ms.save("nebulae_p3_ms.exr")
	#img_ms.save("nebulae_p5_ms.exr")
	#img_ms.save("nebulae_cda_ms.exr")
	exit(1)
	'''


	# boost standard MC using PN-Solution -----------------------------
	volume, light = create_scene_nebulae()
	camera = load_camera("c:/projects/epfl/epfl17/python/pnsolver/results/nebulae/nebulae.scn")

	integrator_groundtruth_ms = renderer.create_simplept_integrator(False, -1)

	pns = renderer.load_pnsolution( "c:/projects/epfl/epfl17/python/pnsolver/results/nebulae/nebulae_p5_2_ms.pns" )
	integrator_pnis = renderer.create_pnispt_integrator(pns, False, -1)

	numSamples = 1


	#img_ms = renderer.render( volume, light, camera, integrator_groundtruth_ms, numSamples )
	#img_ms.save("nebulae_ms_mc.exr")
	img_ms = renderer.render( volume, light, camera, integrator_pnis, numSamples )
	img_ms.save("nebulae_ms_pnis.exr")

	#'''
	depth = integrator_pnis.dbgGet("depth")
	throughput_over_pdf = integrator_pnis.dbgGet("throughput_over_pdf")
	L = integrator_pnis.dbgGet("L")
	phase_over_pdf = integrator_pnis.dbgGet("phase_over_pdf")
	phase_sampling = integrator_pnis.dbgGet("phase_sampling")
	phase_pdf = integrator_pnis.dbgGet("phase_pdf")
	pWS_x = integrator_pnis.dbgGet("pWS_x")
	pWS_y = integrator_pnis.dbgGet("pWS_y")
	pWS_z = integrator_pnis.dbgGet("pWS_z")
	dir_x = integrator_pnis.dbgGet("dir_x")
	dir_y = integrator_pnis.dbgGet("dir_y")
	dir_z = integrator_pnis.dbgGet("dir_z")
	print("depth=")
	print(depth)
	print("throughput_over_pdf=")
	print(throughput_over_pdf)
	print("L=")
	print(L)
	print("phase_over_pdf=")
	print(phase_over_pdf)
	print("phase_sampling=")
	print(phase_sampling)
	print("phase_pdf=")
	print(phase_pdf)
	print("pdf_shsampling=")
	print(integrator_pnis.dbgGet("pdf_shsampling"))
	
	
	#print(direction_pdf)
	#print(pWS_x)
	#print(pWS_y)
	#print(pWS_z)
	#print(dir_x)
	#print(dir_y)
	#print(dir_z)
	#'''

	# render radial fluence distribution at a given worldspace position =======
	'''
	numSamples = 2500
	pWS = np.array([-0.27704304, 0.36083166, -0.22953043])
	camera = renderer.create_sphere_camera( pWS, 128, 64 )
	integrator_groundtruth_ms = renderer.create_simplept_integrator(True, -1)
	img_ms = renderer.render( volume, light, camera, integrator_groundtruth_ms, numSamples )
	img_ms.save("nebulae_dbg_fluence.exr")
	'''




	'''
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


	#integrator = renderer.create_simplept_integrator(True, -1)
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
	'''


