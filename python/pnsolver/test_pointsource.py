import numpy as np
import renderer
import util
import json


# just for plotting
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable



def create_scene_pointlight( power, sigma_t, albedo ):
	volume = renderer.create_volume()
	volume.setBound( np.array([0.0, 0.0, 0.0]), np.array([2.0, 2.0, 2.0]) )
	sigma_t_field = renderer.create_constant_field3d( np.array([sigma_t, sigma_t, sigma_t]) )
	albedo_field = renderer.create_constant_field3d( np.array([albedo, albedo, albedo]) )
	volume.setExtinctionAlbedo( sigma_t_field, albedo_field )
	light = renderer.create_point_light( np.array([1.0, 1.0, 1.0]), np.array([power, power, power]) )
	return volume, light


def validate_pointsource_fluence():
	numSamples = 10000
	# assemble scene to render
	#camera = load_camera( "results/pointsource/result_pointsource.scn", "cam1" )
	integrator = renderer.create_simplept_integrator(True, -1)
	power = 4.0*np.pi
	sigma_t = 8.0
	albedo = 0.9
	volume, light = create_scene_pointlight(power, sigma_t, albedo)

	# get pnsolution
	#pns_p5 = renderer.load_pnsolution( "results/pointsource/pointsource_p5.pns" )
	#pns_p3 = renderer.load_pnsolution( "results/pointsource/pointsource_p3.pns" )
	pns_p1 = renderer.load_pnsolution( "results/pointsource/pointsource_p1_2.pns" )
	pns_p1_ms = renderer.load_pnsolution( "results/pointsource/pointsource_p1_2_ms.pns" )
	#pns_p1_vacuum = renderer.load_pnsolution( "results/pointsource/pointsource_p1_3.pns" )
	#pns_p1_ms64 = renderer.load_pnsolution( "results/pointsource/pointsource_p1_2_ms64.pns" )

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
	fluence_direct = renderer.compute_fluence( volume, light, integrator, points, seed, True, False, numSamples )
	fluence_indirect = renderer.compute_fluence( volume, light, integrator, points, seed, False, True, numSamples )
	#fluence_pns_p5 = power*renderer.compute_fluence_pnsolution( pns_p5, points )
	#fluence_pns_p3 = power*renderer.compute_fluence_pnsolution( pns_p3, points )
	fluence_pns_p1 = power*renderer.compute_fluence_pnsolution( pns_p1, points )
	fluence_pns_p1_ms = power*renderer.compute_fluence_pnsolution( pns_p1_ms, points )
	#fluence_pns_p1_ms64 = power*renderer.compute_fluence_pnsolution( pns_p1_ms64, points )

	r_list2 = np.linspace(1.0e-2, size*0.5, 100)
	points2 = np.zeros( (r_list2.shape[0], 3) )
	points2[:, 0] = r_list2+1.0
	points2[:, 1] = 1.0
	points2[:, 2] = 1.0
	groundtruth_ms = np.array([power*util.grosjean(r, sigma_t, albedo, direct_light=False) for r in r_list2 ])
	groundtruth_ss = np.array([power*util.grosjean(r, sigma_t, albedo, multiple_scattered_light=False) for r in r_list2 ])

	#solution_p5 = np.array([power*pns_p5.evalCoefficient(points2[i], 0) for i in range(points2.shape[0]) ])
	#solution_p3 = np.array([power*pns_p3.evalCoefficient(points2[i], 0) for i in range(points2.shape[0]) ])
	solution_p1 = np.array([power*pns_p1.evalCoefficient(points2[i], 0) for i in range(points2.shape[0]) ])
	solution_p1_ms = np.array([power*pns_p1_ms.evalCoefficient(points2[i], 0) for i in range(points2.shape[0]) ])
	#solution_p1_ms64 = np.array([power*pns_p1_ms64.evalCoefficient(points2[i], 0) for i in range(points2.shape[0]) ])
	#solution_p1_vacuum = np.array([power*pns_p1_vacuum.evalCoefficient(points2[i], 0) for i in range(points2.shape[0]) ])


	fig = plt.figure(figsize=(8,8));
	ax = fig.add_subplot(111)

	#plt.loglog( r_list2, groundtruth_ss, label="Grosjean" )
	#plt.loglog( r_list2, solution_p1_vacuum*np.sqrt(4.0*np.pi), label="P1 vacuum" )

	# pointsource with single emission voxel at center
	plt.loglog( r_list2, groundtruth_ss+groundtruth_ms, label="Grosjean SS+MS" )
	plt.loglog( r_list, fluence_direct+fluence_indirect, label="MonteCarlo SS+MS" )
	plt.loglog( r_list2, solution_p1*np.sqrt(4.0*np.pi), label="P1" )
	plt.loglog( r_list, fluence_pns_p1, label="P1 check", linestyle=" ", marker="." )

	# pointsource with single scattered emissionfield
	#plt.loglog( r_list2, groundtruth_ms, label="Grosjean MS" )
	#plt.loglog( r_list, fluence_indirect, label="MonteCarlo MS" )
	#plt.loglog( r_list2, solution_p1_ms*np.sqrt(4.0*np.pi), label="P1 MS" )
	#plt.loglog( r_list, fluence_pns_p1_ms, label="P1 MS check", linestyle=" ", marker="." )


	#plt.loglog( r_list2, groundtruth_ss, label="Grosjean SS" )
	#plt.loglog( r_list, fluence_direct, label="MonteCarlo SS" )
	#plt.loglog( r_list2, groundtruth_ms, label="Grosjean MS" )
	#plt.loglog( r_list, fluence_indirect, label="MonteCarlo MS" )

	#plt.loglog( r_list2, solution_p5, label="P5" )
	#plt.loglog( r_list, fluence_pns_p5, label="P5 check", linestyle=" ", marker="." )
	#plt.loglog( r_list2, solution_p3, label="P3" )
	#plt.loglog( r_list, fluence_pns_p3, label="P3 check", linestyle=" ", marker="." )

	# when using the 00 coefficient to present vacuum, we have to integrate over solid angle
	#plt.loglog( r_list2, solution_p1_ms*np.sqrt(4.0*np.pi), label="P1 MS" )
	#plt.loglog( r_list2, solution_p1_ms64, label="P1 MS 64" )
	#plt.loglog( r_list, fluence_pns_p1_ms, label="P1 MS check", linestyle=" ", marker="." )

	plt.legend(loc='best')

	plt.show()



if __name__ == "__main__":
	validate_pointsource_fluence()