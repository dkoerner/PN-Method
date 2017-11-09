import numpy as np
import renderer
import util
import json
import random

import render

# just for plotting
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable








def plot_spherical_function( fun, vmin=-0.5, vmax=0.5 ):
    # plot sh functions
    theta = np.arange(0.0, 1.0, 0.01)*np.pi
    phi = np.arange(0.0, 1.0, 0.01)*2.0*np.pi

    f_img = np.zeros((theta.shape[0], phi.shape[0]))
    for j in range(phi.shape[0]):
        for i in range(theta.shape[0]):
            f_img[i, j] = fun(theta[i], phi[j])
    #return plt.imshow(f_img, interpolation="nearest", cmap='jet', vmin=vmin, vmax=vmax)
    vmin = np.min(f_img)
    vmax = np.max(f_img)
    print("vmin={} vmax={}".format(vmin, vmax))
    return plt.imshow(f_img, origin='lower',zorder=1, interpolation="nearest", cmap='jet', vmin=vmin, vmax=vmax, extent=[0.0,2.0*np.pi,0.0,np.pi])






if __name__ == "__main__":

	pWS = np.array([1.1, 1.0, 1.0])

	
	

	'''
	# render groundtruth radial distribution of radiance at pWS into an image
	#numSamples = 5000
	numSamples = 1
	volume, light = render.create_scene_pointlight()
	integrator = renderer.create_simplept_integrator()
	camera = renderer.create_sphere_camera( pWS, 128, 64 )
	img_radiance = renderer.render( volume, light, camera, integrator, numSamples )
	img_radiance.save("test.exr")
	#exit(1)
	'''

	img_groundtruth = renderer.load_image("test.exr")
	img_groundtruth = renderer.blur_image(img_groundtruth, 3.0)
	img_groundtruth_sampler = renderer.create_image_sampler(img_groundtruth)



	
	#pns_p1.evalCoefficients()

	#order = 1
	#block_depth = 7
	#sh = renderer.create_shtest(order, block_depth)
	#res = 128

	#sh.test_basis_integrals()
	#theta, phi, pdf, r1, r2 = sh.sample()
	#sh.sample2(r1, r2, False)

	#img_pdf = np.zeros((res, res))
	'''
	samples2_theta = []
	samples2_phi = []

	for i in range(res):
		for j in range(res):
			debug = False

			if i == 0 and j==0:
				debug = True

			r1 = i/float(res)
			r2 = j/float(res)
			#pdf = sh.pdf( r1, r2 )
			#img_pdf[i, j] = pdf

			theta, phi, pdf, a, b = sh.sample( r1, r2 )
			#samples2_theta.append(theta)
			#samples2_phi.append(phi)
			samples2_theta.append(theta)
			samples2_phi.append(phi)
	'''

	#print( "pdf min={} max={}".format(np.min(img_pdf), np.max(img_pdf)) )


	#'''
	# visualize samples and the sampled sh function ------------------
	samples_theta = []
	samples_phi = []

	numSamples = 10000
	for i in range(numSamples):
		#theta, phi, pdf, r1, r2 = sh.sample()
		#d = -pns_p1.sample(pWS, np.array([random.random(), random.random()]))
		#theta, phi = util.sphericalCoordinates(d)

		#pRaster = img_groundtruth_sampler.sample( random.random(), random.random() )
		#pRaster+= np.array([random.random(), random.random()])
		#uv = img_radiance.rasterToUV(pRaster)
		#theta = uv[1]*np.pi
		#phi = uv[0]*2.0*np.pi
		#samples_theta.append(theta)
		#samples_phi.append(phi)
		pass
	#for i in range(res):
	#	for j in range(res):
	#i = res-10
	#j = res-10
	#r1 = i/float(res)
	#r2 = j/float(res)
	#theta, phi, pdf, a, b = sh.sample2( r1, r2 )
	#if (i==res-10) and (j==res-10):
	#	print(sh.eval(theta, phi))
	#samples_theta.append(theta)
	#samples_phi.append(phi)


	#print("theta min={} max={}".format(np.min(samples_theta), np.max(samples_theta)))
	#print("phi min={} max={}".format(np.min(samples_phi), np.max(samples_phi)))

	def clamped_groundtruth( theta, phi ):
		uv = np.array([phi/(2.0*np.pi), theta/np.pi])
		pRaster = img_groundtruth.uvToRaster(uv)
		value = img_groundtruth.eval( pRaster )[0]
		if value < 0.0:
			value = 0.0
		return value

	
	def clamped_sh( pns, theta, phi ):
		# NB: the negative direction...
		value = pns.eval(pWS, -util.sphericalDirection(theta, phi))
		if value < 0.0:
			value = 0.0
		return value

	for i in range(4):
		pns = renderer.load_pnsolution( "results/pointsource/pointsource_p{}.pns".format(i+1) )

		fig = plt.figure(figsize=(8,8));
		ax = fig.add_subplot(121)

		img_view = plot_spherical_function( lambda theta, phi: clamped_groundtruth(theta, phi) )
		#plt.scatter(samples_phi, samples_theta, c='r', label="sdda",zorder=2, marker='.', s=1.0)
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)
		plt.colorbar(img_view, cax=cax)

		ax = fig.add_subplot(122)
		img_view = plot_spherical_function( lambda theta, phi: clamped_sh(pns, theta, phi) )
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)
		plt.colorbar(img_view, cax=cax)



		plt.show()
	#'''

	# visualize block integrals ------------------
	'''
	block_depth = 1
	fig = plt.figure(figsize=(8,8));
	ax = fig.add_subplot(121)

	#sh.setDebugLM(0, 0)
	img_blocks_org = sh.get_blocks(block_depth, True)
	img_view = plt.imshow(img_blocks_org, origin='lower',zorder=1, interpolation="nearest", cmap='jet', vmin=np.min(img_blocks_org), vmax=np.max(img_blocks_org), extent=[0.0,2.0*np.pi,0.0,np.pi])

	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	plt.colorbar(img_view, cax=cax)

	ax = fig.add_subplot(122)
	img_blocks = sh.get_blocks(block_depth, False)
	img_view = plt.imshow(img_blocks, origin='lower',zorder=1, interpolation="nearest", cmap='jet', vmin=np.min(img_blocks), vmax=np.max(img_blocks), extent=[0.0,2.0*np.pi,0.0,np.pi])

	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	plt.colorbar(img_view, cax=cax)

	plt.show()
	'''

	# visualize block integrals for all lm ------------------
	'''
	for l in range(2):
		for m in range(-l, l+1):
			sh.setDebugLM(l, m)
			fig = plt.figure(figsize=(8,8));
			ax = fig.add_subplot(121)

			img_blocks_org = sh.get_blocks(block_depth, True)
			img_view = plt.imshow(img_blocks_org, origin='lower',zorder=1, interpolation="nearest", cmap='jet', vmin=np.min(img_blocks_org), vmax=np.max(img_blocks_org), extent=[0.0,2.0*np.pi,0.0,np.pi])

			divider = make_axes_locatable(ax)
			cax = divider.append_axes("right", size="5%", pad=0.05)
			plt.colorbar(img_view, cax=cax)

			ax = fig.add_subplot(122)
			img_blocks = sh.get_blocks(block_depth, False)
			img_view = plt.imshow(img_blocks, origin='lower',zorder=1, interpolation="nearest", cmap='jet', vmin=np.min(img_blocks), vmax=np.max(img_blocks), extent=[0.0,2.0*np.pi,0.0,np.pi])

			divider = make_axes_locatable(ax)
			cax = divider.append_axes("right", size="5%", pad=0.05)
			plt.colorbar(img_view, cax=cax)

			plt.show()
	'''




	#img_view = plot_spherical_function( lambda theta, phi: clamped_sh(theta, phi) )
	#divider = make_axes_locatable(ax)
	#cax = divider.append_axes("right", size="5%", pad=0.05)
	#plt.colorbar(img_view, cax=cax)

	#ax = fig.add_subplot(212)
	#plot_spherical_function( lambda theta, phi: theta,  )



	#img_view = plot_spherical_function( lambda theta, phi: clamped_sh(theta, phi) )
	#img_view = plt.imshow(img_pdf, origin='lower',zorder=1, interpolation="nearest", cmap='jet', vmin=np.min(img_pdf), vmax=np.max(img_pdf), extent=[0.0,2.0*np.pi,0.0,np.pi])
	#img_view = plt.imshow(img_blocks, origin='lower',zorder=1, interpolation="nearest", cmap='jet', vmin=np.min(img_blocks), vmax=np.max(img_blocks), extent=[0.0,2.0*np.pi,0.0,np.pi])
	#img_view = plt.imshow(img_pdf, origin='lower',zorder=1, interpolation="nearest", cmap='jet', vmin=np.min(img_pdf), vmax=np.max(img_pdf), extent=[0.0,1.0])

	#plt.scatter(samples_phi, samples_theta, c='r', label="sdda",zorder=2, marker='.', s=5.0)
	#plt.xlim(0.0, 1.0)
	#plt.ylim(0.0, 1.0)


	#divider = make_axes_locatable(ax)
	#cax = divider.append_axes("right", size="5%", pad=0.05)
	#plt.colorbar(img_view, cax=cax)

	#plt.show()
	#fig.savefig("foo.pdf", bbox_inches='tight')

	#fig = plt.figure(figsize=(8,8));
	#ax = fig.add_subplot(111)
	#phi_list = np.linspace(0.0, np.pi, 200)
	#plt.plot(phi_list, [seperated_phi(13, phi) for phi in phi_list])
	#plt.show()

	#print(integrate_seperated_phi(1, 0.0, np.pi))
	#print(integrate_seperated_phi(1, 0.0, np.pi*0.5))
	#print(integrate_seperated_phi(1, np.pi*0.5, np.pi))