import numpy as np
import renderer
import util
import json
import random

import render

from scipy.special import sph_harm

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



def csp( m ):
	if m % 2 == 0:
		return 1.0
	else:
		return -1.0


def sph_real( l, m, theta, phi ):
	# https://en.wikipedia.org/wiki/Spherical_harmonics#Real_form
	# NB: We dont include the condon-shortley-phase as it seems to be already part of sph_harm.
	if m < 0:
		return np.sqrt(2.0)*np.imag(sph_harm(np.abs(m), l, phi, theta))
	elif m == 0:
		return np.real(sph_harm(0, l, phi, theta))
	elif m > 0:
		return np.sqrt(2.0)*np.real(sph_harm(m, l, phi, theta))


def sph_real2( l, m, theta, phi ):
	# https://en.wikipedia.org/wiki/Spherical_harmonics#Real_form
	# NB: We dont include the condon-shortley-phase as it seems to be already part of sph_harm.
	if m < 0:
		return np.real(np.complex(0.0, 1.0)/np.sqrt(2.0)*( csp(m)*sph_harm(-np.abs(m), l, phi, theta)-sph_harm(np.abs(m), l, phi, theta) ))
	elif m == 0:
		return np.real(sph_harm(0, l, phi, theta))
	elif m > 0:
		#return np.sqrt(2.0)*csp(m)*np.real(sph_harm(m, l, phi, theta))
		return np.real(1.0/np.sqrt(2.0)*( csp(m)*sph_harm(-np.abs(m), l, phi, theta)+sph_harm(np.abs(m), l, phi, theta) ))


def sph_real3( l, m, theta, phi ):
	# https://en.wikipedia.org/wiki/Spherical_harmonics#Real_form
	# NB: We dont include the condon-shortley-phase as it seems to be already part of sph_harm.
	if m < 0:
		return np.real(np.complex(0.0, 1.0)/np.sqrt(2.0)*( csp(m)*sph_harm(m, l, phi, theta)-sph_harm(-m, l, phi, theta) ))
	elif m == 0:
		return np.real(sph_harm(0, l, phi, theta))
	elif m > 0:
		#return np.sqrt(2.0)*csp(m)*np.real(sph_harm(m, l, phi, theta))
		return np.real(1.0/np.sqrt(2.0)*( csp(m)*sph_harm(-m, l, phi, theta)+sph_harm(m, l, phi, theta) ))


def clamped_sh( pns, theta, phi ):
    # NB: the negative direction...
    value = pns.eval(pWS, -util.sphericalDirection(theta, phi))
    if value < 0.0:
        value = 0.0
    return value

def plot_spherical_function( fun, vmin=-0.5, vmax=0.5 ):
    # plot sh functions
    theta = np.arange(0.0, 1.0, 0.01)*np.pi
    phi = np.arange(0.0, 1.0, 0.01)*2.0*np.pi

    f_img = np.zeros((theta.shape[0], phi.shape[0]))
    for j in range(phi.shape[0]):
        for i in range(theta.shape[0]):
            f_img[i, j] = fun(theta[i], phi[j])
    #return plt.imshow(f_img, interpolation="nearest", cmap='jet', vmin=vmin, vmax=vmax)
    #vmin = np.min(f_img)
    #vmax = np.max(f_img)
    #print("vmin={} vmax={}".format(vmin, vmax))
    return f_img
    #return plt.imshow(f_img, origin='lower',zorder=1, interpolation="nearest", cmap='jet', vmin=vmin, vmax=vmax, extent=[0.0,2.0*np.pi,0.0,np.pi])


if __name__ == "__main__":

	'''
	theta = np.pi*0.4
	phi = 2.0*np.pi*0.4

	order = 1
	for l in range(order+1):
		for m in range(-l,l+1):
			#print("l={} m={} sph_real={}".format(l, m, sph_real(l, m, theta, phi)))
			print("l={} m={} sph_complex_basis={}".format(l, m, sph_harm(m, l, phi, theta)))
			#a = np.conj(sph_harm(m, l, phi, theta))
			#b = csp(m)*sph_harm(-m, l, phi, theta)
			#print("l={} m={} a={} b={}".format(l, m, a, b))
	print("\n\n")
	'''



	#l = 1
	#m = -1

	#print("sph_complex={}".format(sph_harm(m, l, phi, theta)))

	

	#c = np.complex(0.0, 1.0)/np.sqrt(2.0)*csp(m)*sph_harm(m, l, phi, theta)
	#-np.complex(0.0, 1.0)/np.sqrt(2.0)*csp(m)*sph_harm(-m, l, phi, theta)

	#exit(1)
	#print(c)

	
	

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

	#img_groundtruth = renderer.load_image("test.exr")
	#img_groundtruth = renderer.blur_image(img_groundtruth, 3.0)
	#img_groundtruth_sampler = renderer.create_image_sampler(img_groundtruth)



	
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

	pns = renderer.load_pnsolution( "results/pointsource/pointsource_p5.pns" )
	#pWS = np.array([1.1, 1.0, 1.0])
	for i in range(14):
		random.seed(1233*i)
		pWS = np.array([random.random()*2.0, random.random()*2.0, random.random()*2.0])

		# visualize samples and the sampled sh function ------------------
		samples_theta = []
		samples_phi = []

		numSamples = 10000
		for i in range(numSamples):
			d = -pns.sample(pWS, np.array([random.random(), random.random()]))
			theta, phi = util.sphericalCoordinates(d)
			samples_theta.append(theta)
			samples_phi.append(phi)
			#pass



		fig = plt.figure(figsize=(8,8));
		ax = fig.add_subplot(121)

		img = plot_spherical_function( lambda theta, phi: clamped_sh(pns, theta, phi) )
		img_view = plt.imshow(img, origin='lower',zorder=1, interpolation="nearest", cmap='jet', vmin=0.0, vmax=np.max(img), extent=[0.0,2.0*np.pi,0.0,np.pi])
		plt.scatter(samples_phi, samples_theta, c='r', label="sdda",zorder=2, marker='.', s=1.0)
		#divider = make_axes_locatable(ax)
		#cax = divider.append_axes("right", size="5%", pad=0.05)
		#plt.colorbar(img_view, cax=cax)

		ax = fig.add_subplot(122)

		img = plot_spherical_function( lambda theta, phi: pns.pdf(pWS, -util.sphericalDirection(theta, phi)) )
		img_view = plt.imshow(img, origin='lower',zorder=1, interpolation="nearest", cmap='jet', vmin=0.0, vmax=np.max(img), extent=[0.0,2.0*np.pi,0.0,np.pi])
		#img_view = plt.imshow(pns.getBlocks(pWS, 3), origin='lower',zorder=1, interpolation="nearest", cmap='jet', vmin=0.0, vmax=np.max(img), extent=[0.0,2.0*np.pi,0.0,np.pi])
		#plt.scatter(samples_phi, samples_theta, c='r', label="sdda",zorder=2, marker='.', s=1.0)
		#divider = make_axes_locatable(ax)
		#cax = divider.append_axes("right", size="5%", pad=0.05)
		#plt.colorbar(img_view, cax=cax)



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