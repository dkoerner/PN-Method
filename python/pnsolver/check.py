import os,sys,inspect
import numpy as np
import pnsolver
#import renderer
import util
import stencil

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable





'''
pns = renderer.load_pnsolution( "c:/projects/epfl/epfl17/python/pnsolver/results/nebulae/nebulae_p5_2_ms.pns" )

pWS = np.array([-0.27704304, 0.36083166, -0.22953043])
#sample = np.array([0.3960667149031164, 0.18571028805648737])
#sample = np.array([0.8444218515250481, 0.7579544029403025])
sample = np.array([0.37858849022321517, 0.49207364263387265])

dir, pdf = pns.sample(pWS, sample)
print("")
pdf2 = pns.pdf( pWS, dir )

print(  "pdf={:.02} pdf2={:.02}".format(pdf, pdf2))
'''


'''
a, b, c = pnsolver.test()


fig = plt.figure(figsize=(15,7));
ax = fig.add_subplot(111)
#img_view = ax.imshow(img[:, :, 0].T, cmap='jet', norm=LogNorm(vmin=np.min(img), vmax=np.max(img)), origin='lower')
plt.plot( a, b, label="groundtruth" )
plt.plot( a, c, linestyle=' ', marker='.', label="reconstruction" )
plt.legend(loc='best')
plt.show()
'''

pnsolver.test()


'''
pns = pnsolver.load_solution("C:/projects/epfl/epfl17/python/pnsolver/results/checkerboard2d/checkerboard2d_p1.pns")
test = pns.getCoefficientField(0)
print(test.shape)


img = np.clip( test, 1.0e-8, np.max(test) )

fig = plt.figure(figsize=(15,7));
ax = fig.add_subplot(111)
img_view = ax.imshow(img[:, :, 0].T, cmap='jet', norm=LogNorm(vmin=np.min(img), vmax=np.max(img)), origin='lower')
plt.show()
'''


'''
import stencil


order = 2
pni = stencil.PNInfo3D(order)

S_org = pni.getS()
S_inv_org = pni.getSInv()
S_ours = pnsolver.createComplexToRealConversionMatrix(order)
S_inv_ours = pnsolver.createRealToComplexConversionMatrix(order)


for i in range(pni.num_coeffs()):
	for j in range(pni.num_coeffs()):
		#print(np.abs(np.real(S_org[i, j])-np.real(S_ours[i, j])))
		#print(np.abs(np.imag(S_org[i, j])-np.imag(S_ours[i, j])))
		#print(np.abs(np.real(S_inv_org[i, j])-np.real(S_inv_ours[i, j])))
		#print(np.abs(np.imag(S_inv_org[i, j])-np.imag(S_inv_ours[i, j])))
'''






#filename = "C:/projects/epfl/epfl17/python/pnsolver/test.pns"
#pns = pnsolver.load_solution(filename)
#numCoeffs = pns.getNumCoeffs()

#x = pnsolver.getSolutionVector(pns)
#x = x.reshape((x.shape[0], 1))

'''
for coeff in range(numCoeffs):
    image = util.extract_coefficient_field( x, pns.getResolution(), numCoeffs, coeff )
    image = np.clip(image, 1.0e-8, np.max(image))


    fig = plt.figure(figsize=(8,8));
    ax = fig.add_subplot(111)
    #plt.title("Solution of LSPN0, which just has a single cofficient")
    img_view = ax.imshow(image[:,:,20], cmap='jet', norm=LogNorm(vmin=np.min(image), vmax=np.max(image)), origin='lower')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(img_view, cax=cax)

    plt.show()


'''

'''
def plot_spherical_function( fun, vmin=-0.5, vmax=0.5 ):
    # plot sh functions
    theta = np.arange(0.0, 1.0, 0.01)*np.pi
    phi = np.arange(0.0, 1.0, 0.01)*2.0*np.pi

    f_img = np.zeros((theta.shape[0], phi.shape[0]))
    for j in range(phi.shape[0]):
        for i in range(theta.shape[0]):
            f_img[i, j] = fun(theta[i], phi[j])
    #return plt.imshow(f_img, interpolation="nearest", cmap='jet', vmin=vmin, vmax=vmax)
    return plt.imshow(f_img, interpolation="nearest", cmap='jet', vmin=np.min(f_img), vmax=np.max(f_img))

def fun( theta, phi ):
    dir = util.sphericalDirection(theta, phi)
    pLS = np.array([0.7, 0.5, 0.0])
    pWS = pns.localToWorld(pLS)
    return pns.eval(pWS, dir)
    
fig = plt.figure(figsize=(8,8));
ax = fig.add_subplot(111)

img_view = plot_spherical_function( fun, vmin=-0.5, vmax=0.5 )
    
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(img_view, cax=cax)

plt.show()
'''