import os,sys,inspect
import numpy as np
import pnsolver
import util
import stencil

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable


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






filename = "C:/projects/epfl/epfl17/python/pnsolver/test.pns"
pns = pnsolver.load_solution(filename)
numCoeffs = pns.getNumCoeffs()

x = pnsolver.getSolutionVector(pns)
x = x.reshape((x.shape[0], 1))
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