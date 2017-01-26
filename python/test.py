import pybind_test as pb
import numpy as np
import matplotlib.pyplot as plt
from util import Domain


def pointsource_2d( pWS ):
	r = np.linalg.norm(pWS)
	if r == 0.0:
		return 1.0/(2.0*np.pi*max(r, 1.0e-4))
	return 1.0/(2.0*np.pi*r)

def rasterize_Q( domain ):
	''' rasterizes the emission field (defined by a point light) into a grid
	input: Domain object which defines resolution and voxel to world mapping
	returns: 2d grid with discretized emission field values
	'''
	# find the voxel at the domain center
	pWS = (domain.bound_max+domain.bound_min)/2
	pVS = domain.worldToVoxel(pWS)

	# initialize the result
	q = np.zeros((domain.res, domain.res))

	# we approximate the pointlight with a voxel and therefore
	# we distribute the pointlight unit energy over the voxelarea
	q[int(pVS[1]), int(pVS[0])] = 1.0/(domain.h*domain.h)

	return q

def compute_D( phi, gradphi ):
	''' this function computes the diffusion coefficient from given phi
	and phi gradient values
	'''
	l = np.linalg.norm(gradphi)
	if l>0.0:
		return phi/l
	return 0.0



# --------------------------------------------------

domain = Domain(1.0, 161)
#domain = Domain(1.0, 10)


#N = 0
#N = 1
#N = 600
#N = 12500
N = 20500
#N = 3500





# setup grids ---
phi = np.zeros((domain.res, domain.res)) # scalar field for which we want to solve 
Q = rasterize_Q(domain) # source term/emission field
D = np.zeros((domain.res, domain.res)) # diffusion coefficients
phi_boundary = np.zeros((domain.res, domain.res)) # phi boundary values for dirichlet boundary condition
D_boundary = np.zeros((domain.res, domain.res)) # D boundary values for dirichlet boundary condition 

# compute groundtruth results for evaluation and boundary conditions ---
phi_groundtruth = domain.rasterize( pointsource_2d )
gradphi_groundtruth = domain.gradient(phi_groundtruth)
D_groundtruth = domain.rasterizeVS( lambda i,j:compute_D(phi_groundtruth[i,j], gradphi_groundtruth[i,j]) )


# initialize boundary values from groundtruth results ---
D_boundary = np.copy(D_groundtruth)
phi_boundary = np.copy(phi_groundtruth)


# solve ---
for step in range(N):
	# here we call c++ code for performance reasons
	pb.iterate(phi, D, Q, phi_boundary, D_boundary, domain.h)


# compute the residual to evaluate the solution ---
residual = np.zeros_like(phi)
for i in range(1, domain.res-1):
	for j in range(1, domain.res-1):
		# compute phi gradient at cell center ---
		grad_phi = np.zeros(2)
		grad_phi[0] = (phi[i, j+1]-phi[i, j-1])/(2.0*domain.h)
		grad_phi[1] = (phi[i-1, j]-phi[i+1, j])/(2.0*domain.h)

		# get center diffusion coefficient ---
		dc = D[i, j]

		# compute diffusion coefficients at cell faces ---
		D_xph = (dc + D[i, j+1])*0.5;
		D_xmh = (dc + D[i, j-1])*0.5;
		D_yph = (dc + D[i-1, j])*0.5;
		D_ymh = (dc + D[i+1, j])*0.5;

		# LHS
		lhs = domain.h*Q[i, j]

		# RHS
		rhs = 0.0
		rhs_a = (D_xph + D_xmh + D_yph + D_ymh)*phi[i,j]/(domain.h)
		rhs_b = (D_xph*phi[i, j+1]+D_xmh*phi[i, j-1]+D_yph*phi[i-1, j]+D_ymh*phi[i+1, j])/(domain.h)
		rhs -= (D_xph + D_xmh + D_yph + D_ymh)*phi[i,j]/(domain.h)
		rhs += (D_xph*phi[i, j+1]+D_xmh*phi[i, j-1]+D_yph*phi[i-1, j]+D_ymh*phi[i+1, j])/(domain.h)

		# residual
		residual[i,j] = abs(rhs-lhs)



# plot phi 2d ------------------------------------------------
fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111)
ax.set_title('phi')
plt.imshow(phi, extent = [domain.bound_min[0], domain.bound_max[0], domain.bound_min[1], domain.bound_max[1]])
plt.show()

# plot 1d slice ------------------------------------------------
center_voxel = domain.worldToVoxel((domain.bound_max+domain.bound_min)/2.0)
center_voxel_x = int(center_voxel[0])
center_voxel_y = int(center_voxel[1])
print("center_voxel={} {}".format(center_voxel_x, center_voxel_y))
domain_x = domain.rasterize( lambda pWS: pWS[0] )


plt.figure()

ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='both', direction='out')
ax.get_xaxis().tick_bottom()   # remove unneeded ticks 
ax.get_yaxis().tick_left()
plt.ylim([10.0e-2,10.0e1])
plt.loglog( domain_x[center_voxel_y, center_voxel_x:], phi[center_voxel_y, center_voxel_x:], label="solution", color = 'g', linestyle=' ', marker='.', markersize=10 )
plt.loglog( domain_x[center_voxel_y, center_voxel_x:], phi_groundtruth[center_voxel_y, center_voxel_x:], label="groundtruth", color = 'g' )
#plt.semilogy( domain[res_half:], residual[res_half, res_half:], label="" )
plt.title("point light in vacuum (2d)")
plt.xlabel(r'$r \left[m\right]$', fontsize=18)
plt.ylabel(r'$\phi \left[\frac{W}{m}\right]$', fontsize=18)
plt.grid(True, linestyle='-',color='0.75')
plt.legend(loc='best')
plt.draw()
plt.show()


#gradphi = domain.gradient(phi)
#plt.figure()
#ax = plt.gca()

#plt.plot( domain_x[center_voxel_y, center_voxel_x:], phi[center_voxel_y, center_voxel_x:], label="", color = 'g', linestyle=' ', marker='.', markersize=10 )
#plt.plot( domain_x[center_voxel_y, center_voxel_x:], phi_groundtruth[center_voxel_y, center_voxel_x:], label="", color = 'g' )


#plt.plot( domain_x[center_voxel_y, center_voxel_x:], gradphi[center_voxel_y, center_voxel_x:, 0], label="", color = 'g', linestyle=' ', marker='.', markersize=10 )
#plt.plot( domain_x[center_voxel_y, center_voxel_x:], gradphi_groundtruth[center_voxel_y, center_voxel_x:, 0], label="", color = 'g' )

#plt.plot( domain_x[center_voxel_y, center_voxel_x:], D[center_voxel_y, center_voxel_x:], label="", color = 'g', linestyle=' ', marker='.', markersize=10 )
#plt.plot( domain_x[center_voxel_y, center_voxel_x:], D_groundtruth[center_voxel_y, center_voxel_x:], label="", color = 'g' )

#plt.semilogy( domain_x[center_voxel_y, center_voxel_x:], residual[center_voxel_y, center_voxel_x:], label="" )

#plt.plot( domain_x[center_voxel_y, center_voxel_x:], residual[center_voxel_y, center_voxel_x:], label="" )

#plt.title("gradient phi")
#plt.legend(loc='best')
#plt.draw()
#plt.show()

#fig = plt.figure(figsize=(6, 6))
#ax = fig.add_subplot(111)
#ax.set_title('Q')
##ax.set_ylim([-domain_extent, domain_extent])
##ax.set_xlim([-domain_extent, domain_extent])
#plt.imshow(Q, extent = [-domain_extent, domain_extent, -domain_extent, domain_extent])
#plt.show()

