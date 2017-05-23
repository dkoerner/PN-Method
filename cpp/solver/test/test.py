import solver
import numpy as np
import matplotlib.pyplot as plt
from util import *





res = 350
domain = Domain(1.0, res)

phi = np.zeros((domain.res, domain.res))

def geometry_term( p0, n0, p1, n1 ):
	d = p1-p0
	dist = np.linalg.norm(d)

	d = d/dist

	g = 0.0
	g = 1.0/dist
	if n0 is not None:
		g*= np.abs(np.dot(n0, d))
	if n1 is not None:
		g*= np.abs(np.dot(n1, -d))
	return g

# now we compute groundtruth solution
def compute_phi_groundtruth(pWS):
	center = np.array([0.3, 0.0])
	power = 1.0

	# one bounce at reflector
	n = np.array([-1.0, 0.0])
	reflection_point = np.array([0.0, 0.0])
	#reflection_point[0] = domain.voxelToWorld(np.array([domain.res, domain.res/2.0]))[0]
	reflection_point[0] = 0.51
	reflection_point[1] = pWS[1]*0.5+center[1]*0.5
	reflection_point_ls = domain.worldToLocal(reflection_point)

	if pWS[0] > reflection_point[0]:
		return 0.0
	#print(reflection_point_ls)


	radiance = power/(2.0*np.pi)

	L0 = radiance
	L0 *= geometry_term(center, None, pWS, None)

	L1 = 0.0
	if reflection_point_ls[1] > 0.0 or reflection_point_ls[1] < 1.0:
		L1 = radiance
		# geometry term center<->reflection_point
		L1 *= geometry_term(center, None, reflection_point, n)
		L1 *= geometry_term(reflection_point, n, pWS, None)

	return L0+L1

phi_groundtruth = domain.rasterize(compute_phi_groundtruth)

#solver.toSRGB(phi_groundtruth)

fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111)
plt.imshow(phi_groundtruth, interpolation="nearest", extent = [domain.bound_min[0], domain.bound_max[0], domain.bound_min[1], domain.bound_max[1]], cmap='jet')
#plt.imshow(phi_groundtruth, interpolation="bilinear", extent = [domain.bound_min[0], domain.bound_max[0], domain.bound_min[1], domain.bound_max[1]], cmap='gray')

plt.show()


'''
n = np.array([-1.0, 0.0])


N = 10
for i in range(N):
	reflection_point = np.array([0.0, 0.0])
	reflection_point[0] randomly select a point on the reflection surface
	reflection_point[1] = 
	p0 = np.array( [np.random.uniform()*0.89, np.random.uniform()] )
	# find p1 by performing the reflection
	wi = p0-p
	wo = -wi + 2.0*n*np.dot(n, wi)
	p1 = p+wo

	# now we find the y coordinate from p0 and p1 and compare it with the correct solution
	y = p0[1]*0.5+p1[1]*0.5
	print("check: y={}  found={}".format(p[1], y))
'''

exit()






def rasterize_Q( domain, center ):
	''' rasterizes the emission field (defined by a point light) into a grid
	input: Domain object which defines resolution and voxel to world mapping
	returns: 2d grid with discretized emission field values
	'''
	# find the voxel at the domain center
	#pWS = (domain.bound_max+domain.bound_min)/2
	pWS = center
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

# reflection matrix...
def reflection_matrix( n ):
	#n = np.array([-1.0, 0.0])
	return np.identity(2) - 2.0*np.outer(n, n)


# --------------------------------------------------
#pointsource_center = np.array([0.4, 0.0])
pointsource_center = np.array([0.0, 0.0])
res = 161
#res = 2 
domain = Domain(1.0, res)
#domain = Domain(1.0, 10)


#N = 0
#N = 1
#N = 10
#N = 100
#N = 250
#N = 300
#N = 325
#N = 350
#N = 400
#N = 500
#N = 1000
#N = 2000
#N = 5000
#N = 12500
N = 20500
#N = 3500





# setup grids ---
phi = np.zeros((domain.res, domain.res)) # scalar field for which we want to solve 
Q = np.zeros((domain.res, domain.res)) # emission field
Q = rasterize_Q(domain, pointsource_center) # source term/emission field
D = np.zeros((domain.res, domain.res)) # diffusion coefficients
phi_boundary = np.zeros((domain.res, domain.res)) # phi boundary values for dirichlet boundary condition
D_boundary = np.zeros((domain.res, domain.res)) # D boundary values for dirichlet boundary condition 


# In each voxel we store two matrices. One at the left and another at the bottom cell face (staggered grid).
# M_x holds all matrices on the cell faces in x direction and M_y in y direction
M_x = np.zeros((domain.res+1, domain.res+1, 2, 2))
M_y = np.zeros((domain.res+1, domain.res+1, 2, 2))
# initialize each entry in M to be the identity matrix
for i in range(0, domain.res+1):
	for j in range(0, domain.res+1):
		M_x[i, j] = np.identity(2)
		M_y[i, j] = np.identity(2)
		'''
		if i==80 and j == 80:
			M_x[i, j] = np.array([[0.0, 1.0], [2.0, 3.0]])
			M_y[i, j] = np.array([[0.0, 1.0], [2, 3.0]])
		if i==80 and j == 81:
			M_x[i, j] = np.array([[3.0, 2.0], [1.0, 0.0]])
		if i==81 and j == 80:
			M_y[i, j] = np.array([[3.0, 2.0], [1.0, 0.0]])
		'''

# now lets get crazy...initialize the rightmost anisotropy matrices to reflection matrices
R = reflection_matrix(np.array([-1.0, 0.0]))
for i in range(domain.res+1):
	#M_x[i, domain.res] = R
	#M_x[i, domain.res-1] = R
	pass


# compute groundtruth results for evaluation and boundary conditions ---
phi_groundtruth = domain.rasterize( lambda pWS: pointsource_2d(pWS, pointsource_center) )
gradphi_groundtruth = domain.gradient(phi_groundtruth)
D_groundtruth = domain.rasterizeVS( lambda i,j:compute_D(phi_groundtruth[i,j], gradphi_groundtruth[i,j]) )

for i in range(domain.res):
	phi_boundary[i, domain.res-1] = 10.0


# initialize boundary values from groundtruth results ---
#D_boundary = np.copy(D_groundtruth)
#phi_boundary = np.copy(phi_groundtruth)


# solve ---
for step in range(N):
	# here we call c++ code for performance reasons
	debug = False
	if step>20400:
		debug = True
	solver.iterate_2d(phi, D, Q, phi_boundary, D_boundary, domain.h, debug)
	#solver.iterate_2d_anisotropic(phi, D, Q, M_x, M_y, phi_boundary, D_boundary, domain.h)

phi[76, 93] = 10.0
'''
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
'''


# plot phi 2d ------------------------------------------------
fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111)
ax.set_title('phi')
plt.imshow(phi, interpolation="nearest", extent = [domain.bound_min[0], domain.bound_max[0], domain.bound_min[1], domain.bound_max[1]])
plt.show()

# plot 1d slice ------------------------------------------------
center_voxel = domain.worldToVoxel((domain.bound_max+domain.bound_min)/2.0)
center_voxel_x = int(center_voxel[0])
center_voxel_y = int(center_voxel[1])
print("center_voxel={} {}".format(center_voxel_x, center_voxel_y))
center_voxel = domain.worldToVoxel((domain.bound_max+domain.bound_min)/2.0)
domain_x = domain.rasterize( lambda pWS: pWS[0] )

phi_groundtruth[center_voxel_y, center_voxel_x+15] = 100
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

