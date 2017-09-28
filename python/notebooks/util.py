import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

from scipy.special import sph_harm, lpmv


# Spherical Harmonics related =====================================================

def num_sh_coeffs(order):
    return (order + 1) * (order + 1)

def sh_index( l, m):
	if l<0 or m < -l or m > l:
		return None
	return l * (l + 1) + m

def sh_sum( order, coeffs, theta, phi ):
    result = 0.0
    for l in range(order+1):
        for m in range(-l, l+1):
            result+=coeffs[sh_index(l,m)]*sph_harm(m, l, phi, theta)
    return result



class PNInfo2D(object):
	def __init__(self, order, staggered=True):
		self.order = order

		# the following dictionaries help mapping lm to shindex and vice versa ---
		# NB: sh_index refers to the complex valued SH coefficient
		# find out how much coefficients we have in 2d
		self.index_to_lm = [] # this will map equation index to l,m indices
		self.lm_to_index = {} # this dict will map l,m indices to the equation index
		# NB: we dont use l(l+1)+m because in 2d, we skipp odd (l+m) entries, which changes the sequence.
		# iterate all sh coefficients for given truncation order
		for l in range(0, self.order+1):
			for m in range(-l, l+1):
				# in 2d, we only need to solve for moments where l+m is even
				if (l+m) % 2 == 0:
					self.index_to_lm.append( (l, m) )
					self.lm_to_index[(l,m)] = len(self.index_to_lm)-1
		self.numCoeffs = len(self.index_to_lm)

		self.build_S()

		# the following dictionaries help mapping lm to uindex and vice versa ---
		# NB: u_index refers to the real-valued coefficient of the solution vector u which
		# represents real and imaginary part (both are real numbers) of sh coefficients lm with m>=0
		self.u_index_to_lm = []
		self.lm_to_u_index = {}
		local_u_index = 0 # this counter represents the current row in the solution vector for a single voxel
		for l in range(0, self.order+1):
			# NB: the starmap code builds the solution vector from reverse order, we do it the same,
			# as we want our conversion matrix S (and therefore the structure of Mx, My) to be indentical
			for m in range(l, -1, -1):

				# In 2d, we only need to solve for moments where l+m is even.
				# This is why all other coefficients are not even part of the global system Ax=b.
				if (l+m) % 2 != 0:
					continue

				# handle real valued part of the current sh coefficient
				key = (l,m,0)
				self.u_index_to_lm.append(key)
				self.lm_to_u_index[key] = local_u_index
				local_u_index += 1

				# handle imaginary part of the current sh coefficient
				# NB: if m==0, then the imaginary part will always be zero, this is why we skip it
				if m>0:
					key = (l,m,1)
					self.u_index_to_lm.append(key)
					self.lm_to_u_index[key] = local_u_index
					local_u_index += 1


		# we keep track of where the unknowns are placed
		self.unknown_info = [ {} for i in range(self.numCoeffs)]

		# by default we place all unknowns at the cell centers
		for i in range(self.numCoeffs):
			self.place_unknown(i, (1, 1))
		# and we use full voxel central differences
		self.stencil_half_steps = 2

		if staggered == True:
			# TODO: generalize this to arbitrary order
			# see starmap python implementation on how to do it
			# but how to do it without M matrix?
			self.place_unknown( 0, (1,1) )
			if self.numCoeffs > 1:
				self.place_unknown( 1, (0,1) )
			if self.numCoeffs > 2:
				self.place_unknown( 2, (1,0) )
			self.stencil_half_steps = 1


	def build_S(self):
		'''builds the S matrix, which converts from complex-valued to real valued coefficients'''
		# build S matrix ( we iterate over l, m to make sure that the order is correct)

		self.S = np.zeros((self.numCoeffs, self.numCoeffs),dtype=complex)
		count = 0
		for l in range(0, self.order+1):
			for m in range(l, -1, -1):
				# in 2D, we skip coefficients for which l+m is odd
				if (l+m) % 2 != 0:
					continue
				
				# build S matrix, which converts solution from complex to real values

				# computes the real part coefficients for a row (defined by l,m) in the S matrix
				# (see bottom of p.5 in the starmap paper)
				if m == 0:
					self.S[count, self.lm_to_index[(l,m)]] = 1.0
				else:
					self.S[count, self.lm_to_index[(l,m)]] = np.power(-1.0, m)/np.sqrt(2)
					if (l,-m) in self.lm_to_index:
						self.S[count, self.lm_to_index[(l,-m)]] = np.power(-1.0, 2.0*m)/np.sqrt(2)
				count+=1

				# computes the imaginary part coefficients for a row (defined by l,m) in the S matrix
				# (see bottom of p.5 in the starmap paper)
				if m > 0:
					self.S[count, self.lm_to_index[(l,m)]] = np.power(-1.0, m)/np.sqrt(2)*1j
					if (l,-m) in self.lm_to_index:
						self.S[count, self.lm_to_index[(l,-m)]] = -np.power(-1.0, 2*m)/np.sqrt(2)*1j
					count+=1
		self.S_inv = np.linalg.inv(self.S)

	def to_complex( self, x_real):
	    # use this to convert the solution from complex valued to real valued
	    numVoxels = int(x_real.shape[0]/self.numCoeffs)
	    x_complex = np.zeros( (numVoxels*self.numCoeffs, 1), dtype=complex )
	    for i in range(numVoxels):
	        block_i = i*self.numCoeffs
	        x_complex[block_i:block_i + self.numCoeffs, :] = np.real(self.S_inv @ x_real[block_i:block_i + self.numCoeffs, :])

	    return x_complex

	def to_real(self, x_complex):
		raise ValueError("not correctly implemented yet")
		# use this to convert the solution from real valued to complex valued
		numVoxels = self.domain.res_x*self.domain.res_y
		x_real = np.zeros( (numVoxels*self.numCoeffs), dtype=float )
		for voxel_x in range(self.domain.res_x):
			for voxel_y in range(self.domain.res_y):
				block_i = self.get_global_index(voxel_x, voxel_y, 0)
				#if block_i <= 6628 and block_i+self.numCoeffs >= 6628:
				#if block_i == 6627:
				#	print(block_i)
				#	print(x_complex[block_i:block_i + self.numCoeffs])
				#	print(np.real(self.S.dot(x_complex[block_i:block_i + self.numCoeffs])))
				x_real[block_i:block_i + self.numCoeffs] = np.real(self.S.dot(x_complex[block_i:block_i + self.numCoeffs]))
		return x_real

	def getS(self):
		return self.S

	def getSInv(self):
		return self.S_inv

	def num_coeffs(self):
		return self.numCoeffs

	def u_index( self, l, m, part ):
		# the third argument identifies the real(0) or imaginary(1) part of the complex number
		key = (l,m,part)
		if key in self.lm_to_u_index:
			return self.lm_to_u_index[key]
		return None

	def lmp_index(self, coeff_index):
		return self.u_index_to_lm[coeff_index]

	def sh_index(self, l, m):
		key = (l,m)
		#NB: we dont use l(l+1)+m because in 2d, we skipp odd (l+m) entries, which changes the sequence.
		if key in self.lm_to_index:
			return self.lm_to_index[key]
		return None

	def lm_index(self, coeff_index):
		return self.index_to_lm[coeff_index]

	def getOrder(self):
		return self.order


	def place_unknown( self, coeff_index, grid_id ):
		self.unknown_info[coeff_index]['grid_id'] = grid_id
		self.unknown_info[coeff_index]['offset'] = np.array( [grid_id[0], grid_id[1]] , dtype=int)

	def unknown_offset(self, coeff_index):
		return self.unknown_info[coeff_index]['offset']

	def getOffset(self, coeff_index):
		return self.unknown_info[coeff_index]['offset']

	#def getLocation(self, voxel, coeff_index):
	#	return GridLocation2D(voxel, self.unknown_offset(coeff_index)) 



def load_pn_system( filename ):
    #print("loading PN solution from {}".format(filename))
    data = scipy.io.loadmat(filename)
   
    result = {}
    if "sigma_t" in data:
        result["sigma_t"] = data["sigma_t"]
    if "sigma_a" in data:
        result["sigma_a"] = data["sigma_a"]
    if "sigma_s" in data:
        result["sigma_s"] = data["sigma_s"]
    if "q00" in data:
        result["q00"] = data["q00"]
    if "x" in data:
        result["x"] = data["x"]
    if "b" in data:
        result["b"] = data["b"]
    if "A" in data:
        result["A"] = data["A"]
    if "info" in data:
        info = data["info"]
        result["order"] = info["order"][0][0][0][0]
        result["numCoeffs"] = info["numCoeffs"][0][0][0][0]
        result["resolution"] = np.array(info["resolution"][0][0][0])
    else:
        result["order"] = 1
        result["numCoeffs"] = 3
        result["resolution"] = np.array([70, 70])
        
    #print("\torder={}  numCoeffs={}  resolution={} {}".format(result["order"], result["numCoeffs"], result["resolution"][0], result["resolution"][1]))
        
    return result


   
def visualize_solution_vector( x, res, numCoeffs ):

    coeff = 0 # the coefficient we want to visualize
    res_x = res[0]
    res_y = res[1]
    
    u0 = np.zeros( (res_x, res_y) )
    for voxel_i in range(res_x):
        for voxel_j in range(res_y):
            i = (voxel_j*res_x + voxel_i)*numCoeffs + coeff
            value = x[i, 0]
            u0[voxel_i, voxel_j] = np.real(value)

    #u0 = np.abs(u0)
    #u0 = -u0
    u0 = np.clip(u0,1.0e-8, np.max(u0))

    vmin = np.min(u0)
    vmax = np.max(u0)
    
    if vmin >= vmax:
        vmin = vmax
    #print("vmin={} vmax={}".format(vmin, vmax))

    if vmin==vmax or vmin < 0.0:
        img_view = plt.imshow(u0.T, interpolation="nearest", cmap='jet', vmin=vmin, vmax=vmax, origin='lower')
    else:
        img_view = plt.imshow(u0.T, interpolation="nearest", cmap='jet', norm=LogNorm(vmin=vmin, vmax=vmax), origin='lower')
    return img_view



def visualize_solution(result):
	fig = plt.figure(figsize=(15, 15));
	ax = plt.subplot(131)
	plt.title("x>0")

	img_view = visualize_solution_vector( result["x"], result["resolution"], result["numCoeffs"] )
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	plt.colorbar(img_view, cax=cax)

	ax = plt.subplot(132)
	plt.title("-x>0")
	img_view = visualize_solution_vector( -result["x"], result["resolution"], result["numCoeffs"] )
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	plt.colorbar(img_view, cax=cax)

	ax = plt.subplot(133)
	plt.title("abs(x)")
	img_view = visualize_solution_vector( np.abs(result["x"]), result["resolution"], result["numCoeffs"] )
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	plt.colorbar(img_view, cax=cax)

	plt.show()


def visualize_solution_simple(result):
	fig = plt.figure(figsize=(6, 6));
	ax = plt.subplot(111)
	img_view = visualize_solution_vector( np.abs(result["x"]), result["resolution"], result["numCoeffs"] )
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	plt.colorbar(img_view, cax=cax)

	plt.show()

def pointsource_2d_fluence( pWS, center = np.array([0.0, 0.0]), power = 1.0 ):
	''' returns fluence at position pWS for a unit power point light located at origin'''
	r = np.linalg.norm(pWS-center)
	if r == 0.0:
		return power/(2.0*np.pi*max(r, 1.0e-4))
	return power/(2.0*np.pi*r)

def geometry_term_2d( p0, n0, p1, n1 ):
    '''Geometry term in 2d'''
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


def integrate( func, numSamples ):
	tmin = 0.0
	tmax = 2.0*np.pi
	P = tmax-tmin
	dt = P/numSamples
	s = 0.0
	for i in range(numSamples):
		t = dt*i
		s += func(t)*dt
	return s

def normalize( v ):
	length = np.linalg.norm(v)
	if length>0.0:
		return v/length

class SGGX2D:
	#def __init__(self, S_xx = None, S_xy = None, S_yy = None, ):
	#	self.S_xx = S_xx;
	#	self.S_xy = S_xy;
	#	self.S_yy = S_yy;

	def __init__(self, frame_s, projectedDistances ):

		frame_t = np.array([-frame_s[1], frame_s[0]])

		S_11 = projectedDistances[0]*projectedDistances[0]
		S_22 = projectedDistances[1]*projectedDistances[1]

		self.S_xx = S_11*frame_s[0]*frame_s[0] + S_22*frame_t[0]*frame_t[0];
		self.S_xy = S_11*frame_s[0]*frame_s[1] + S_22*frame_t[0]*frame_t[1];
		self.S_yx = self.S_xy
		self.S_yy = S_11*frame_s[1]*frame_s[1] + S_22*frame_t[1]*frame_t[1];


	def projectedDistance( self, d ):
		sigma_squared = self.S_xx*d[0]*d[0] + self.S_yy*d[1]*d[1] + 2.0*self.S_xy*d[0]*d[1];
		if sigma_squared>0.0:
			return np.sqrt(sigma_squared)
		# conditional to avoid numerical errors
		return 0.0

	def projectedDistanceAdapter( self, t ):
		return self.projectedDistance( np.array([np.cos(t), np.sin(t)]) )

	def get_moment_0(self):
		# compute fourier coefficient a_0
		a_0 = integrate( self.projectedDistanceAdapter, 10000 )/(2.0*np.pi)
		return a_0

	def get_moment_2(self):
		def func_a(t):
			return self.projectedDistanceAdapter(t)*np.cos(2*t)
		def func_b(t):
			return self.projectedDistanceAdapter(t)*np.sin(2*t)
		# compute fourier coefficients a_2 and b2
		a_2 = integrate( func_a, 10000 )/np.pi
		b_2 = integrate( func_b, 10000 )/np.pi
		# turn this into moments (see notebook: moment_expansion_2d for details)
		return np.array([[a_2, b_2], [b_2, -a_2]])


class Domain3D:
	def __init__(self, size, res, center=None):
		self.res = res
		self.size = size
		self.h = size/float(res)
		if center == 'origin':
			# center of the bounding box is origin
			self.bound_min = np.array([-size*0.5, -size*0.5, -size*0.5])
			self.bound_max = np.array([size*0.5, size*0.5, size*0.5])
		else:
			# origin is at lower left
			self.bound_min = np.array([0.0, 0.0, 0.0])
			self.bound_max = np.array([size, size, size])

		#self.voxelToWorldTransform = 

	def voxelToLocal( self, pVS ):
		return np.array([pVS[0]/self.res, pVS[1]/self.res, pVS[2]/self.res])

	def localToVoxel(self, pLS):
		return np.array([pLS[0]*self.res, pLS[1]*self.res, pLS[2]*self.res])

	def localToWorld( self, pLS ):
		return np.array([pLS[0]*self.size + self.bound_min[0], pLS[1]*self.size + self.bound_min[1], pLS[2]*self.size + self.bound_min[2]])

	def worldToLocal( self, pWS ):
		return np.array([(pWS[0]-self.bound_min[0])/self.size, (pWS[1]-self.bound_min[1])/self.size, (pWS[2]-self.bound_min[2])/self.size])

	def voxelToWorld(self, pVS):
		return self.localToWorld(self.voxelToLocal(pVS))

	def worldToVoxel(self, pWS):
		return self.localToVoxel(self.worldToLocal(pWS))

	def voxelToIndex(self, pVS):
		return (int(pVS[1]), int(pVS[0]), int(pVS[2]))

	def worldToIndex(self, pWS):
		pVS = self.worldToVoxel(pWS)
		return (int(pVS[1]), int(pVS[0]), int(pVS[2]))

	def rasterize( self, fun, shape = None ):
		if shape == None:
			shape = (self.res, self.res)
		field = np.zeros(shape)
		for i in range(0, self.res):
			for j in range(0, self.res):
				for k in range(0, self.res):
					pVS = np.array([j+0.5, i+0.5, k+0.5])
					pWS = self.voxelToWorld(pVS)
					field[i, j, k] = fun(pWS)
		return field

	def rasterizeVS( self, fun ):
		field = np.zeros((self.res, self.res, self.res))
		for i in range(0, self.res):
			for j in range(0, self.res):
				for k in range(0, self.res):
					field[i, j, k] = fun(i, j, k)
		return field

	def gradient( self, field ):
		grad_field = np.zeros((self.res, self.res, self.res, 3))
		for i in range(0, self.res):
			for j in range(0, self.res):
				for k in range(0, self.res):
					if j==0:
						# forward differences for x
						grad_field[i,j,k,0] = (field[i,j+1,k]-field[i,j,k])/self.h
					elif j==self.res-1:
						#backward differences for x
						grad_field[i,j,k,0] = (field[i,j,k]-field[i,j-1,k])/self.h
					else:
						# central differences for x
						grad_field[i,j,k,0] = (field[i,j+1,k]-field[i,j-1,k])/(self.h*2.0)

					if i==0:
						# forward differences for y
						grad_field[i,j,k,1] = (field[i,j,k]-field[i+1,j,k])/self.h
					elif i==self.res-1:
						#backward differences for y
						grad_field[i,j,k,1] = (field[i-1,j,k]-field[i,j,k])/self.h
					else:
						# central differences for y
						grad_field[i,j,k,1] = (field[i-1,j,k]-field[i+1,j,k])/(self.h*2.0)

					if k==0:
						# forward differences for z
						grad_field[i,j,k,2] = (field[i,j,k]-field[i,j,k+1])/self.h
					elif k==self.res-1:
						#backward differences for z
						grad_field[i,j,k,2] = (field[i,j,k-1]-field[i,j,k])/self.h
					else:
						# central differences for z
						grad_field[i,j,k,2] = (field[ix,j,k-1]-field[i,j,k+1])/(self.h*2.0)

		return grad_field

class Domain2D:
	def __init__(self, size, res, center=None):
		self.res_x = res
		self.res_y = res
		self.size_x = size
		self.size_y = size
		self.h_x = size/float(self.res_x)
		self.h_y = size/float(self.res_y)
		# NB: voxelsize channels are switched because i component is in y-axis, while
		# j component is in x-axis
		self.voxelsize = np.array([self.h_y, self.h_x])
		self.numVoxels = self.res_x*self.res_y
		if center == 'origin':
			# center of the bounding box is origin
			self.bound_min = np.array([-self.size_x*0.5, -self.size_y*0.5])
			self.bound_max = np.array([self.size_x*0.5, self.size_y*0.5])
		else:
			# origin is at lower left
			self.bound_min = np.array([0.0, 0.0])
			self.bound_max = np.array([self.size_x, self.size_y])
		self.center = (self.bound_min + self.bound_max)*0.5

		#self.voxelToWorldTransform = 

	def voxelToLocal( self, pVS ):
		return np.array([pVS[0]/self.res_x, pVS[1]/self.res_y])

	def localToVoxel(self, pLS):
		return np.array([pLS[0]*self.res_x, pLS[1]*self.res_y])

	def localToWorld( self, pLS ):
		return np.array([pLS[0]*self.size_x + self.bound_min[0], pLS[1]*self.size_y + self.bound_min[1]])

	def worldToLocal( self, pWS ):
		return np.array([(pWS[0]-self.bound_min[0])/self.size_x, (pWS[1]-self.bound_min[1])/self.size_y])

	def voxelToWorld(self, pVS):
		return self.localToWorld(self.voxelToLocal(pVS))

	def worldToVoxel(self, pWS):
		return self.localToVoxel(self.worldToLocal(pWS))

	def voxelToIndex(self, pVS):
		return (int(pVS[1]), int(pVS[0]))

	def worldToIndex(self, pWS):
		pVS = self.worldToVoxel(pWS)
		return (int(pVS[1]), int(pVS[0]))

	'''
	def rasterize( self, fun, shape = None ):
		if shape == None:
			shape = (self.res_x, self.res_y)
		field = np.zeros(shape)
		for i in range(0, self.res_y):
			for j in range(0, self.res_x):
				pVS = np.array([j+0.5, i+0.5])
				pWS = self.voxelToWorld(pVS)
				field[i, j] = fun(pWS)
		return field
	'''

	def rasterize( self, fun, shape = None ):
		if shape == None:
			shape = (self.res_x, self.res_y)
		field = np.zeros(shape)
		for i in range(0, self.res_y):
			for j in range(0, self.res_x):
				pVS = np.array([i+0.5, j+0.5])
				pWS = self.voxelToWorld(pVS)
				field[i, j] = fun(pWS)
		return field

	def rasterizeVS( self, fun ):
		field = np.zeros((self.res_x, self.res_y))
		for i in range(0, self.res_y):
			for j in range(0, self.res_x):
				field[i, j] = fun(i, j)
		return field

	def gradient( self, field ):
		grad_field = np.zeros((self.res_x, self.res_y, 2))
		for i in range(0, self.res_y):
			for j in range(0, self.res_x):
				if j==0:
					# forward differences for x
					grad_field[i,j,0] = (field[i,j+1]-field[i,j])/self.h_x
				elif j==self.res-1:
					#backward differences for x
					grad_field[i,j,0] = (field[i,j]-field[i,j-1])/self.h_x
				else:
					# central differences for x
					grad_field[i,j,0] = (field[i,j+1]-field[i,j-1])/(self.h_x*2.0)

				if i==0:
					# forward differences for y
					grad_field[i,j,1] = (field[i,j]-field[i+1,j])/self.h_y
				elif i==self.res-1:
					#backward differences for y
					grad_field[i,j,1] = (field[i-1,j]-field[i,j])/self.h_y
				else:
					# central differences for y
					grad_field[i,j,1] = (field[i-1,j]-field[i+1,j])/(self.h_y*2.0)
		return grad_field


class Domain1D:
	def __init__(self, size, res):
		self.res = res
		self.size = size
		self.h = size/float(res)
		self.bound_min = -size*0.5
		self.bound_max = size*0.5

	def voxelToLocal( self, pVS ):
		return pVS/self.res

	def localToVoxel(self, pLS):
		return pLS*self.res

	def localToWorld( self, pLS ):
		return pLS*self.size + self.bound_min

	def worldToLocal( self, pWS ):
		return (pWS-self.bound_min)/self.size

	def voxelToWorld(self, pVS):
		return self.localToWorld(self.voxelToLocal(pVS))

	def worldToVoxel(self, pWS):
		return self.localToVoxel(self.worldToLocal(pWS))

	def voxelToIndex(self, pVS):
		return int(pVS)

	def rasterize( self, fun ):
		field = np.zeros((self.res))
		for i in range(0, self.res):
				pVS = i+0.5
				pWS = self.voxelToWorld(pVS)
				field[i] = fun(pWS)
		return field

	def rasterizeVS( self, fun ):
		field = np.zeros((self.res))
		for i in range(0, self.res):
				field[i] = fun(i)
		return field

	def gradient( self, field ):
		grad_field = np.zeros((self.res))
		for i in range(0, self.res):
				if i==0:
					# forward differences for x
					grad_field[i] = (field[i+1]-field[i])/self.h
				elif i==self.res-1:
					#backward differences for x
					grad_field[i] = (field[i]-field[i-1])/self.h
				else:
					# central differences for x
					grad_field[i] = (field[i+1]-field[i-1])/(self.h*2.0)
		return grad_field














if __name__ == "__main__":

	h = 0.5
	i = 0
	j = 0
	k = 0

	# ---------------
	stencil_dx = CentralDifference(h, 0)
	stencil_dy = CentralDifference(h, 1)
	stencil_dz = CentralDifference(h, 2)

	stencil_dxdx = CentralDifference(h, 0)*CentralDifference(h, 0)
	stencil_dxdy = CentralDifference(h, 0)*CentralDifference(h, 1)
	#stencil_dy = CentralDifference(h, 1)
	#stencil_dz = CentralDifference(h, 2)


	stencil_points = stencil_dxdy.getPoints(i, j, k)
	for p in stencil_points:
		print( "stencilpoint {} {} {} {}".format(p.coord[0], p.coord[1], p.coord[2], p.weight) )
