import numpy as np
import util
#import solver
import time
import shtools




# these functions are the coefficients for the recursive relation of the sh basis function
# (see p. 4 in the starmap paper)
def a_lm( l, m ):
    return np.sqrt((l-m+1)*(l+m+1)/((2*l+3)*(2*l+1)))
def b_lm( l, m ):
    return np.sqrt((l-m)*(l+m)/((2*l+1)*(2*l-1)))
def c_lm( l, m ):
    return np.sqrt((l+m+1)*(l+m+2)/((2*l+3)*(2*l+1)))
def d_lm( l, m ):
    return np.sqrt((l-m)*(l-m-1)/((2*l+1)*(2*l-1)))
def e_lm( l, m ):
    return np.sqrt((l-m+1)*(l-m+2)/((2*l+3)*(2*l+1)))
def f_lm( l, m ):
    return np.sqrt((l+m)*(l+m-1)/((2*l+1)*(2*l-1)))

def expm1div( x ):
    y = 1+x*0.5+np.power( x, 2.0)/6.0
    indices = np.abs(x) > 2.0e-4
    x_selected = x[np.where(indices)]
    y[indices] = np.divide( np.exp(x_selected)-1.0, x_selected)
    return y

def rasterize( fun, X, Y ):
    if not X.shape == Y.shape:
        print("error, shape needs to match")
    result = np.zeros( X.shape )
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            result[i,j] = fun(X[i,j], Y[i,j])
    return result

def rasterize3d( fun, X, Y, Z ):
	if not X.shape == Y.shape or not X.shape == Z.shape or not Y.shape == Z.shape:
		print("error, shape needs to match")
	result = np.zeros( X.shape )
	for i in range(X.shape[0]):
		for j in range(X.shape[1]):
			for k in range(X.shape[2]):
				result[i,j,k] = fun(X[i,j,k], Y[i,j,k], Z[i,j,k])
	return result


def rasterize_moment( fun, X, Y, m ):
    if not X.shape == Y.shape:
        print("error, shape needs to match")
    result = np.zeros( X.shape )
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            result[i,j] = fun(m, X[i,j], Y[i,j])
    return result

def interp( grid, axis ):
    # same as np.diff but interpolates between neighbour gridpoints along given axes
    # (instead of computing the difference)
    # there is a more elegent way for this using np.rollaxis (https://goo.gl/XQxRdD)
    if axis == 0:
        result = np.zeros((grid.shape[0]-1, grid.shape[1]))
        for i in range(grid.shape[0]-1):
            result[i] = (grid[i] + grid[i+1])*0.5
    elif axis == 1:
        result = np.zeros((grid.shape[0], grid.shape[1]-1))
        for j in range(grid.shape[1]-1):
            result[:, j] = (grid[:,j] + grid[:,j])*0.5
    return result
    

class StaggeredGridAssignment2D(object):
	'''This is a helper class for assigning moment coefficients to staggered grid locations (in 2D)'''
	def __init__(self):
		# c will hold the list of components for each staggered grid location
		self.c = [[[] for j in range(2)] for i in range(2)]
		# a list of (i,j) pairs which identify the indiviual grids (used for iterating over all grids)
		self.grids = []
		# the grid offset for each grid (in voxelspace)
		self.offsets = [[[] for j in range(2)] for i in range(2)]
		for i in range(2):
			for j in range(2):
				self.grids.append((i,j))
				self.offsets[i][j] = ((1-i)*0.5, (1-j)*0.5)

		# the offsets in 2D (hardcoded)
		self.offsets[0][0] = (0.5,0.5)
		self.offsets[0][1] = (0.5,0.0)
		self.offsets[1][0] = (0.0,0.5)
		self.offsets[1][1] = (0.0,0.0)
        
    
	def assign( self, component, i, j ):
		if not component in self.c[i%2][j%2]:
			self.c[i%2][j%2].append(component)
	def get_grid_index( self, component ):
		'''returns the grid (identified by (i,j) pair) which is associated with the given coefficient'''
		for i in range(2):
			for j in range(2):
				if component in self.c[i][j]:
					return (i, j)
		return None
	def get_grid_components(self, i, j):
		return self.c[i][j]
	def get_grids(self):
		return self.grids
	def get_num_grids(self):
		return len(self.grids)
	def get_offset(self, i, j):
		return self.offsets[i][j]
	def get_grid_resolution( self, i, j, res ):
		res_x = res
		res_y = res
		if j == 1:
			res_y = res_y+1
		if i == 1:
			res_x = res_x+1
		return (res_x, res_y)
	def is_even( self, c ):
		if c in self.c[0][0] or c in self.c[1][1]:
			return True
		return False
	def is_odd( self, c ):
		if self.is_even(c):
			return False
		return True

class Starmap2D(object):
	''' This class implements the starmap matlab code in python.

	The class is initialized with the approximation order and the discretization domain (including bounding box and gridsize).
	The run method takes all problem specific input sich as source and absorption/scattering coefficient functions and runs
	the problem for a given number of time steps.
	'''

	def __init__(self, order = 3, domain = None):
		'''The init function takes approximation order and domain (including discretization info)'''

		self.N = order
		self.domain = domain

		self.index_to_lm = [] # this will map equation index to l,m indices
		self.lm_to_index = {} # this dict will map l,m indices to the equation index
		# NB: we dont use l(l+1)+m because in 2d, we skipp odd (l+m) entries, which changes the sequence.
		# iterate all sh coefficients for given truncation order
		for l in range(0, self.N+1):
			for m in range(-l, l+1):
				# in 2d, we only need to solve for moments where l+m is even
				if (l+m) % 2 == 0:
					self.index_to_lm.append( (l, m) )
					self.lm_to_index[(l,m)] = len(self.index_to_lm)-1
		self.numCoeffs = len(self.index_to_lm)

		self.build_S()
		self.build_M()
		self.build_staggered_grid()

	def shIndex(self, l, m):
		#NB: we dont use l(l+1)+m because in 2d, we skipp odd (l+m) entries, which changes the sequence.
		return self.lm_to_index[(l,m)]


	def build_S(self):
		'''builds the S matrix, which converts from complex-valued to real valued coefficients'''
		# build S matrix ( we iterate over l, m to make sure that the order is correct)

		self.S = np.zeros((self.numCoeffs, self.numCoeffs),dtype=complex)
		count = 0
		for l in range(0, self.N+1):
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

	def build_M(self):
		'''builds the complex-valued M matrices which model the dependencies between moments.
		The M matrices encode the dependencies between coefficients. This is the same everywhere in the domain.
		This dependency has a specific structure, which is exploited in the staggered grid approach used by starmap.
		Each component of the solution vector is assigned to one of a the four existing staggered grids.
		This assignment is driven by the dependencies between components and by some initial assignment of the first component.
		'''
		self.Mx_complex = np.zeros((self.numCoeffs, self.numCoeffs),dtype=complex)
		self.My_complex = np.zeros((self.numCoeffs, self.numCoeffs),dtype=complex)

		# iterate all moment equations and build rows for (complex) Mx and My matrices
		for i in range(self.numCoeffs):
			l = self.index_to_lm[i][0]
			m = self.index_to_lm[i][1]

			# build complex-valued Mx and My matrices
			if (l-1,m-1) in self.lm_to_index:
				self.Mx_complex[i, self.lm_to_index[(l-1,m-1)]] = -0.5*c_lm(l-1, m-1)
				self.My_complex[i, self.lm_to_index[(l-1,m-1)]] = 0.5*1j*c_lm(l-1, m-1)
			if (l+1,m-1) in self.lm_to_index:
				if i == 0:
					tt = 0.5*d_lm(l+1, m-1)
					print("test={}".format(tt))
				self.Mx_complex[i, self.lm_to_index[(l+1,m-1)]] = 0.5*d_lm(l+1, m-1)
				self.My_complex[i, self.lm_to_index[(l+1,m-1)]] = -0.5*1j*d_lm(l+1, m-1)   
			if (l-1,m+1) in self.lm_to_index:
				self.Mx_complex[i, self.lm_to_index[(l-1,m+1)]] = 0.5*e_lm(l-1, m+1)
				self.My_complex[i, self.lm_to_index[(l-1,m+1)]] = 0.5*1j*e_lm(l-1, m+1)
			if (l+1,m+1) in self.lm_to_index:
				self.Mx_complex[i, self.lm_to_index[(l+1,m+1)]] = -0.5*f_lm(l+1, m+1)
				self.My_complex[i, self.lm_to_index[(l+1,m+1)]] = -0.5*1j*f_lm(l+1, m+1)

		# compute real valued M matrices (see middle of p.6 in starmap paper)
		self.Mx_real = np.real(self.S.dot(self.Mx_complex.dot(self.S_inv)))
		self.My_real = np.real(self.S.dot(self.My_complex.dot(self.S_inv)))

		# the Ix/Iy arrays hold the valid indices for each row in Mx/My and are used later for the update step
		# in other words: Ix/Iy hold for each moment, the indices of moments on which the moment depends
		self.Ix = [[] for i in range(self.numCoeffs)]
		self.Iy = [[] for i in range(self.numCoeffs)]

		for row in range(self.numCoeffs):
			for col in range(self.numCoeffs):
				if np.abs(self.Mx_real[row, col]) > 0:
					# keep track of valud indices in Mx for each row
					if not col in self.Ix[row]:
						self.Ix[row].append(col)
				if np.abs(self.My_real[row, col]) > 0:
					# keep track of valud indices in My for each row
					if not col in self.Iy[row]:
						self.Iy[row].append(col)



	def build_staggered_grid(self):
		'''assigns components to staggered grid locations and builds some helper structures'''
		self.sga = StaggeredGridAssignment2D()

		# we place the first component by hand
		self.sga.assign(0, 0, 0)

		# we start propagating all other components
		todo = [0]
		done = []
		count = 0
		maxCount = 20
		while todo:# and count<maxCount:
			count = count +1
			row = todo[0]
			todo.remove(row)
			done.append(row)

			# get the grid index which is associated with the current component/row
			(i, j) = self.sga.get_grid_index(row)

			# work out which components the component associated with the current row depends on according to Mx
			# we perform the dot product between valid elements in the row of Mx (stored in Ix) and the solution vector
			for col in self.Ix[row]:
				# the component depends on the x-derivative of the component col in the solution vector.
				# Due to our central differences discretization, we will register col accordingly
				self.sga.assign(col, i+1, j)
				self.sga.assign(col, i-1, j)

				if row == 0:
					print("coefficient {} depends on coefficient {} {}".format(row, " dx ", col) )

				# do this next
				if not col in done and not col in todo:
					todo.append(col)

			# work out which components the component associated with the current row depends on according to My
			# we perform the dot product between valid elements in the row of My (stored in Iy) and the solution vector
			for col in self.Iy[row]:
				# the component depends on the y-derivative of the component col in the solution vector.
				# Due to our central differences discretization, we will register col accordingly
				self.sga.assign(col, i, j+1)
				self.sga.assign(col, i, j-1)

				if row == 0:
					print("coefficient {} depends on coefficient {} {}".format(row, " dy ", col) )

				# do this next
				if not col in done and not col in todo:
					todo.append(col)

		if self.domain != None:
			# now compute X, Y: two arrays which hold the x and y locations of all staggered grid points in the domain
			# this is used later for rasterization of sigma_a and sigma_s, etc. into the domain
			# X and Y hold the grid positions in worldspace ---
			self.X = [[0 for j in range(2)] for i in range(2)]
			self.Y = [[0 for j in range(2)] for i in range(2)]

			for (grid_i, grid_j) in self.sga.get_grids():
				(res_x, res_y) = self.sga.get_grid_resolution(grid_i, grid_j, self.domain.res)

				# compute the spatial coordinates for each staggered grid
				offset = self.sga.get_offset(grid_i, grid_j)
				self.X[grid_i][grid_j] = np.array([np.arange(res_x) for y in range(res_y)]).T*self.domain.h + offset[0]*self.domain.h
				self.Y[grid_i][grid_j] = np.array([np.arange(res_y) for x in range(res_x)])*self.domain.h + offset[1]*self.domain.h

	def run( self, dt, numTimeSteps = 1, sigma_a = None, sigma_s = None, source = None ):
		'''The run method takes all problem specific inputs, such as absorption- and scattering coefficient fields and source field.'''

		# solution U ---
		# u is an array of grids. One for each moment coefficient
		u = [0 for i in range(self.numCoeffs)]

		# source Q ---
		# Q is an array of grids. One for each moment coefficient
		Q = [0 for i in range(self.numCoeffs)]

		# sigma_a ---
		# the absorption coefficient is isotropic and therefore we only need the zero moment
		# however, for the update step, we need it at all staggered grid locations
		# therefore we have the zero moment of sigma_a at each staggered grid
		# rasterization is done in the loop below
		sa = [[0 for j in range(2)] for i in range(2)]

		# sigma_s ---
		# the scattering coefficient is isotropic and therefore we only need the zero moment
		# however, for the update step, we need it at all staggered grid locations
		# therefore we have the zero moment of sigma_a at each staggered grid
		# rasterization is done in the loop below
		ss = [[0 for j in range(2)] for i in range(2)]

		# sigma_t = sigma_a + sigma_s ---
		# the extinction coefficient coeffcients are simply the sum of coefficients from scattering and absorption
		st = [[0 for j in range(2)] for i in range(2)]
		self.solve = {'st':st}

		# et is the expm1div function applied to all elements of st
		et = [[0 for j in range(2)] for i in range(2)]


		# phase funtion ---
		# TODO
		# currently we assume isotropic phase function and therfore its first moment is 1
		# and all higher moments are 0


		# now setup grid for each component. take staggered discretization into account
		for (grid_i, grid_j) in self.sga.get_grids():
			# rasterize zero moment of sigma_a at the different grid locations
			sa[grid_i][grid_j] = rasterize(sigma_a, self.X[grid_i][grid_j], self.Y[grid_i][grid_j])

			# rasterize zero moment of sigma_s at the different grid locations
			ss[grid_i][grid_j] = rasterize(sigma_s, self.X[grid_i][grid_j], self.Y[grid_i][grid_j])

			# compute st by summing absorption and scattering moments
			st[grid_i][grid_j] = sa[grid_i][grid_j] + ss[grid_i][grid_j]

			# compute the decay terms
			et[grid_i][grid_j] = expm1div(-st[grid_i][grid_j]*dt*0.5)

			# now iterate over all higher moment coefficients which are associated with the current
			# staggered grid
			for g in self.sga.get_grid_components(grid_i, grid_j):
				(res_x, res_y) = self.sga.get_grid_resolution(grid_i, grid_j, self.domain.res)
				u[g] = np.zeros( (res_x, res_y) )

				if g == 0:
					Q[g] = rasterize(source, self.X[grid_i][grid_j], self.Y[grid_i][grid_j])
				else:
					# we assume isotropic sources for now...all higher moments are zero
					Q[g] = np.zeros((res_x, res_y))

		# ea is the matrix sa (at grid00) at timestep 0 (thats why -dt*0.5)
		# with exp1mdiv applied to all elements
		# we use sa[0][0], because we know that G00 is the grid associated with c=0
		# (and we only use ea with c=0)
		ea = expm1div(-sa[0][0]*dt*0.5)


		c00 = self.sga.get_grid_components(0,0)
		c11 = self.sga.get_grid_components(1,1)
		c01 = self.sga.get_grid_components(0,1)
		c10 = self.sga.get_grid_components(1,0)
		components_even = c00+c11
		components_odd = c01+c10

		# derivatives for each component of u in x and y
		dxU = [[] for i in range(self.numCoeffs)]
		dyU = [[] for i in range(self.numCoeffs)]

		# for each time step
		t = 0.0
		for tt in range(numTimeSteps):
			#print('timestep {} t_old={}  dt={}  t_new={}'.format(tt, t, dt, t+dt))
			t = t+dt

			for step in [1, 2, 1]:
			#for step in [1]:
				if step == 1:
					# update odd grids using single half-step ---

					# we compute derivatives of even grids at odd gridpoint positions
					for c in c00:
						# u_bc will be the grid for coefficient c _including_ boundaries

						# x-derivative ------------------------------------
						# setup u_bc and add boundary in x (i index) direction (+2 rows)
						u_bc = np.zeros( (u[c].shape[0]+2, u[c].shape[1]) )

						# now copy the content from u to the non-boundary cells of u_bc
						u_bc[1:u[c].shape[0]+1, :] = u[c]

						# duplicate the last inner row onto the boundary row at index i=0
						u_bc[0, :] = u[c][0, :]

						# duplicate the last inner row onto the boundary row at index i=res_x
						u_bc[u[c].shape[0]+1, :] = u[c][u[c].shape[0]-1, :]

						# now compute the derivative in x (using central differencing)
						dxU[c] = np.diff(u_bc, 1, 0)/self.domain.h

						# y-derivative ------------------------------------
						# setup u_bc and add boundary in y (j index) direction (+2 columns)
						u_bc = np.zeros( (u[c].shape[0], u[c].shape[1]+2) )

						# now copy the content from u to the non-boundary cells of u_bc
						u_bc[:, 1:u[c].shape[1]+1] = u[c]

						# duplicate the last inner column onto the boundary row at index j=0
						u_bc[:, 0] = u[c][:, 0]

						# duplicate the last inner column onto the boundary row at index j=res_y
						u_bc[:, u[c].shape[1]+1] = u[c][:, u[c].shape[1]-1]

						# now compute the derivative in y (using central differencing)
						dyU[c] = np.diff(u_bc, 1, 1)/self.domain.h

					for c in c11:
						u_bc = u[c]
						dxU[c] = np.diff(u_bc, 1, 0)/self.domain.h
						dyU[c] = np.diff(u_bc, 1, 1)/self.domain.h


		            # iterate all odd components and do update step...
					for c in components_odd:
						(grid_i, grid_j) = self.sga.get_grid_index(c)

						W = np.zeros( u[c].shape )
						# we now iterate over all valid(>0) elements in Mx
						# these are the weights with which we will add the
						# dxU of the moment which is associated with the current index
						for c2 in self.Ix[c]:
							W -= self.Mx_real[ c, c2 ]*dxU[c2]
						for c2 in self.Iy[c]:
							W -= self.My_real[ c, c2 ]*dyU[c2]

						# half step
						# it is unclear to me, why we use the same st/et for higher moment coefficients
						# when there is no anisotropic phase function...
						u[c] = u[c] + dt*0.5*np.multiply( W + Q[c] - np.multiply( st[grid_i][grid_j], u[c] ), et[grid_i][grid_j])

				elif step == 2:
					# update even grids using two half-steps ---

					# compute derivatives...
					for c in c10:
						# derivative in x
						u_bc = u[c]
						dxU[c] = np.diff(u_bc, 1, 0)/self.domain.h
						#derivative in y
						u_bc = np.zeros( (u[c].shape[0], u[c].shape[1]+2) )
						u_bc[:, 1:u[c].shape[1]+1] = u[c]
						u_bc[:, 0] = u[c][:, 0]
						u_bc[:, u[c].shape[1]+1] = u[c][:, u[c].shape[1]-1]
						dyU[c] = np.diff(u_bc, 1, 1)/self.domain.h

					for c in c01:
						# derivative in x
						u_bc = np.zeros( (u[c].shape[0]+2, u[c].shape[1]) )
						u_bc[1:u[c].shape[0]+1, :] = u[c]
						u_bc[0, :] = u[c][0, :]
						u_bc[u[c].shape[0]+1, :] = u[c][u[c].shape[0]-1, :]
						dxU[c] = np.diff(u_bc, 1, 0)/self.domain.h
						#derivative in y
						u_bc = u[c]
						dyU[c] = np.diff(u_bc, 1, 1)/self.domain.h


					# now iterate all even components and do update step...
					for c in components_even:
						W = np.zeros( u[c].shape )
						# we now iterate over all valid(>0) elements in Mx
						# these are the weights with which we will add the
						# dxU of the moment which is associated with the current index
						for c2 in self.Ix[c]:
							W -= self.Mx_real[ c, c2 ]*dxU[c2]
						for c2 in self.Iy[c]:
							W -= self.My_real[ c, c2 ]*dyU[c2]

						# perform two half-steps
						for k in range(2):
							# the zero moment update only depends on absorption
							# this is because the zero moment of the phase function is 1 and
							# the scattering coefficient from scattering term cancels out with
							# the scattering coefficient from sigma_t (see note in the python notebook)
							if c == 0:
								# we use sa[0][0] because we know, that G00 is the grid associated with c=0
								u[c] = u[c] + dt*0.5*np.multiply( W + Q[c] - np.multiply( sa[0][0], u[c] ), ea)
							else:
								(grid_i, grid_j) = self.sga.get_grid_index(c)
								# it is unclear to me, why we use the same st/et for higher moment coefficients
								# when there is no anisotropic phase function...
								u[c] = u[c] + dt*0.5*np.multiply( W + Q[c] - np.multiply( st[grid_i][grid_j], u[c] ), et[grid_i][grid_j])
		return u

	def compute_sh_coefficients_at_cell_centers(self, u):
		'''This method returns the radiance field coefficients L_lm at the cell centers of the discretization grid.
		The method interpolates the solution u to cell centers and uses the inverse of S to produce the final complex-valued
		coefficients.'''

		# first we want to evaluate all the moment coefficients at the cell centers
		# so we unstagger the grids
		c00 = self.sga.get_grid_components(0,0)
		c11 = self.sga.get_grid_components(1,1)
		c01 = self.sga.get_grid_components(0,1)
		c10 = self.sga.get_grid_components(1,0)

		u_center = np.zeros((self.numCoeffs, self.domain.res, self.domain.res))

		for c in c00:
			# this is already at the cell centers
			u_center[c] = u[c]

		for c in c11:
			u_center[c] = interp(interp(u[c], 0), 1)

		for c in c01:
			u_center[c] = interp(u[c], 1)

		for c in c10:
			u_center[c] = interp(u[c], 0)

		# convert from real to complex valued sh coefficients over -m, m
		L_lm = np.zeros((self.numCoeffs, self.domain.res, self.domain.res), dtype=complex)
		for i in range(self.domain.res):
			for j in range(self.domain.res):
				L_lm[:, i, j] = self.S_inv.dot(u_center[:, i, j])

		return L_lm




class StaggeredGridAssignment3D(object):
	'''This is a helper class for assigning moment coefficients to staggered grid locations (in 2D)'''
	def __init__(self):
		# c will hold the list of components for each staggered grid location
		self.c = [[[[] for k in range(2)] for j in range(2)] for i in range(2)]
		# a list of (i,j,k) pairs which identify the indiviual grids (used for iterating over all grids)
		self.grids = []
		# the grid offset for each grid (in voxelspace)
		self.offsets = [[[[] for k in range(2)] for j in range(2)] for i in range(2)]
		for i in range(2):
			for j in range(2):
				for k in range(2):
					self.grids.append((i,j,k))
					#self.offsets[i][j][k] = ((1-i)*0.5, (1-j)*0.5, (1-k)*0.5)

		# the offsets in 3D (hardcoded)
		self.offsets[0][0][0] = (0.5,0.5,0.5) # even (voxel center)
		self.offsets[0][0][1] = (0.5,0.5,0.0) # odd (voxel face)
		self.offsets[0][1][0] = (0.5,0.0,0.5) # odd (voxel face)
		self.offsets[0][1][1] = (0.5,0.0,0.0) # even (voxel edge)
		self.offsets[1][0][0] = (0.0,0.5,0.5) # odd (voxel face)
		self.offsets[1][0][1] = (0.0,0.5,0.0) # even (voxel edge)
		self.offsets[1][1][0] = (0.0,0.0,0.5) # even (voxel edge)
		self.offsets[1][1][1] = (0.0,0.0,0.0) # odd (voxel corner)

        
    
	def assign( self, component, i, j, k ):
		if not component in self.c[i%2][j%2][k%2]:
			self.c[i%2][j%2][k%2].append(component)
	def get_grid_index( self, component ):
		'''returns the grid (identified by (i,j) pair) which is associated with the given coefficient'''
		for i in range(2):
			for j in range(2):
				for k in range(2):
					if component in self.c[i][j][k]:
						return (i, j, k)
		return None
	def get_grid_components(self, i, j, k):
		return self.c[i][j][k]
	def get_grids(self):
		return self.grids
	def get_num_grids(self):
		return len(self.grids)
	def get_offset(self, i, j, k):
		return self.offsets[i][j][k]
	def get_grid_resolution( self, i, j, k, res ):
		res_x = res
		res_y = res
		res_z = res
		if k == 1:
			res_z = res_z+1
		if j == 1:
			res_y = res_y+1
		if i == 1:
			res_x = res_x+1
		return (res_x, res_y, res_z)
	def is_assigned(self, c ):
		for i in range(2):
			for j in range(2):
				for k in range(2):
					if c in self.c[i][j][k]:
						return True
		return False
	def is_even( self, c ):
		if not self.is_assigned(c):
			raise ValueError('c is not assigned to any grid')
		if c in self.c[0][0][0] or c in self.c[0][1][1] or c in self.c[1][0][1] or c in self.c[1][1][0]:
			return True
		return False
	def is_odd( self, c ):
		if self.is_even(c):
			return False
		return True

class Starmap3D(object):
	''' This class implements the starmap matlab code in python.

	The class is initialized with the approximation order and the discretization domain (including bounding box and gridsize).
	The run method takes all problem specific input sich as source and absorption/scattering coefficient functions and runs
	the problem for a given number of time steps.
	'''

	def __init__(self, order = 3, domain = None):
		'''The init function takes approximation order and domain (including discretization info)'''

		self.N = order
		self.domain = domain

		self.index_to_lm = [] # this will map equation index to l,m indices
		self.lm_to_index = {} # this dict will map l,m indices to the equation index
		# iterate all sh coefficients for given truncation order
		for l in range(0, self.N+1):
			for m in range(-l, l+1):
				self.index_to_lm.append( (l, m) )
				self.lm_to_index[(l,m)] = len(self.index_to_lm)-1
		self.numCoeffs = len(self.index_to_lm)

		self.build_S()
		self.build_M()
		self.build_staggered_grid()



	def build_S(self):
		'''builds the S matrix, which converts from complex-valued to real valued coefficients'''
		print("Starmap3D::build_S")
		# build S matrix ( we iterate over l, m to make sure that the order is correct)
		self.S = np.zeros((self.numCoeffs, self.numCoeffs),dtype=complex)
		count = 0
		for l in range(0, self.N+1):
			# note that the starmap code builds the solution vector from reverse order
			for m in range(l, -1, -1):
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


	def build_M(self):
		'''builds the complex-valued M matrices which model the dependencies between moments.
		The M matrices encode the dependencies between coefficients. This is the same everywhere in the domain.
		This dependency has a specific structure, which is exploited in the staggered grid approach used by starmap.
		Each component of the solution vector is assigned to one of a the four existing staggered grids.
		This assignment is driven by the dependencies between components and by some initial assignment of the first component.
		'''
		print("Starmap3D::build_M")
		self.Mx_complex = np.zeros((self.numCoeffs, self.numCoeffs),dtype=complex)
		self.My_complex = np.zeros((self.numCoeffs, self.numCoeffs),dtype=complex)
		self.Mz_complex = np.zeros((self.numCoeffs, self.numCoeffs),dtype=complex)

		# iterate all moment equations and build rows for (complex) Mx and My matrices
		for i in range(self.numCoeffs):
			l = self.index_to_lm[i][0]
			m = self.index_to_lm[i][1]

			# build complex-valued Mx and My and Mz matrices
			if (l-1,m-1) in self.lm_to_index:
				self.Mx_complex[i, self.lm_to_index[(l-1,m-1)]] = -0.5*c_lm(l-1, m-1)
				self.My_complex[i, self.lm_to_index[(l-1,m-1)]] = 0.5*1j*c_lm(l-1, m-1)
			if (l+1,m-1) in self.lm_to_index:
				self.Mx_complex[i, self.lm_to_index[(l+1,m-1)]] = 0.5*d_lm(l+1, m-1)   
				self.My_complex[i, self.lm_to_index[(l+1,m-1)]] = -0.5*1j*d_lm(l+1, m-1)   
			if (l-1,m+1) in self.lm_to_index:
				self.Mx_complex[i, self.lm_to_index[(l-1,m+1)]] = 0.5*e_lm(l-1, m+1)
				self.My_complex[i, self.lm_to_index[(l-1,m+1)]] = 0.5*1j*e_lm(l-1, m+1)
			if (l+1,m+1) in self.lm_to_index:
				self.Mx_complex[i, self.lm_to_index[(l+1,m+1)]] = -0.5*f_lm(l+1, m+1)
				self.My_complex[i, self.lm_to_index[(l+1,m+1)]] = -0.5*1j*f_lm(l+1, m+1)
			if (l-1,m) in self.lm_to_index:
				self.Mz_complex[i, self.lm_to_index[(l-1,m)]] = a_lm(l-1, m)
			if (l+1,m) in self.lm_to_index:
				self.Mz_complex[i, self.lm_to_index[(l+1,m)]] = b_lm(l+1, m)

		# compute real valued M matrices (see middle of p.6 in starmap paper)
		self.Mx_real = np.real(self.S.dot(self.Mx_complex.dot(self.S_inv)))
		self.My_real = np.real(self.S.dot(self.My_complex.dot(self.S_inv)))
		self.Mz_real = np.real(self.S.dot(self.Mz_complex.dot(self.S_inv)))

		# the Ix/Iy arrays hold the valid indices for each row in Mx/My and are used later for the update step
		# in other words: Ix/Iy hold for each moment, the indices of moments on which the moment depends
		self.Ix = [[] for i in range(self.numCoeffs)]
		self.Iy = [[] for i in range(self.numCoeffs)]
		self.Iz = [[] for i in range(self.numCoeffs)]

		for row in range(self.numCoeffs):
			for col in range(self.numCoeffs):
				if np.abs(self.Mx_real[row, col]) > 0:
					# keep track of valud indices in Mx for each row
					if not col in self.Ix[row]:
						self.Ix[row].append(col)
				if np.abs(self.My_real[row, col]) > 0:
					# keep track of valud indices in My for each row
					if not col in self.Iy[row]:
						self.Iy[row].append(col)
				if np.abs(self.Mz_real[row, col]) > 0:
					# keep track of valud indices in Mz for each row
					if not col in self.Iz[row]:
						self.Iz[row].append(col)

	def build_staggered_grid(self):
		'''assigns components to staggered grid locations and builds some helper structures'''
		print("Starmap3D::build_staggered_grid")

		self.sga = StaggeredGridAssignment3D()

		# we place the first component by hand
		self.sga.assign(0, 0, 0, 0)


		# we start propagating all other components
		todo = [0]
		done = []
		count = 0
		maxCount = 20
		while todo:# and count<maxCount:
			count = count +1
			row = todo[0]
			todo.remove(row)
			done.append(row)

			# get the grid index which is associated with the current component/row
			(i, j, k) = self.sga.get_grid_index(row)

			#print("assigning row={} i={} j={} k={}".format(row, i, j, k))
			#if self.sga.is_even(row):
			#	print("\tis even")
			#else:
			#	print("\tis odd")

			# work out which components the component associated with the current row depends on according to Mx
			# we perform the dot product between valid elements in the row of Mx (stored in Ix) and the solution vector
			for col in self.Ix[row]:
				# the component depends on the x-derivative of the component col in the solution vector.
				# Due to our central differences discretization, we will register col accordingly
				self.sga.assign(col, i+1, j, k)
				self.sga.assign(col, i-1, j, k)
				# do this next
				if not col in done and not col in todo:
					todo.append(col)

			# work out which components the component associated with the current row depends on according to My
			# we perform the dot product between valid elements in the row of My (stored in Iy) and the solution vector
			for col in self.Iy[row]:
				# the component depends on the y-derivative of the component col in the solution vector.
				# Due to our central differences discretization, we will register col accordingly
				self.sga.assign(col, i, j+1, k)
				self.sga.assign(col, i, j-1, k)
				# do this next
				if not col in done and not col in todo:
					todo.append(col)

			# work out which components the component associated with the current row depends on according to Mz
			# we perform the dot product between valid elements in the row of Mz (stored in Iz) and the solution vector
			for col in self.Iz[row]:
				# the component depends on the y-derivative of the component col in the solution vector.
				# Due to our central differences discretization, we will register col accordingly
				self.sga.assign(col, i, j, k+1)
				self.sga.assign(col, i, j, k-1)
				# do this next
				if not col in done and not col in todo:
					todo.append(col)

		# as a sanity check: make sure all components are assigned
		for c in range(self.numCoeffs):
			if not self.sga.is_assigned(c):
				raise ValueError("component not assigned to staggered grid c={}".format(c))

		# as a sanity check: verify that even components only depend on odd components and vice versa
		self.check_coupling(self.Mx_real)
		self.check_coupling(self.My_real)
		self.check_coupling(self.Mz_real)


		if self.domain != None:
			print("Starmap3D::build_staggered_grid creating grid position arrays")

			# now compute X, Y, Z: two arrays which hold the x and y and z locations of all staggered grid points in the domain
			# this is used later for rasterization of sigma_a and sigma_s, etc. into the domain
			# X and Y and Z hold the grid positions in worldspace ---
			self.X = [[[[] for k in range(2)] for j in range(2)] for i in range(2)]
			self.Y = [[[[] for k in range(2)] for j in range(2)] for i in range(2)]
			self.Z = [[[[] for k in range(2)] for j in range(2)] for i in range(2)]

			for (grid_i, grid_j, grid_k) in self.sga.get_grids():
				(res_x, res_y, res_z) = self.sga.get_grid_resolution(grid_i, grid_j, grid_k, self.domain.res)
				print("\tgrid={} {} {}".format(grid_i, grid_j, grid_k))

				# compute the spatial coordinates for each staggered grid
				offset = self.sga.get_offset(grid_i, grid_j, grid_k)

				#self.X[grid_i][grid_j][grid_k] = [[[ i for k in range(res_z)] for j in range(res_y)] for i in range(res_x)]
				#self.X[grid_i][grid_j][grid_k] = np.array(self.X[grid_i][grid_j][grid_k])*self.domain.h + offset[0]*self.domain.h
				self.X[grid_i][grid_j][grid_k] = solver.create_position_array(res_x, res_y, res_z, 0)*self.domain.h  + offset[0]*self.domain.h

				#self.Y[grid_i][grid_j][grid_k] = [[[ j for k in range(res_z)] for j in range(res_y)] for i in range(res_x)]
				#self.Y[grid_i][grid_j][grid_k] = np.array(self.Y[grid_i][grid_j][grid_k])*self.domain.h + offset[1]*self.domain.h
				self.Y[grid_i][grid_j][grid_k] = solver.create_position_array(res_x, res_y, res_z, 1)*self.domain.h  + offset[1]*self.domain.h

				#self.Z[grid_i][grid_j][grid_k] = [[[ k for k in range(res_z)] for j in range(res_y)] for i in range(res_x)]
				#self.Z[grid_i][grid_j][grid_k] = np.array(self.Z[grid_i][grid_j][grid_k])*self.domain.h + offset[2]*self.domain.h
				self.Z[grid_i][grid_j][grid_k] = solver.create_position_array(res_x, res_y, res_z, 2)*self.domain.h  + offset[2]*self.domain.h

	def check_coupling(self, M, threshold = 0.0):
		violations = np.zeros(M.shape)
		for row in range(self.numCoeffs):
			for col in range(self.numCoeffs):
				if np.abs(M[row, col]) > threshold:
					# row depends on col
					# check if they are both even or both odd...that would be wrong for the method
					if self.sga.is_even(row) and self.sga.is_even(col):
						violations[row, col] = 1.0
						#print("even/odd dependeny violated row={} col={}".format(row, col))
					if self.sga.is_odd(row) and self.sga.is_odd(col):
						violations[row, col] = 1.0
						#print("odd/even dependeny violated row={} col={}".format(row, col))
		return violations

	# u_bc will be the grid for coefficient c _including_ boundaries
	def extend_x(self, u):
		# setup u_bc and add boundary in x (i index) direction (+2 rows)
		u_bc = np.zeros( (u.shape[0]+2, u.shape[1], u.shape[2]) )
		# now copy the content from u to the non-boundary cells of u_bc
		u_bc[1:u.shape[0]+1, :, :] = u
		# duplicate the last inner yz-slice onto the boundary sice at index i=0
		u_bc[0, :, :] = u[0, :, :]
		# duplicate the last inner yz slice onto the boundary slice at index i=res_x
		u_bc[u.shape[0]+1, :, :] = u[u.shape[0]-1, :, :]
		return u_bc
	def extend_y(self, u):
		u_bc = np.zeros( (u.shape[0], u.shape[1]+2, u.shape[2]) )

		# now copy the content from u to the non-boundary cells of u_bc
		u_bc[:, 1:u.shape[1]+1, :] = u

		# duplicate the last inner inner y-z slice onto the boundary row at index j=0
		u_bc[:, 0,:] = u[:, 0, :]

		# duplicate the last inner inner y-z slice onto the boundary row at index j=res_y
		u_bc[:, u.shape[1]+1, :] = u[:, u.shape[1]-1, :]
		return u_bc
	def extend_z(self, u):
		# setup u_bc and add boundary in z (k index) direction (+2 columns)
		u_bc = np.zeros( (u.shape[0], u.shape[1], u.shape[2]+2) )

		# now copy the content from u to the non-boundary cells of u_bc
		u_bc[:, :, 1:u.shape[2]+1] = u

		# duplicate the last inner inner xy slice onto the boundary row at index j=0
		u_bc[:, :, 0] = u[:, :, 0]

		# duplicate the last inner inner y-z slice onto the boundary row at index j=res_y
		u_bc[:, :, u.shape[2]+1] = u[:, :, u.shape[2]-1]
		return u_bc
	def run( self, dt, numTimeSteps = 1, sigma_a = None, sigma_s = None, source = None ):
		'''The run method takes all problem specific inputs, such as absorption- and scattering coefficient fields and source field.'''
		print("Starmap3D::run")


		# solution U ---
		# u is an array of grids. One for each moment coefficient
		u = [0 for i in range(self.numCoeffs)]

		# source Q ---
		# Q is an array of grids. One for each moment coefficient
		Q = [0 for i in range(self.numCoeffs)]

		# sigma_a ---
		# the absorption coefficient is isotropic and therefore we only need the zero moment
		# however, for the update step, we need it at all staggered grid locations
		# therefore we have the zero moment of sigma_a at each staggered grid
		# rasterization is done in the loop below
		sa = [[[[] for k in range(2)] for j in range(2)] for i in range(2)]

		# sigma_s ---
		# the scattering coefficient is isotropic and therefore we only need the zero moment
		# however, for the update step, we need it at all staggered grid locations
		# therefore we have the zero moment of sigma_a at each staggered grid
		# rasterization is done in the loop below
		ss = [[[[] for k in range(2)] for j in range(2)] for i in range(2)]

		# sigma_t = sigma_a + sigma_s ---
		# the extinction coefficient coeffcients are simply the sum of coefficients from scattering and absorption
		st = [[[[] for k in range(2)] for j in range(2)] for i in range(2)]

		# et is the expm1div function applied to all elements of st
		et = [[[[] for k in range(2)] for j in range(2)] for i in range(2)]


		# phase funtion ---
		# TODO
		# currently we assume isotropic phase function and therfore its first moment is 1
		# and all higher moments are 0


		# now setup grid for each component. take staggered discretization into account
		print("Starmap3D::run setting up grids")
		for (grid_i, grid_j, grid_k) in self.sga.get_grids():
			print("\tgrid={} {} {}".format(grid_i, grid_j, grid_k))
			(res_x, res_y, res_z) = self.sga.get_grid_resolution(grid_i, grid_j, grid_k, self.domain.res)
			#print("grid={} {} {} res={} {} {}".format(grid_i, grid_j, grid_k, res_x, res_y, res_z ))

			# rasterize zero moment of sigma_a at the different grid locations
			print("\trasterizing sigma_a")
			#sa[grid_i][grid_j][grid_k] = rasterize3d(sigma_a, self.X[grid_i][grid_j][grid_k], self.Y[grid_i][grid_j][grid_k], self.Z[grid_i][grid_j][grid_k])
			sa[grid_i][grid_j][grid_k] = solver.rasterize( "sigma_a", self.X[grid_i][grid_j][grid_k], self.Y[grid_i][grid_j][grid_k], self.Z[grid_i][grid_j][grid_k] )
			#print(sa[grid_i][grid_j][grid_k].shape)

			# rasterize zero moment of sigma_s at the different grid locations
			print("\trasterizing sigma_s")
			#ss[grid_i][grid_j][grid_k] = rasterize3d(sigma_s, self.X[grid_i][grid_j][grid_k], self.Y[grid_i][grid_j][grid_k], self.Z[grid_i][grid_j][grid_k])
			ss[grid_i][grid_j][grid_k] = solver.rasterize( "sigma_s", self.X[grid_i][grid_j][grid_k], self.Y[grid_i][grid_j][grid_k], self.Z[grid_i][grid_j][grid_k] )
			#print(ss[grid_i][grid_j][grid_k].shape)

			# compute st by summing absorption and scattering moments
			st[grid_i][grid_j][grid_k] = sa[grid_i][grid_j][grid_k] + ss[grid_i][grid_j][grid_k]

			# compute the decay terms
			print("\trasterizing computing decay term")
			et[grid_i][grid_j][grid_k] = expm1div(-st[grid_i][grid_j][grid_k]*dt*0.5)

			# now iterate over all higher moment coefficients which are associated with the
			# current staggered grid
			for g in self.sga.get_grid_components(grid_i, grid_j, grid_k):	
				print("\tinitializing component c={}".format(g))
				u[g] = np.zeros( (res_x, res_y, res_z) )

				if g == 0:
					#Q[g] = rasterize3d(source, self.X[grid_i][grid_j][grid_k], self.Y[grid_i][grid_j][grid_k], self.Z[grid_i][grid_j][grid_k])
					Q[g] = solver.rasterize( "q", self.X[grid_i][grid_j][grid_k], self.Y[grid_i][grid_j][grid_k], self.Z[grid_i][grid_j][grid_k] )
				else:
					# we assume isotropic sources for now...all higher moments are zero
					Q[g] = np.zeros((res_x, res_y, res_z))

		# just for debugging: write rasterized RTE fields to disk ---
		solver.write_scalarfield( "io/check_sigma_a.bgeo", sa[0][0][0] )
		solver.write_scalarfield( "io/check_sigma_s.bgeo", ss[0][0][0] )
		solver.write_scalarfield( "io/check_q.bgeo", Q[0] )


		# ea is the matrix sa (at grid00) at timestep 0 (thats why -dt*0.5)
		# with exp1mdiv applied to all elements
		# we use sa[0][0], because we know that G00 is the grid associated with c=0
		# (and we only use ea with c=0)
		ea = expm1div(-sa[0][0][0]*dt*0.5)


		c000 = self.sga.get_grid_components(0,0,0)
		c100 = self.sga.get_grid_components(1,0,0)
		c010 = self.sga.get_grid_components(0,1,0)
		c110 = self.sga.get_grid_components(1,1,0)

		c001 = self.sga.get_grid_components(0,0,1)
		c101 = self.sga.get_grid_components(1,0,1)
		c011 = self.sga.get_grid_components(0,1,1)
		c111 = self.sga.get_grid_components(1,1,1)

		components_even = c000 + c101 + c011 + c110
		components_odd = c111 + c100 + c001 + c010


		# derivatives for each component of u in x and y and z
		dxU = [[] for i in range(self.numCoeffs)]
		dyU = [[] for i in range(self.numCoeffs)]
		dzU = [[] for i in range(self.numCoeffs)]

		# for each time step
		t = 0.0
		for tt in range(numTimeSteps):
			solver.write_scalarfield( "io/solution_u.{0:03d}.bgeo".format(tt), np.real(u[0]) )
			start = time.time()
			print('timestep {} t_old={}  dt={}  t_new={}'.format(tt, t, dt, t+dt))
			t = t+dt

			for step in [1, 2, 1]:
			#for step in [2]:
				#print("step")
				if step == 1:
					# update odd grids using single half-step ---

					# we compute derivatives of even grids at odd gridpoint positions
					for c in c000:
						dxU[c] = np.diff(self.extend_x(u[c]), 1, 0)/self.domain.h
						dyU[c] = np.diff(self.extend_y(u[c]), 1, 1)/self.domain.h
						dzU[c] = np.diff(self.extend_z(u[c]), 1, 2)/self.domain.h

					for c in c110:
						dxU[c] = np.diff(u[c], 1, 0)/self.domain.h
						dyU[c] = np.diff(u[c], 1, 1)/self.domain.h
						dzU[c] = np.diff(self.extend_z(u[c]), 1, 2)/self.domain.h

					for c in c101:
						dxU[c] = np.diff(u[c], 1, 0)/self.domain.h
						dyU[c] = np.diff(self.extend_y(u[c]), 1, 1)/self.domain.h
						dzU[c] = np.diff(u[c], 1, 2)/self.domain.h

					for c in c011:
						dxU[c] = np.diff(self.extend_x(u[c]), 1, 0)/self.domain.h
						dyU[c] = np.diff(u[c], 1, 1)/self.domain.h
						dzU[c] = np.diff(u[c], 1, 2)/self.domain.h

		            # iterate all odd components and do update step...
					#print(dzU)
					for c in components_odd:
						(grid_i, grid_j, grid_k) = self.sga.get_grid_index(c)

						W = np.zeros( u[c].shape )

						# we now iterate over all valid(>0) elements in Mx
						# these are the weights with which we will add the
						# dxU of the moment which is associated with the current index
						for c2 in self.Ix[c]:
							W -= self.Mx_real[ c, c2 ]*dxU[c2]
						for c2 in self.Iy[c]:
							W -= self.My_real[ c, c2 ]*dyU[c2]
						for c2 in self.Iz[c]:
							W -= self.Mz_real[ c, c2 ]*dzU[c2]

						# half step
						# it is unclear to me, why we use the same st/et for higher moment coefficients
						# when there is no anisotropic phase function...
						u[c] = u[c] + dt*0.5*np.multiply( W + Q[c] - np.multiply( st[grid_i][grid_j][grid_k], u[c] ), et[grid_i][grid_j][grid_k])

				elif step == 2:
					# update even grids using two half-steps ---
					for c in c100:
						dxU[c] = np.diff(u[c], 1, 0)/self.domain.h
						dyU[c] = np.diff(self.extend_y(u[c]), 1, 1)/self.domain.h
						dzU[c] = np.diff(self.extend_z(u[c]), 1, 2)/self.domain.h
					for c in c010:
						dxU[c] = np.diff(self.extend_x(u[c]), 1, 0)/self.domain.h
						dyU[c] = np.diff(u[c], 1, 1)/self.domain.h
						dzU[c] = np.diff(self.extend_z(u[c]), 1, 2)/self.domain.h

					for c in c001:
						dxU[c] = np.diff(self.extend_x(u[c]), 1, 0)/self.domain.h
						dyU[c] = np.diff(self.extend_y(u[c]), 1, 1)/self.domain.h
						dzU[c] = np.diff(u[c], 1, 2)/self.domain.h

					for c in c111:
						dxU[c] = np.diff(u[c], 1, 0)/self.domain.h
						dyU[c] = np.diff(u[c], 1, 1)/self.domain.h
						dzU[c] = np.diff(u[c], 1, 2)/self.domain.h


					# now iterate all even components and do update step...
					for c in components_even:
						(grid_i, grid_j, grid_k) = self.sga.get_grid_index(c)

						W = np.zeros( u[c].shape )
						# we now iterate over all valid(>0) elements in Mx
						# these are the weights with which we will add the
						# dxU of the moment which is associated with the current index
						for c2 in self.Ix[c]:
							W -= self.Mx_real[ c, c2 ]*dxU[c2]
						for c2 in self.Iy[c]:
							W -= self.My_real[ c, c2 ]*dyU[c2]
						for c2 in self.Iz[c]:
							W -= self.Mz_real[ c, c2 ]*dzU[c2]

						# perform two half-steps
						for k in range(2):
							# the zero moment update only depends on absorption
							# this is because the zero moment of the phase function is 1 and
							# the scattering coefficient from scattering term cancels out with
							# the scattering coefficient from sigma_t (see note in the python notebook)
							if c == 0:
								# we use sa[0][0] because we know, that G00 is the grid associated with c=0
								u[c] = u[c] + dt*0.5*np.multiply( W + Q[c] - np.multiply( sa[0][0][0], u[c] ), ea)
							else:
								# it is unclear to me, why we use the same st/et for higher moment coefficients
								# when there is no anisotropic phase function...
								u[c] = u[c] + dt*0.5*np.multiply( W + Q[c] - np.multiply( st[grid_i][grid_j][grid_k], u[c] ), et[grid_i][grid_j][grid_k])
			end = time.time()
			print(end - start)

		return u





if __name__ == "__main__":

	#order = 4
	#starmap = Starmap3D(order, util.Domain3D(7.0, 100))
	#dt = 0.01609501
	#numTimeSteps = 500
	#numTimeSteps = 2
	#u = starmap.run( dt, numTimeSteps )

	#print(shtools.numSHCoeffs(4))

	sm = Starmap2D(1)

	## the offsets in 2D (hardcoded)
	#self.offsets[0][0] = (0.5,0.5)
	#self.offsets[0][1] = (0.5,0.0)
	#self.offsets[1][0] = (0.0,0.5)
	#self.offsets[1][1] = (0.0,0.0)

	for i in range(2):
		for j in range(2):
			print("i={} j={}".format(i, j))
			print(sm.sga.get_grid_components(i, j))



