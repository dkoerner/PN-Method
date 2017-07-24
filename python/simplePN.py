import numpy as np
import util
import time
import shtools
import scipy.io




class SimplePN2D(object):
	''' This class implements PN using the first order form of the transport equation.

	The class is initialized with the approximation order and the discretization domain (including bounding box and gridsize).
	The run method takes all problem specific input sich as source and absorption/scattering coefficient functions and runs
	the problem.
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
		print("numCoeffs={}".format(self.numCoeffs))

		# X and Y are scalar fields which store the x and y components of the center for each voxel
		offset = (0.5,0.5)
		self.X = np.array([np.arange(self.domain.res_x) for y in range(self.domain.res_y)]).T*self.domain.h_x + offset[0]*self.domain.h_x
		self.Y = np.array([np.arange(self.domain.res_y) for x in range(self.domain.res_x)])*self.domain.h_y + offset[1]*self.domain.h_y


		self.build_S()
		self.build_M()


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

	def matrix_prod(self, A, B):
		C = np.zeros( (A.shape[1], B.shape[0]), dtype=complex )
		for i in range(C.shape[0]):
			# for each column in the destination matrix
			for j in range(C.shape[1]):
				C[i,j] = complex(0.0)
				# compute dot product
				for k in range(A.shape[1]):
					C[i,j] += A[i,k]*B[k,j]
		return C
	def build_M(self):
		'''The M matrices encode the dependencies between coefficients.
		This is the same everywhere in the domain and therfore is build once.
		Thhis can be understood as some sort of local matrix for each voxel element.
		'''
		self.Mx_complex = np.zeros((self.numCoeffs, self.numCoeffs),dtype=complex)
		self.My_complex = np.zeros((self.numCoeffs, self.numCoeffs),dtype=complex)

		# iterate all moment equations and build rows for (complex) Mx and My matrices
		for i in range(self.numCoeffs):
			l = self.index_to_lm[i][0]
			m = self.index_to_lm[i][1]

			# build complex-valued Mx and My matrices
			if (l-1,m-1) in self.lm_to_index:
				self.Mx_complex[i, self.lm_to_index[(l-1,m-1)]] = -0.5*shtools.c_lm(l-1, m-1)
				self.My_complex[i, self.lm_to_index[(l-1,m-1)]] = 0.5*1j*shtools.c_lm(l-1, m-1)
			if (l+1,m-1) in self.lm_to_index:
				self.Mx_complex[i, self.lm_to_index[(l+1,m-1)]] = 0.5*shtools.d_lm(l+1, m-1)   
				self.My_complex[i, self.lm_to_index[(l+1,m-1)]] = -0.5*1j*shtools.d_lm(l+1, m-1)   
			if (l-1,m+1) in self.lm_to_index:
				self.Mx_complex[i, self.lm_to_index[(l-1,m+1)]] = 0.5*shtools.e_lm(l-1, m+1)
				self.My_complex[i, self.lm_to_index[(l-1,m+1)]] = 0.5*1j*shtools.e_lm(l-1, m+1)
			if (l+1,m+1) in self.lm_to_index:
				self.Mx_complex[i, self.lm_to_index[(l+1,m+1)]] = -0.5*shtools.f_lm(l+1, m+1)
				self.My_complex[i, self.lm_to_index[(l+1,m+1)]] = -0.5*1j*shtools.f_lm(l+1, m+1)

		# compute real valued M matrices (see middle of p.6 in starmap paper)
		self.Mx_real = np.real(self.S.dot(self.Mx_complex.dot(self.S_inv)))
		self.My_real = np.real(self.S.dot(self.My_complex.dot(self.S_inv)))


		# ----------------------------------------

		# compute matrix product by hand
		D = np.zeros(self.Mx_complex.shape, dtype=complex)
		E = np.zeros(self.Mx_complex.shape)
		A = self.S
		B = self.Mx_complex
		C = self.S_inv
		# for each row in the destination matrix
		for i in range(D.shape[0]):
			# for each column in the destination matrix
			for j in range(D.shape[1]):
				D[i,j] = complex(0.0)
				# compute dot product
				for k in range(B.shape[0]):
					D[i,j] += B[i,k]*C[k,j]

		# for each row in the destination matrix
		for i in range(E.shape[0]):
			# for each column in the destination matrix
			for j in range(E.shape[1]):
				sum = complex(0.0, 0.0)
				for k in range(E.shape[0]):
					for l in range(E.shape[0]):
						sum += A[i,k]*B[k,l]*C[l,j]
				E[i,j] = np.real(sum)

		D = self.matrix_prod(B, C)
		#E = self.matrix_prod(A, D)

		data = {}
		data["S_inv"] = self.S_inv
		data["S"] = self.S
		data["Mx_complex"] = self.Mx_complex
		data["Mx_real"] = self.Mx_real
		data["Mx_complex_dot_S_inv"] = self.Mx_complex.dot(self.S_inv)
		data["D"] = D
		data["E"] = np.real(E)
		#data["b"] = b.reshape((numVoxels*self.numCoeffs, 1))
		#data["b"] = b
		#scipy.io.savemat("C:/projects/epfl/epfl17/python/simplepn/data.mat", data)
		scipy.io.savemat("C:/projects/epfl/epfl17/python/sopn/data_test.mat", data)

		# ----------------------------------------

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



	#def assemble_global_matrix(self):
	#	'''This function assembles the global matrix which expresses the PN equations for all voxels
	#	and all SH coefficients in one _big_ matrix.'''

	def get_global_index( self, voxel_i, voxel_j, coeff ):
		'''Returns the equation index for the given h coeffient at the given voxel.'''
		voxel = voxel_j*self.domain.res_x + voxel_i
		return voxel*self.numCoeffs + coeff


	def run( self, sigma_a = None, sigma_s = None, source = None ):
		'''The run method takes all problem specific inputs, such as absorption-
		and scattering coefficient fields and source field and assembles a global
		matrix which is used to solve for all sh coefficients of all voxels.'''


		numVoxels = self.domain.res_x*self.domain.res_y

		# coefficient matrix A and rhs b of our global problem
		A = np.zeros( (numVoxels*self.numCoeffs, numVoxels*self.numCoeffs) )
		b = np.zeros( (numVoxels*self.numCoeffs) )

		print( "size of global coefficient matrix: {}x{}".format(A.shape[0], A.shape[1]) )

		# now we build the global coefficient matrix A and rhs b by iterating over all elements (voxels)
		# and within each voxel, we iterate over all sh coefficients
		# each row within the global system is associated with a specific sh coefficient of
		# a specific voxel
		for voxel_i in range(self.domain.res_x):
			for voxel_j in range(self.domain.res_y):

				# get location of voxel center
				pWS = np.array([self.X[voxel_i, voxel_j], self.Y[voxel_i, voxel_j]])
				#print("voxel={} {} center={} {}".format(voxel_i, voxel_j, pWS[0], pWS[1]))


				# here we iterate the local matrix vector product between M and u
				# M represents the coupling between coefficients at the same location
				# u contains the sh coefficients for the current voxel
				# the spatial derivative will cause dependencies on u's from other voxels
				# which we can express easily because we have a global system
				for coeff in range(self.numCoeffs):

					# evaluate q (the RHS) for current coefficient at current voxel center
					q_lm = 0.0
					if coeff == 0:
						# we assume isotropic sources for now...all higher moments of q_lm are zero
						q_lm = source( pWS[0], pWS[1] )

					# evaluate phase function moment at current voxel center
					# todo support for spatial varying phase function moments
					# TODO: work out/verify lambda_l coefficient
					f_lm = 0.0
					if coeff == 0:
						# we assume isotropic phase function for now...all higher moments of f_lm are zero
						f_lm = 1.0

					ss = sigma_s( pWS[0], pWS[1] )
					st = sigma_a( pWS[0], pWS[1] )+ss


					# find the equation index within our global system
					i = self.get_global_index(voxel_i, voxel_j, coeff)

					# set the RHS in the global system
					# TODO: apply to real transformation, once we use non-isotropic sources
					# in the isotropic source case, the zero moment is real and is not touched
					# by the transformation in S. Therefore S has no effect when source is isotropic.
					b[i] = q_lm

					# now we execute the dot product between the current row in M (coeff)
					# and the solution vector u
					# k runs over the (non-zero) columns in M (identified by Ix/Iy) and the rows in u

					#'''
					# Mx
					for k in self.Ix[coeff]:
						M_ij = self.Mx_real[coeff, k]
						k_im1 = self.get_global_index(voxel_i-1, voxel_j, k)
						k_ip1 = self.get_global_index(voxel_i+1, voxel_j, k)
						
						if voxel_i > 0:
							A[i, k_im1] = -1.0/(2.0*self.domain.h_x)*M_ij
						if voxel_i < self.domain.res_x-1:
							A[i, k_ip1] = 1.0/(2.0*self.domain.h_x)*M_ij
					# My
					for k in self.Iy[coeff]:
						M_ij = self.My_real[coeff, k]
						k_jm1 = self.get_global_index(voxel_i, voxel_j-1, k)
						k_jp1 = self.get_global_index(voxel_i, voxel_j+1, k)
						if voxel_j > 0:
							A[i, k_jm1] = -1.0/(2.0*self.domain.h_y)*M_ij
						if voxel_j < self.domain.res_y-1:
							A[i, k_jp1] = 1.0/(2.0*self.domain.h_y)*M_ij
					#'''

					# C (is assumed to be diagonal, which implies isotropic medium and 1d phase function)
					# NB: this here is the problem why this method doesnt work in voids
					# iterative methods rely upon isolating the unknown of each equation on the left side,
					# which requires multiplying by the inverse of its coefficient, which in case
					# of vacuum would be a division by zero
					A[i, i] = st-ss*f_lm
					#A[i, i] = st
					#A[i, i] = -ss*f_lm

		data = {}
		data["A"] = A
		data["b"] = b.reshape((numVoxels*self.numCoeffs, 1))
		data["numCoeffs"] = self.numCoeffs
		#data["b"] = b
		#scipy.io.savemat("C:/projects/epfl/epfl17/python/simplepn/data.mat", data)
		scipy.io.savemat("C:/projects/epfl/epfl17/python/sopn/data_simplepn.mat", data)
		#alpha = 1.0
		#R = np.identity(numVoxels*self.numCoeffs) - alpha*A
		#ev = np.linalg.eigvals(R)
		#print(np.max(np.abs(ev)))


		# now we do simple richardson iteration








# define problem ---
def source( x, y ):
    if x > 3.0 and x < 4.0 and y > 3.0 and y < 4.0:
        return 1.0
    return 0.0

def sigma_a( x, y ):
    cx = np.ceil(x)
    cy = np.ceil(y)
    g = 0
    if np.ceil((x+y)/2.0)*2.0 == (cx+cy) and cx > 1.0 and cx < 7.0 and cy > 1.0 and cy-2.0*np.abs(cx-4.0) < 4:
        g = 1
    return (1.0-g)*0 + g*10

def sigma_s( x, y ):
    cx = np.ceil(x)
    cy = np.ceil(y)
    g = 0
    if np.ceil((x+y)/2.0)*2.0 == (cx+cy) and cx > 1.0 and cx < 7.0 and cy > 1.0 and cy-2.0*np.abs(cx-4.0) < 4:
        g = 1
    return (1.0-g)*1 + g*0

if __name__ == "__main__":

	order = 1
	#solver = SimplePN2D(order, util.Domain2D(7.0, 50))
	solver = SimplePN2D(order, util.Domain2D(7.0, 70))
	u = solver.run(sigma_a, sigma_s, source)
