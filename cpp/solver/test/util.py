import numpy as np

def pointsource_2d_fluence( pWS, center = np.array([0.0, 0.0]), power = 1.0 ):
	''' returns fluence at position pWS for a unit power point light located at origin'''
	r = np.linalg.norm(pWS-center)
	if r == 0.0:
		return power/(2.0*np.pi*max(r, 1.0e-4))
	return power/(2.0*np.pi*r)

class Domain2D:
	def __init__(self, size, res):
		self.res = res
		self.size = size
		self.h = size/float(res)
		self.bound_min = np.array([-size*0.5, -size*0.5])
		self.bound_max = np.array([size*0.5, size*0.5])
		#self.voxelToWorldTransform = 

	def voxelToLocal( self, pVS ):
		return np.array([pVS[0]/self.res, pVS[1]/self.res])

	def localToVoxel(self, pLS):
		return np.array([pLS[0]*self.res, pLS[1]*self.res])

	def localToWorld( self, pLS ):
		return np.array([pLS[0]*self.size + self.bound_min[0], pLS[1]*self.size + self.bound_min[1]])

	def worldToLocal( self, pWS ):
		return np.array([(pWS[0]-self.bound_min[0])/self.size, (pWS[1]-self.bound_min[1])/self.size])

	def voxelToWorld(self, pVS):
		return self.localToWorld(self.voxelToLocal(pVS))

	def worldToVoxel(self, pWS):
		return self.localToVoxel(self.worldToLocal(pWS))

	def voxelToIndex(self, pVS):
		return (int(pVS[1]), int(pVS[0]))

	def rasterize( self, fun ):
		field = np.zeros((self.res, self.res))
		for i in range(0, self.res):
			for j in range(0, self.res):
				pVS = np.array([j+0.5, i+0.5])
				pWS = self.voxelToWorld(pVS)
				field[i, j] = fun(pWS)
		return field

	def rasterizeVS( self, fun ):
		field = np.zeros((self.res, self.res))
		for i in range(0, self.res):
			for j in range(0, self.res):
				field[i, j] = fun(i, j)
		return field

	def gradient( self, field ):
		grad_field = np.zeros((self.res, self.res, 2))
		for i in range(0, self.res):
			for j in range(0, self.res):
				if j==0:
					# forward differences for x
					grad_field[i,j,0] = (field[i,j+1]-field[i,j])/self.h
				elif j==self.res-1:
					#backward differences for x
					grad_field[i,j,0] = (field[i,j]-field[i,j-1])/self.h
				else:
					# central differences for x
					grad_field[i,j,0] = (field[i,j+1]-field[i,j-1])/(self.h*2.0)

				if i==0:
					# forward differences for y
					grad_field[i,j,1] = (field[i,j]-field[i+1,j])/self.h
				elif i==self.res-1:
					#backward differences for y
					grad_field[i,j,1] = (field[i-1,j]-field[i,j])/self.h
				else:
					# central differences for y
					grad_field[i,j,1] = (field[i-1,j]-field[i+1,j])/(self.h*2.0)
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
				elif j==self.res-1:
					#backward differences for x
					grad_field[i] = (field[i]-field[i-1])/self.h
				else:
					# central differences for x
					grad_field[i] = (field[i+1]-field[i-1])/(self.h*2.0)
		return grad_field
