import numpy as np
import pnb


offset = np.array([1,1])

resolution = np.array([3,3])
voxelsize = np.array([0.1,0.1])

voxels = np.zeros([resolution[0],resolution[1]], dtype=complex)

voxels[2, 0] = 1.23j
voxels[0, 1] = 1.23

domain = pnb.Domain( resolution, voxelsize )

test = pnb.CoefficientGrid(voxels, domain, offset)


voxel = np.array([0,1])
l = pnb.GridLocation(voxel, offset)

result = test(l)

print(result)