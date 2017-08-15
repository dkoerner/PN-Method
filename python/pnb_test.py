import numpy as np
import pnb


offset = np.array([1,1])


voxels = np.zeros([3,3], dtype=complex)

voxels[2, 0] = 1.23j
voxels[0, 1] = 1.23

test = pnb.CoefficientGrid(voxels, offset)


voxel = np.array([0,1])
l = pnb.Location(voxel, offset)

result = test(l)

print(result)