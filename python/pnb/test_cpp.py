import numpy as np
import pnbuilder as pnb


size = np.array([7.0, 7.0])
resolution = np.array([70, 70])
offset = np.array([0.0, 0.0])
domain = pnb.Domain( size, resolution, offset )

voxel = np.array([0,0])
offset = np.array([0,0])
l = pnb.GridLocation( domain, voxel, offset )
'''
t = l.getOffset()

def check(t):
	print(type(t))
	print(t.shape)
	print(t.dtype)
	print(t)

check(size)
check(t)
'''
a = -1
b= 2

l.test( a, b )


print(divmod( a, b ))
print(a//b)