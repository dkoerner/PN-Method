import numpy as np





gridoffset_fine = np.array([1,1,0])
gridoffset_coarse = np.array([2,2,0])


grids_3d = []
grids_3d.append((0, 0, 1))
grids_3d.append((1, 0, 1))
grids_3d.append((1, 1, 1))
grids_3d.append((0, 1, 1))
grids_3d.append((0, 0, 0))
grids_3d.append((1, 0, 0))
grids_3d.append((1, 1, 0))
grids_3d.append((0, 1, 0))

grids_2d = []
grids_2d.append((0, 0, 0))
grids_2d.append((1, 0, 0))
grids_2d.append((1, 1, 0))
grids_2d.append((0, 1, 0))


class Voxel(object):
	def __init__(self, coord, grid=None, parent=None):
		self.coord = coord
		self.grid = grid
		self.parent = parent
		#self.flags = []
	def createSubVoxels(self):
		self.sub_voxels = {}
		for a in range(2):
			for b in range(2):
				for c in range(2):
					coord = (self.coord[0]*2+a, self.coord[1]*2+b, self.coord[2]*2+c)
					grid = (a, b, c)
					v = Voxel( coord, grid, self )
					self.sub_voxels[grid] = v
	def __str__(self):
		string = ""
		string += "Voxel {} {} {}".format(self.coord[0], self.coord[1], self.coord[2])
		#string += "\t"
		return string
	def test(self, grid, level):
		# returns coordinate of the finest voxel which holds the grid location
		c = self.sub_voxels[grid].coord
		m = 2**level
		return (c[0]*m, c[1]*m, c[2]*m)

		





# construct a 3x3 neighbourhood of coarse voxels and its 2level subdivision
grid_coarse = {}
voxel_grid = {}
for i_coarse in range(-1,2):
	for j_coarse in range(-1,2):
		for k_coarse in range(-1,2): # 3d
		#for k_coarse in range(1): # 2d
			coord_coarse = (i_coarse, j_coarse, k_coarse)
			v_coarse = Voxel(coord_coarse)
			grid_coarse[v_coarse.coord] = v_coarse
			v_coarse.createSubVoxels()
			# create second level as well
			for coord, voxel in v_coarse.sub_voxels.items():
				voxel.createSubVoxels()
				# and construct our finest level grid lookup table
				for coord2, voxel2 in voxel.sub_voxels.items():
					voxel_grid[voxel2.coord] = voxel2




v_coarse = grid_coarse[(0,0,0)]

print( "std::vector<std::vector<std::pair<V3i, double>>> stamps(8);" )

grids = grids_3d # 3d
#grids = grids_2d # 2d


for grid_index in range(len(grids)):
	grid = grids[grid_index]
	v = v_coarse.sub_voxels[grid]
	center = v_coarse.test(grid, 1)

	# now search within a 3x3x3 search radius
	for offset_i in range(-3, 4):
		for offset_j in range(-3, 4):
			for offset_k in range(-3, 4): #3d
			#for offset_k in range(1): # 2d
				offset = (offset_i, offset_j, offset_k)
				coord = (center[0]+offset_i, center[1]+offset_j, center[2]+offset_k)
				v = voxel_grid[coord]
				# check if the component we are looking for is located
				# on the finest grid
				if v.grid == grid:
					# it is. find the coordinate/offset of the fine voxel
					offset_fine = v.parent.coord
					weight = [1,1,1]
					for i in range(3):
						weight[i] = [1.0, 0.75, 0.5, 0.25][np.abs(offset[i])]

					#print( "offset={} {} {}  weight={}*{}*{}".format(offset_fine[0], offset_fine[1], offset_fine[2], weight[0], weight[1], weight[2]) )
					print("stamps[{}].push_back(std::make_pair(V3i({}, {}, {}), {}*{}*{}));".format(grid_index, offset_fine[0], offset_fine[1], offset_fine[2], weight[0], weight[1], weight[2]))
