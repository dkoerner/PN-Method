# This module uses the pnsolver to compute the solution of the checkerboard problem.
# It assumes that the cpp stencil had been generated and compiled with the pnsolver module.

import numpy as np
import problems
import util
import pnsolver


if __name__ == "__main__":
	# define problem ----------------------------
	problem = problems.checkerboard()

	# initialize the solver -----------------------------------
	# NB: the truncation order has been baked into the cpp code as it is linked directly with the stencil
	sys = pnsolver.PNSystem(problem["domain"])

	#index = 6392
	#(voxel, coeff) = sys.getVoxelAndCoefficient(index)
	#(l,m) = sys.getLM(coeff)
	#print("index={} voxel={} {} coeff={} (l={} m={})".format(index, voxel[0], voxel[1], coeff, l, m))
	#sys.setDebugVoxel(voxel)


	# set the RTE parameter fields ---------------------------------------------
	sys.setField( "sigma_t", problem["sigma_t"] )
	sys.setField( "sigma_a", problem["sigma_a"] )
	sys.setField( "sigma_s", problem["sigma_s"] )

	# if only one coefficient field is given for f_p, we assume its meant to
	# define the zero coefficient field
	if len(problem["f_p"]) == 1:
		sys.setField( "f_p", problem["f_p"][0] )
	else:
		# TODO: set f_p field for given shindex (or l,m)
		raise ValueError("no higher order phase functions supported yet")

	# if only one coefficient field is given for q, we assume its meant to
	# define the zero coefficient field
	if len(problem["q"]) == 1:
		sys.setField( "q", problem["q"][0] )
	else:
		# TODO: set f_p field for given shindex (or l,m)
		raise ValueError("no higher order source functions supported yet")


	# build the system Ax=b using the cpp stencil -------------------
	sys.build()

	# solve for x -------------------
	x = None
	x = sys.solve()

	# write the result to disk for visualization ---------------
	path = "C:/projects/epfl/epfl17/python/pnsolver/results/terms_new"
	postfix = ""
	#postfix = "_nonrasterized"
	#postfix = "_not3"
	filename = "{}/{}{}.mat".format(path, problem["id"], postfix)
	util.write_pn_system(filename, sys, problem, x)








