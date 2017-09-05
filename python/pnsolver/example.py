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
	x = sys.solve()

	# write the result to disk for visualization ---------------
	util.write_pn_system(sys, problem, "pns_highres2", x)








