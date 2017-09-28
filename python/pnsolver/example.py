# This module uses the pnsolver to compute the solution of the checkerboard problem.
# It assumes that the cpp stencil had been generated and compiled with the pnsolver module.

import numpy as np
import problems
import util
import pnsolver
import itertools


def setProblem( sys, problem ):
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


def solve( stencil_name, problem, filename, do_neumannBC = False ):
	# initialize the solver -----------------------------------
	# NB: the truncation order has been baked into the cpp code as it is linked directly with the stencil
	sys = pnsolver.PNSystem(stencil_name, problem["domain"])

	# set the RTE parameter fields ---------------------------------------------
	setProblem(sys, problem)
	sys.setNeumannBoundaryConditions(do_neumannBC)

	# build the system Ax=b using the cpp stencil -------------------
	sys.build()

	#A = sys.get_boundary_A_real()
	#b = sys.get_boundary_b_real()

	# solve for x -------------------
	x = None

	try:
		#x = sys.solve()
		x = sys.solve_boundary()

	except:
		print("Solve failed...")

	# write the result to disk for visualization ---------------
	util.write_pn_system(filename, sys, problem, x)

def groundtruth( problem, numSamples, filename ):
	stencil_name = "noop"
	sys = pnsolver.PNSystem(stencil_name, problem["domain"])

	# set the RTE parameter fields ---------------------------------------------
	setProblem(sys, problem)

	# compute groundtruth ---------------------
	x = sys.computeGroundtruth(numSamples);

	# xport ---------------------
	util.write_pn_system(filename, sys, problem, x)





if __name__ == "__main__":
	path = "C:/projects/epfl/epfl17/python/pnsolver/results/studies"

	# define problem ----------------------------
	problem = problems.checkerboard()
	#problem = problems.vacuum()

	
	#filename = "{}/{}_groundtruth.mat".format(path, problem["id"] )
	#groundtruth( problem, 8000000, filename )


	#'''
	staggered_id = {True:"sg", False:"cg"}
	bc_id = {True:"_nbc", False:""}


	rte_forms = ["fopn", "sopn"]
	order = [0,1,2,3,4,5]
	staggered = [True, False]
	boundary_conditions = [True, False]

	#rte_forms = ["fopn"]
	#order = [1]
	#staggered = [False]
	#boundary_conditions = [False]

	test = itertools.product(rte_forms, order, staggered, boundary_conditions)
	for c in test:
		rte_form_name = c[0]
		order = c[1]
		is_staggered = c[2]
		do_neumannBC = c[3]

		stencil_name = "stencil_{}_p{}_{}".format(rte_form_name, order, staggered_id[is_staggered] )
		
		filename = "{}/{}{}{}.mat".format(path, problem["id"], stencil_name, bc_id[do_neumannBC])
		#filename = "{}/{}_test.mat".format(path, problem["id"], stencil_name)
		print("clear;filename=\"{}\";compute_condest;".format(filename))
		#solve(stencil_name, problem, filename, do_neumannBC=do_neumannBC)
	#'''


	#postfix = ""
	#postfix = "_allvacuum"
	#postfix = "_sigmas=0"
	#postfix = "_sigmat=const"
	#postfix = "_blurred"
	#postfix = ""
	#postfix = "_nonrasterized"
	#postfix = "_sopn_cg"
	#postfix = "_sopn_sg_x0"
	#postfix = "_sopn_sg_vacuum2"
	#postfix = "_sopn_sg_u"
	#postfix = "_fopn_sg_u"
	#postfix = "_p1_complex_sg"

	#postfix = "_sopn_p0_sg"
	#postfix = "_sopn_p1_sg"
	#postfix = "_sopn_p2_sg"
	#postfix = "_sopn_p3_sg"

	#postfix = "_fopn_p0_sg"
	#postfix = "_fopn_p1_sg"
	#postfix = "_fopn_p2_sg"
	#postfix = "_fopn_p3_sg"
	#postfix = "_fopn_p4_sg"
	#postfix = "_fopn_p5_sg"

	#postfix = "_fopn_p1_sg1"
	#postfix = "_fopn_p2_sg1"
	#postfix = "_fopn_p3_sg1"
	#postfix = "_fopn_p4_sg1"

	#postfix = "_fopn_p1_cg"
	#postfix = "_fopn_p2_cg"
	#postfix = "_fopn_p3_cg"
	#postfix = "_fopn_p4_cg"
	#postfix = "_fopn_p5_cg"


	#postfix = "_fopn_p3_sg1"
	#postfix = "_fopn_p1"
	#postfix = "_fopn_p2"
	#postfix = "_fopn_p3"
	#postfix = "_fopn_p4"
	#postfix = "_fopn_p5"
	#postfix = "_sopn_p1"
	#postfix = "_sopn_p1_not13"
	#postfix = "_sopn_p2"
	#postfix = "_sopn_p3"
	#postfix = "_cg_u"
	
	

	'''
	# compare matrices A
	#filename = "c:/projects/epfl/epfl17/python/pnsolver/results/terms_original/checkerboard_sopn_cg_term1.mat"
	filename = "c:/projects/epfl/epfl17/python/pnsolver/results/terms_original/checkerboard_sopn_cg.mat"
	data = util.load_pn_system(filename)
	print("comparing A ----")
	util.compare_matrices( sys.get_A_real(), data["A"], "new", "original" )
	print("comparing b ----")
	util.compare_matrices( sys.get_b_real(), data["b"], "new", "original" )
	'''








