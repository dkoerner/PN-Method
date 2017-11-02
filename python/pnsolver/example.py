# This module uses the pnsolver to compute the solution of the checkerboard problem.
# It assumes that the cpp stencil had been generated and compiled with the pnsolver module.

import numpy as np
import problems
import util
import pnsolver
import itertools
import scipy.io




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
	sys = pnsolver.PNSystem(stencil_name, problem["domain"], do_neumannBC)

	# set the RTE parameter fields ---------------------------------------------
	setProblem(sys, problem)
	#sys.setNeumannBoundaryConditions(do_neumannBC)

	# build the system Ax=b using the cpp stencil -------------------
	# build is now done in the solve functions
	#sys.build()

	#A = sys.get_boundary_A_real()
	#b = sys.get_boundary_b_real()

	# solve for x -------------------
	x = None

	#try:
	#x = sys.solve()
	#b = sys.get_b_real_test()
	#x = sys.solve_cg( b )
	#x = sys.solve_mg()
	#x, convergence, timestamps = pnsolver.solve_sparseLU( sys )
	#x = x.reshape((x.shape[0], 1))


	data = {}

	numLevels = 7
	#x, convergence, timestamps = pnsolver.solve_multigrid( sys, numLevels )
	#data["convergence_mg"] = convergence
	#data["timestamps_mg"] = timestamps


	#x, convergence, timestamps = pnsolver.solve_multigrid2( sys, numLevels )
	#data["convergence_mg2"] = convergence
	#data["timestamps_mg2"] = timestamps
	#data["x_mg"] = x.reshape((x.shape[0], 1))

	#x, convergence, timestamps = pnsolver.solve_gs( sys )
	#data["convergence_gs"] = convergence
	#data["timestamps_gs"] = timestamps

	#x, convergence, timestamps = pnsolver.solve_blockgs( sys )
	#data["convergence_mg2"] = convergence
	#data["timestamps_mg2"] = timestamps


	x, convergence, timestamps = pnsolver.solve_lscg( sys )
	data["convergence_cg"] = convergence
	data["timestamps_cg"] = timestamps


	pnsolver.save_solution("test.pns", sys, x)
	pnsolver.save_fields("test.fields", sys) # stores the voxelgrids of the RTE parameters
	#x_test = pnsolver.getSolutionVector(pnsolver.load_solution("test.pns"))
	#util.compare_vectors( x.reshape((x.shape[0], 1)), x_test.reshape((x.shape[0], 1)), "org", "new" )

	data["x"] = x.reshape((x.shape[0], 1))



	#data["x_cg"] = x.reshape((x.shape[0], 1))

	#x, convergence, timestamps = pnsolver.solve_cg_eigen( sys )
	#data["convergence_cg_eigen"] = convergence
	#data["timestamps_cg_eigen"] = timestamps

	#x, convergence, timestamps = pnsolver.solve_sparseLU( sys )
	#data["convergence_cg_eigen"] = convergence
	#data["timestamps_cg_eigen"] = timestamps


	data["resolution"] = sys.getResolution()
	data["numCoeffs"] = sys.getNumCoefficients()

	#a,b,c = sys.get_debug();
	#data["debug_x"] = a.reshape((a.shape[0], 1))
	#data["debug_x_upsampled_downsampled"] = c.reshape((c.shape[0], 1))

	#filename = "c:/projects/epfl/temp/multigrid/multigrid-master/checkerboard_test.mat"
	#filename = "c:/projects/epfl/temp/multigrid/multigrid-master/homogeneous_test.mat"
	#scipy.io.savemat(filename, data)

	#pass
	#except:
	#	print("Solve failed...")

	# write the result to disk for visualization ---------------
	util.write_pn_system(filename, sys, problem, x, convergence, timestamps)

def groundtruth( problem, numSamples, filename ):
	stencil_name = "noop"
	sys = pnsolver.PNSystem(stencil_name, problem["domain"], False)

	# set the RTE parameter fields ---------------------------------------------
	setProblem(sys, problem)

	# compute groundtruth ---------------------
	x = sys.computeGroundtruth(numSamples);

	# xport ---------------------
	util.write_pn_system(filename, sys, problem, x)



def memory():
	import os
	import psutil
	process = psutil.Process(os.getpid())
	return process.memory_info().rss

def multigrid_debug_test():
	#stencil_name = "noop"
	#domain = pnsolver.Domain( np.array([1.0, 1.0, 1]), np.array([32, 32, 1]), np.array([0.0, 0.0, 0.0]))
	#do_neumannBC = False
	#sys = pnsolver.PNSystem(stencil_name, domain, do_neumannBC)
	
	numLevels = 9
	mgtest = pnsolver.Solver(numLevels)

	#'''
	for l in range(numLevels):
		filename = "c:/projects/epfl/temp/multigrid/multigrid-master/poisson_problem_lvl{}.mat".format(l)
		data = scipy.io.loadmat(filename)
		mgtest.setMultigridLevel(l, data["L"], data["u"], data['restrict_m'], data['interp_m'])
	#'''


	filename = "c:/projects/epfl/temp/multigrid/multigrid-master/poisson_problem_ref.mat"
	data = scipy.io.loadmat(filename)
	#mgtest.setRef( data["u_exact"] )
	mgtest.setb( data["f_gg"] )

	data = {}

	x, convergence, timestamps = mgtest.solve(100)
	data["convergence_mg"] = convergence
	data["timestamps_mg"] = timestamps
	x, convergence, timestamps = mgtest.solve_gs()
	data["convergence_gs"] = convergence
	data["timestamps_gs"] = timestamps

	x, convergence, timestamps = mgtest.solve_cg()
	data["convergence_cg"] = convergence
	data["timestamps_cg"] = timestamps


	filename = "c:/projects/epfl/temp/multigrid/multigrid-master/poisson_test.mat"
	scipy.io.savemat(filename, data)

	#pnsolver.solve_multigrid(sys)
	'''
	

	'''





if __name__ == "__main__":


	#multigrid_debug_test()
	#exit(0)

	path = "C:/projects/epfl/epfl17/python/pnsolver/results/studies"

	# define problem ----------------------------
	#problem = problems.checkerboard()
	#problem = problems.homogeneous()
	#util.write_problem(path+"/checkerboard_problem.mat", problem)
	#exit(1)
	problem = problems.checkerboard3d()
	#problem = problems.pointsource3d()
	#problem = problems.vacuum()

	
	#filename = "{}/{}_groundtruth.mat".format(path, problem["id"] )
	#groundtruth( problem, 8000000, filename )


	#'''
	staggered_id = {True:"sg", False:"cg"}
	bc_id = {True:"_nbc", False:""}


	#rte_forms = ["fopn", "sopn"]
	#order = [0,1,2,3,4,5]
	#order = [0,1,2,3,4]
	#staggered = [True, False]
	#boundary_conditions = [True, False]

	
	#'''
	#rte_forms = ["fopn"]
	#order = [1]
	#staggered = [True]
	#boundary_conditions = [False]


	rte_forms = ["fopn"]
	order = [1]
	staggered = [True]
	boundary_conditions = [False]


	#rte_forms = ["fopn"]
	#order = [1,2,3,4,5]
	#staggered = [True]
	#boundary_conditions = [False]

	#rte_forms = []
	#order = []
	#staggered = []
	#boundary_conditions = []

	test = itertools.product(rte_forms, order, staggered, boundary_conditions)
	for c in test:
		rte_form_name = c[0]
		order = c[1]
		is_staggered = c[2]
		do_neumannBC = c[3]

		stencil_name = "stencil_{}_p{}_{}".format(rte_form_name, order, staggered_id[is_staggered] )
		#stencil_name = "noop"
		#print(stencil_name)
		
		#filename = "{}/{}_{}{}.mat".format(path, problem["id"], stencil_name, bc_id[do_neumannBC])
		filename = "{}/{}_test.mat".format(path, problem["id"], stencil_name)
		#print("clear;filename=\"{}\";compute_condest;".format(filename))
		solve(stencil_name, problem, filename, do_neumannBC=do_neumannBC)

	#stencil_name = "stencil_cda"
	#do_neumannBC = False
	#filename = "{}/{}_{}.mat".format(path, problem["id"], stencil_name)
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
	filename_a = "c:/projects/epfl/epfl17/python/pnsolver/results/studies/checkerboard_test.mat"
	filename_b = "c:/projects/epfl/epfl17/python/pnsolver/results/studies/checkerboard_test_org.mat"
	data_a = util.load_pn_system(filename_a)
	data_b = util.load_pn_system(filename_b)
	A_a = data_a["A"]
	A_b = data_b["A"]
	b_a = data_a["b"]
	b_b = data_b["b"]

	#print("comparing A ----")
	#util.compare_matrices( A_a, A_b, "new", "original" )
	#gg = 2131
	#print(A_a[gg,gg])
	#print(A_b[gg,gg])
	#print("comparing b ----")
	#util.compare_matrices( b_a, b_b, "new", "original" )
	#print(b_a[6603,0])
	#print(b_b[6603,0])

	'''








