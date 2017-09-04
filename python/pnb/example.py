# This module represents the application which uses pnbuilder to solve some problem.

import numpy as np
import problems
import pnbuilder
import util
import lspn
import meh
import itertools
import pnsolver
import scipy.io
import stencil



def print_matrix( A ):
	for i in range(A.shape[0]):
		for j in range(A.shape[1]):
			v = A[i,j]
			if np.abs(np.real(v)) > 1.0e-8 or np.abs(np.imag(v)) > 1.0e-8:
				print("  ({}, {}) {}".format(i, j, v))



def compare( A0, A1, id0, id1 ):
	if A0.shape != A1.shape:
		raise ValueError("unmatched shape")

	if A0.dtype != A1.dtype:
		raise ValueError("unmatched dtype")

	diff = A0-A1
	abs_real_diff = np.abs(np.real(diff))
	abs_imag_diff = np.abs(np.imag(diff))
	#diff_real_flatten = abs_real_diff.flatten()
	#diff_imag_flatten = abs_imag_diff.flatten()

	max_index = np.unravel_index(np.argmax(abs_real_diff), abs_real_diff.shape)
	mm = abs_real_diff[max_index]
	print("real: max={} max_index={} {}".format(mm,max_index[0], max_index[1]))
	max_index = np.unravel_index(np.argmax(abs_imag_diff), abs_imag_diff.shape)
	mm = abs_imag_diff[max_index]
	print("imag: max={} max_index={} {}".format(mm,max_index[0], max_index[1]))


	#su = diff_real_flatten.sum(axis=1)
	#print("real: sum={}".format(su[0,0]))
	#mm = diff_imag_flatten.max(axis=1)
	#su = diff_imag_flatten.sum(axis=1)

	#print("imag: max={}".format(mm[0,0]))
	#print("imag: sum={}".format(su[0,0]))

	'''
	print("{} result:".format(id0))
	for i in range(A0.shape[0]):
		for j in range(A0.shape[1]):
			v = A0[i,j]
			if np.abs(np.real(v)) > 1.0e-8 or np.abs(np.imag(v)) > 1.0e-8:
				print("  ({}, {}) {}".format(i, j, v))


	print("{} result:".format(id1))
	for i in range(A1.shape[0]):
		for j in range(A1.shape[1]):
			v = A1[i,j]
			if np.abs(np.real(v)) > 1.0e-8 or np.abs(np.imag(v)) > 1.0e-8:
				print("  ({}, {}) {}".format(i, j, v))
	'''



def test_pnsystem():
	pass

	#global_index=2131 -> voxel=10 10  coeff=1
	#global_index = 2131
	#(voxel, coeff) = sys.getVoxelAndCoefficient(global_index)
	#print( "global_index={} -> voxel={} {}  coeff={}".format(global_index, voxel[0], voxel[1], coeff) )
	#return


	#debugVoxel = np.array([10,10])
	#index_i = sys.getGlobalIndex(debugVoxel, 0)
	#sys.setDebugVoxel(debugVoxel)


	'''
	
	#return

	A_sys = sys.get_A_real()
	b_sys = sys.get_b_real()
	#A_sys = sys.get_A()
	#b_sys = sys.get_b()


	# load 
	filename = "C:/projects/epfl/epfl17/python/sopn/system2_{}.mat".format("checkerboard")
	data = scipy.io.loadmat(filename)
	A_builder = data["A"]#[index_i:index_i+sys.getNumCoefficients(),:]
	b_builder = data["b"]#[index_i:index_i+sys.getNumCoefficients(),:]

	if 'index_i' in locals():
		A_sys = A_sys[index_i:index_i+sys.getNumCoefficients(),:]
		b_sys = b_sys[index_i:index_i+sys.getNumCoefficients(),:]
		A_builder = A_builder[index_i:index_i+sys.getNumCoefficients(),:]
		b_builder = b_builder[index_i:index_i+sys.getNumCoefficients(),:]

	print("A ----------------------")
	compare(A_sys, A_builder, "PNSystem", "PNBuilder")
	print("b ----------------------")
	compare(b_sys, b_builder, "PNSystem", "PNBuilder")
	'''



if __name__ == "__main__":
	# define problem ----------------------------
	problem = problems.checkerboard()

	# initialize the solver -----------------------------------
	# NB: the truncation order has been baked into the cpp code as it is linked directly with the stencil
	sys = pnsolver.PNSystem(problem["domain_cpp"])

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








