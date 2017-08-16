

import numpy as np
import pnbuilder
import problems
import lspn
import util
import field
import scipy.io


class Test(object):
	def __init__( self, location, problem, term ):
		self.location = location
		self.problem = problem
		self.term = term
	def __call__( self, angle ):
		return self.term(self.location, angle[0], angle[1], self.problem)

if __name__ == "__main__":

	debug_voxel_x = 40
	debug_voxel_y = 35
	voxels_min = np.array([debug_voxel_x, debug_voxel_y])
	voxels_max = np.array([debug_voxel_x+1, debug_voxel_y+1])

	#voxels_min = np.array([20, 20])
	#voxels_max = np.array([50, 50])

	id = "checkerboard_blur10.0_term1"
	data = util.load_pn_solution("C:/projects/epfl/epfl17/python/sopn/solution_{}.mat".format(id))
	pnb = data["pnb"]
	x_real = data["x_real"]

	# here we setup the radiance field from the pn solution
	#problem = {}
	problem = problems.blurred(problems.checkerboard2d(), 10.0)
	problem["L"] = field.from_pn_solution(pnb.to_complex(x_real), pnb)

	terms = []
	terms.append((0, lspn.term0))
	#terms.append((1, lspn.term1))


	data = {}
	for (term_index, term) in terms:
		#problem = problems.blurred(problems.checkerboard2d(), 10.0)

		# get the matrix A for the current term
		A_term = scipy.io.loadmat("C:/projects/epfl/epfl17/python/sopn/debug_terms/debug_A_term{}.mat".format(term_index))["A"]
		x_term_real = A_term.dot(x_real)

		x_term_complex_gt = np.zeros(pnb.domain_cpp.numVoxels()*pnb.num_coeffs(), dtype=complex)
		for voxel_i in range(voxels_min[0], voxels_max[0]):
			#print("voxel_i={}".format(voxel_i))
			for voxel_j in range(voxels_min[1], voxels_max[1]):
				print("voxel_i={} voxel_j={}".format(voxel_i, voxel_j))
				# here we project the term under investigation into SH. This is done at the respective locations
				# of all the unknowns
				#for i in range(pnb.num_coeffs()):
				for i in range(1):
					(l,m) = pnb.lm_index(i)
					location = pnb.unknown_location(voxel_i, voxel_j, i)
					pWS = location.getPWS()
					global_i = pnb.global_index(voxel_i, voxel_j, i)
					# NB: we take into account, that for 2d, pnb will have different index and lm ordering
					#print("---------------")
					#gg = lambda theta, phi: term(pWS, theta, phi, problem)
					#print(gg(np.pi*0.3, 0.0))
					coeff = util.project_sh_coeff(lambda theta, phi: term(pWS, theta, phi, problem), l, m)
					#coeff = shtools.project_sh_coeff_x(lambda angle: term(pWS, angle[0], angle[1], problem), l, m)
					#coeff = util.project_sh_coeff(Test(location, problem, term), l, m)
					#print("done")
					x_term_complex_gt[global_i] = coeff
		x_term_real_gt = pnb.to_real(x_term_complex_gt)

		debug_i = pnb.global_index(debug_voxel_x,debug_voxel_y,0)
		print(x_term_real[debug_i])
		print(x_term_real_gt[debug_i])
		#data['x_term{}'.format(term_index)] = x_term_real_gt.reshape((pnb.domain.numVoxels*pnb.numCoeffs, 1))
		#scipy.io.savemat("C:/projects/epfl/epfl17/python/pnb/debug_terms/data_term{}.mat".format(term_index), data)


