// This file was generated by stencil.py

#include <PNSystem.h>

void stencil_fopn_p0_cg(PNSystem::VoxelSystem& sys,
					PNSystem::Fields& fields)
{
	V2i vi = sys.getVoxel();
	V2d vd = sys.getVoxel().cast<double>();
	V2d h_inv( 1.0/(2*sys.getVoxelSize()[0]), 1.0/(2*sys.getVoxelSize()[1]) );

	Eigen::Matrix<std::complex<double>, 1, 1> S;
	S.coeffRef(0, 0) = std::complex<double>(1.0, 0.0);
	Eigen::Matrix<std::complex<double>, 1, 1> SInv;
	SInv.coeffRef(0, 0) = std::complex<double>(1.0, 0.0);

	//Producing complex-valued matrices =============
	//M_0dxL + M_1dyL + M_2dzL + M_3L = b

	//M_0 ---
	// all components vanish

	//M_1 ---
	// all components vanish

	//M_2 ---
	// all components vanish

	//M_3 ---
	Eigen::Matrix<std::complex<double>, 1, 1> M_3;
	M_3(0, 0) = (fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5)))));
	Eigen::Matrix<double, 1, 1> M_3_real = (S*M_3*SInv).real();

	//b ---
	Eigen::Matrix<std::complex<double>, 1, 1> b;
	b(0, 0) = fields.q->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5)));
	Eigen::Matrix<double, 1, 1> b_real = (S*b).real();

	// Assembling global system =============
	sys.coeff_A( 0, vi + V2i(0,0), 0 ) += M_3_real.coeffRef(0, 0);
	sys.coeff_b( 0 ) += b_real.coeffRef(0, 0);
}
REGISTER_STENCIL(stencil_fopn_p0_cg, 0, 0)
