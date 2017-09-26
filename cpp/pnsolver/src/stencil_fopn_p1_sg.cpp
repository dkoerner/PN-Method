// This file was generated by stencil.py

#include <PNSystem.h>

void stencil_fopn_p1_sg(PNSystem::VoxelSystem& sys,
					PNSystem::Fields& fields)
{
	V2i vi = sys.getVoxel();
	V2d vd = sys.getVoxel().cast<double>();
	V2d h_inv( 1.0/(1*sys.getVoxelSize()[0]), 1.0/(1*sys.getVoxelSize()[1]) );

	Eigen::Matrix<std::complex<double>, 3, 3> S;
	S.coeffRef(0, 0) = std::complex<double>(1.0, 0.0);
	S.coeffRef(1, 1) = std::complex<double>(0.7071067811865475, 0.0);
	S.coeffRef(1, 2) = std::complex<double>(-0.7071067811865475, 0.0);
	S.coeffRef(2, 1) = std::complex<double>(-0.0, -0.7071067811865475);
	S.coeffRef(2, 2) = std::complex<double>(-0.0, -0.7071067811865475);
	Eigen::Matrix<std::complex<double>, 3, 3> SInv;
	SInv.coeffRef(0, 0) = std::complex<double>(1.0, 0.0);
	SInv.coeffRef(1, 1) = std::complex<double>(0.7071067811865476, 0.0);
	SInv.coeffRef(1, 2) = std::complex<double>(0.0, 0.7071067811865476);
	SInv.coeffRef(2, 1) = std::complex<double>(-0.7071067811865476, 0.0);
	SInv.coeffRef(2, 2) = std::complex<double>(-0.0, 0.7071067811865476);

	//Producing complex-valued matrices =============
	//M_0dxL + M_1dyL + M_2dzL + M_3L = b

	//M_0 ---
	// is constant
	Eigen::Matrix<double, 3, 3> M_0_real;
	M_0_real(0, 0) = 0.0;
	M_0_real(0, 1) = 0.5773502691896258;
	M_0_real(0, 2) = 0.0;
	M_0_real(1, 0) = 0.5773502691896257;
	M_0_real(1, 1) = 0.0;
	M_0_real(1, 2) = 0.0;
	M_0_real(2, 0) = 0.0;
	M_0_real(2, 1) = 0.0;
	M_0_real(2, 2) = 0.0;

	//M_1 ---
	// is constant
	Eigen::Matrix<double, 3, 3> M_1_real;
	M_1_real(0, 0) = 0.0;
	M_1_real(0, 1) = 0.0;
	M_1_real(0, 2) = 0.5773502691896258;
	M_1_real(1, 0) = 0.0;
	M_1_real(1, 1) = 0.0;
	M_1_real(1, 2) = 0.0;
	M_1_real(2, 0) = 0.5773502691896257;
	M_1_real(2, 1) = 0.0;
	M_1_real(2, 2) = 0.0;

	//M_2 ---
	// all components vanish

	//M_3 ---
	Eigen::Matrix<std::complex<double>, 3, 3> M_3;
	M_3(0, 0) = (fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(1, 1) = (fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))+
			-(fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5)))));
	M_3(2, 2) = (fields.sigma_t->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))+
			-(fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0)))));
	Eigen::Matrix<double, 3, 3> M_3_real = (S*M_3*SInv).real();

	//b ---
	Eigen::Matrix<std::complex<double>, 3, 1> b;
	b(0, 0) = fields.q->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(1, 0) = fields.q->eval(1, -1, sys.voxelToWorld(vd+V2d(0.0, 0.5)));
	b(2, 0) = fields.q->eval(1, 1, sys.voxelToWorld(vd+V2d(0.5, 0.0)));
	Eigen::Matrix<double, 3, 1> b_real = (S*b).real();

	// Assembling global system =============
	sys.coeff_A( 1, vi + V2i(-1,0), 0 ) += -(h_inv[0]*0.57735026919);
	sys.coeff_A( 1, vi + V2i(0,0), 0 ) += (h_inv[0]*0.57735026919);
	sys.coeff_A( 0, vi + V2i(0,0), 1 ) += -(h_inv[0]*0.57735026919);
	sys.coeff_A( 0, vi + V2i(1,0), 1 ) += (h_inv[0]*0.57735026919);
	sys.coeff_A( 2, vi + V2i(0,-1), 0 ) += -(h_inv[1]*0.57735026919);
	sys.coeff_A( 2, vi + V2i(0,0), 0 ) += (h_inv[1]*0.57735026919);
	sys.coeff_A( 0, vi + V2i(0,0), 2 ) += -(h_inv[1]*0.57735026919);
	sys.coeff_A( 0, vi + V2i(0,1), 2 ) += (h_inv[1]*0.57735026919);
	sys.coeff_A( 0, vi + V2i(0,0), 0 ) += M_3_real.coeffRef(0, 0);
	sys.coeff_A( 1, vi + V2i(0,0), 1 ) += M_3_real.coeffRef(1, 1);
	sys.coeff_A( 2, vi + V2i(0,0), 2 ) += M_3_real.coeffRef(2, 2);
	sys.coeff_b( 0 ) += b_real.coeffRef(0, 0);
	sys.coeff_b( 1 ) += b_real.coeffRef(1, 0);
	sys.coeff_b( 2 ) += b_real.coeffRef(2, 0);
}
REGISTER_STENCIL(stencil_fopn_p1_sg, 1)