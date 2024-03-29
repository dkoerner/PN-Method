// This file was generated by stencil.py

#include <PNSystem.h>

void stencil_fopn_p1_cg(PNSystem::Stencil::Context& ctx)
{
	V3i vi = ctx.getVoxelCoord();
	V3d vd = vi.cast<double>();
	const Domain& domain = ctx.getDomain();
	const PNVolume& problem = ctx.getProblem();
	V3d h_inv( 1.0/(2*domain.getVoxelSize()[0]), 1.0/(2*domain.getVoxelSize()[1]), 1.0/(2*domain.getVoxelSize()[2]) );
	int color_channel = 0;

	Eigen::Matrix<std::complex<double>, 3, 3> S;
	S.coeffRef(0, 0) = std::complex<double>(1.0, 0.0);
	S.coeffRef(1, 1) = std::complex<double>(0.0, 0.7071067811865475);
	S.coeffRef(1, 2) = std::complex<double>(0.0, 0.7071067811865475);
	S.coeffRef(2, 1) = std::complex<double>(0.7071067811865475, 0.0);
	S.coeffRef(2, 2) = std::complex<double>(-0.7071067811865475, 0.0);
	Eigen::Matrix<std::complex<double>, 3, 3> SInv;
	SInv.coeffRef(0, 0) = std::complex<double>(1.0, 0.0);
	SInv.coeffRef(1, 1) = std::complex<double>(0.0, -0.7071067811865476);
	SInv.coeffRef(1, 2) = std::complex<double>(0.7071067811865476, 0.0);
	SInv.coeffRef(2, 1) = std::complex<double>(0.0, -0.7071067811865476);
	SInv.coeffRef(2, 2) = std::complex<double>(-0.7071067811865476, -0.0);

	//Producing complex-valued matrices =============
	//M_0dxL + M_1dyL + M_2dzL

	//M_0 ---
	// is constant

	//M_1 ---
	// is constant

	//M_2 ---
	// all components vanish

	// Assembling global system =============
	ctx.coeff_A( 2, vi + V3i(-1,0,0), 0 ) += -(h_inv[0]*0.57735026919);
	ctx.coeff_A( 2, vi + V3i(1,0,0), 0 ) += (h_inv[0]*0.57735026919);
	ctx.coeff_A( 0, vi + V3i(-1,0,0), 2 ) += -(h_inv[0]*0.57735026919);
	ctx.coeff_A( 0, vi + V3i(1,0,0), 2 ) += (h_inv[0]*0.57735026919);
	ctx.coeff_A( 1, vi + V3i(0,-1,0), 0 ) += -(h_inv[1]*-0.57735026919);
	ctx.coeff_A( 1, vi + V3i(0,1,0), 0 ) += (h_inv[1]*-0.57735026919);
	ctx.coeff_A( 0, vi + V3i(0,-1,0), 1 ) += -(h_inv[1]*-0.57735026919);
	ctx.coeff_A( 0, vi + V3i(0,1,0), 1 ) += (h_inv[1]*-0.57735026919);
}
V3i stencil_fopn_p1_cg_get_offset(int coeff)
{
	switch(coeff)
	{
		case 0:return V3i(1, 1, 1);break;
		case 1:return V3i(1, 1, 1);break;
		case 2:return V3i(1, 1, 1);break;
		default:throw std::runtime_error("unexpected coefficient index");break;
	};
}
REGISTER_STENCIL(stencil_fopn_p1_cg, 1, 3, 1)
