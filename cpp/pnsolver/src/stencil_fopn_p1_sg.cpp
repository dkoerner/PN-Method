// This file was generated by stencil.py

#include <PNSystem.h>

void stencil_fopn_p1_sg(PNSystem::Stencil::Context& ctx)
{
	V2i vi = ctx.getVoxel();
	V2d vd = vi.cast<double>();
	const Domain& domain = ctx.getDomain();
	const PNSystem::Fields& fields = ctx.getFields();
	V2d h_inv( 1.0/(1*domain.getVoxelSize()[0]), 1.0/(1*domain.getVoxelSize()[1]) );

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

	//M_1 ---
	// is constant

	//M_2 ---
	// all components vanish

	//M_3 ---
	Eigen::Matrix<std::complex<double>, 3, 3> M_3;
	M_3(0, 0) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(0, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(1, 1) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.0, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(1, 0, domain.voxelToWorld(vd+V2d(0.0, 0.5)))));
	M_3(2, 2) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.0)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(1, 0, domain.voxelToWorld(vd+V2d(0.5, 0.0)))));
	Eigen::Matrix<double, 3, 3> M_3_real = (S*M_3*SInv).real();

	//b ---
	Eigen::Matrix<std::complex<double>, 3, 1> b;
	b(0, 0) = fields.q->eval(0, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(1, 0) = fields.q->eval(1, -1, domain.voxelToWorld(vd+V2d(0.0, 0.5)));
	b(2, 0) = fields.q->eval(1, 1, domain.voxelToWorld(vd+V2d(0.5, 0.0)));
	Eigen::Matrix<double, 3, 1> b_real = (S*b).real();

	// Assembling global system =============
	ctx.coeff_A( 1, vi + V2i(-1,0), 0 ) += -(h_inv[0]*0.57735026919);
	ctx.coeff_A( 1, vi + V2i(0,0), 0 ) += (h_inv[0]*0.57735026919);
	ctx.coeff_A( 0, vi + V2i(0,0), 1 ) += -(h_inv[0]*0.57735026919);
	ctx.coeff_A( 0, vi + V2i(1,0), 1 ) += (h_inv[0]*0.57735026919);
	ctx.coeff_A( 2, vi + V2i(0,-1), 0 ) += -(h_inv[1]*0.57735026919);
	ctx.coeff_A( 2, vi + V2i(0,0), 0 ) += (h_inv[1]*0.57735026919);
	ctx.coeff_A( 0, vi + V2i(0,0), 2 ) += -(h_inv[1]*0.57735026919);
	ctx.coeff_A( 0, vi + V2i(0,1), 2 ) += (h_inv[1]*0.57735026919);
	ctx.coeff_A( 0, vi + V2i(0,0), 0 ) += M_3_real.coeffRef(0, 0);
	ctx.coeff_A( 1, vi + V2i(0,0), 1 ) += M_3_real.coeffRef(1, 1);
	ctx.coeff_A( 2, vi + V2i(0,0), 2 ) += M_3_real.coeffRef(2, 2);
	ctx.coeff_b( 0 ) += b_real.coeffRef(0, 0);
	ctx.coeff_b( 1 ) += b_real.coeffRef(1, 0);
	ctx.coeff_b( 2 ) += b_real.coeffRef(2, 0);
}
V2i stencil_fopn_p1_sg_get_offset(int coeff)
{
	switch(coeff)
	{
		case 0:return V2i(1, 1);break;
		case 1:return V2i(0, 1);break;
		case 2:return V2i(1, 0);break;
		default:throw std::runtime_error("unexpected coefficient index");break;
	};
}
REGISTER_STENCIL(stencil_fopn_p1_sg, 1, 1)
