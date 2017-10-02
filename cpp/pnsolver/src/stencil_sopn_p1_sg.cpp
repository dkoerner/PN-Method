// This file was generated by stencil.py

#include <PNSystem.h>

void stencil_sopn_p1_sg(PNSystem::Stencil::Context& ctx)
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
	//M_0dxdxL + M_1dxdyL + M_2dxdzL + M_3dydxL + M_4dydyL + M_5dydzL + M_6dzdxL + M_7dzdyL + M_8dzdzL + M_9L + M_10dxL + M_11dyL + M_12dzL = b

	//M_0 ---
	// is constant

	//M_1 ---
	// is constant

	//M_2 ---
	// all components vanish

	//M_3 ---
	// is constant

	//M_4 ---
	// is constant

	//M_5 ---
	// all components vanish

	//M_6 ---
	// all components vanish

	//M_7 ---
	// all components vanish

	//M_8 ---
	// is constant

	//M_9 ---
	Eigen::Matrix<std::complex<double>, 3, 3> M_9;
	M_9(0, 0) = (std::pow(fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5))), 2)+
			-(fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(0, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_9(0, 1) = (-(0.408248290464*((-h_inv[0]*fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(h_inv[0]*fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(1.0, 0.5))))))+
			(std::complex<double>(0.0, 0.408248290463863)*((-h_inv[1]*fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(h_inv[1]*fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 1.0))))))+
			(0.408248290464*((-h_inv[0]*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(h_inv[0]*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(1.0, 0.5)))))*fields.f_p->eval(1, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5))))+
			(0.408248290464*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*((-h_inv[0]*fields.f_p->eval(1, 0, domain.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(h_inv[0]*fields.f_p->eval(1, 0, domain.voxelToWorld(vd+V2d(1.0, 0.5))))))+
			-(std::complex<double>(0.0, 0.408248290463863)*((-h_inv[1]*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(h_inv[1]*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 1.0)))))*fields.f_p->eval(1, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5))))+
			-(std::complex<double>(0.0, 0.408248290463863)*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*((-h_inv[1]*fields.f_p->eval(1, 0, domain.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(h_inv[1]*fields.f_p->eval(1, 0, domain.voxelToWorld(vd+V2d(0.5, 1.0)))))));
	M_9(0, 2) = ((0.408248290464*((-h_inv[0]*fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(h_inv[0]*fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(1.0, 0.5))))))+
			(std::complex<double>(0.0, 0.408248290463863)*((-h_inv[1]*fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(h_inv[1]*fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 1.0))))))+
			-(0.408248290464*((-h_inv[0]*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(h_inv[0]*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(1.0, 0.5)))))*fields.f_p->eval(1, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5))))+
			-(0.408248290464*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*((-h_inv[0]*fields.f_p->eval(1, 0, domain.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(h_inv[0]*fields.f_p->eval(1, 0, domain.voxelToWorld(vd+V2d(1.0, 0.5))))))+
			-(std::complex<double>(0.0, 0.408248290463863)*((-h_inv[1]*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(h_inv[1]*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 1.0)))))*fields.f_p->eval(1, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5))))+
			-(std::complex<double>(0.0, 0.408248290463863)*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*((-h_inv[1]*fields.f_p->eval(1, 0, domain.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(h_inv[1]*fields.f_p->eval(1, 0, domain.voxelToWorld(vd+V2d(0.5, 1.0)))))));
	M_9(1, 0) = (-(0.408248290464*((-h_inv[0]*fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(-0.5, 0.5))))+
			(h_inv[0]*fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5))))))+
			-(std::complex<double>(0.0, 0.408248290463863)*((-h_inv[1]*fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(h_inv[1]*fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.0, 1.0))))))+
			(0.408248290464*((-h_inv[0]*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(-0.5, 0.5))))+
			(h_inv[0]*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))))*fields.f_p->eval(0, 0, domain.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(0.408248290464*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.0, 0.5)))*((-h_inv[0]*fields.f_p->eval(0, 0, domain.voxelToWorld(vd+V2d(-0.5, 0.5))))+
			(h_inv[0]*fields.f_p->eval(0, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5))))))+
			(std::complex<double>(0.0, 0.408248290463863)*((-h_inv[1]*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(h_inv[1]*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.0, 1.0)))))*fields.f_p->eval(0, 0, domain.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(std::complex<double>(0.0, 0.408248290463863)*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.0, 0.5)))*((-h_inv[1]*fields.f_p->eval(0, 0, domain.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(h_inv[1]*fields.f_p->eval(0, 0, domain.voxelToWorld(vd+V2d(0.0, 1.0)))))));
	M_9(1, 1) = (std::pow(fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.0, 0.5))), 2)+
			-(fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(1, 0, domain.voxelToWorld(vd+V2d(0.0, 0.5)))));
	M_9(2, 0) = ((0.408248290464*((-h_inv[0]*fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(h_inv[0]*fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(1.0, 0.0))))))+
			-(std::complex<double>(0.0, 0.408248290463863)*((-h_inv[1]*fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, -0.5))))+
			(h_inv[1]*fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5))))))+
			-(0.408248290464*((-h_inv[0]*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(h_inv[0]*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(1.0, 0.0)))))*fields.f_p->eval(0, 0, domain.voxelToWorld(vd+V2d(0.5, 0.0))))+
			-(0.408248290464*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.0)))*((-h_inv[0]*fields.f_p->eval(0, 0, domain.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(h_inv[0]*fields.f_p->eval(0, 0, domain.voxelToWorld(vd+V2d(1.0, 0.0))))))+
			(std::complex<double>(0.0, 0.408248290463863)*((-h_inv[1]*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, -0.5))))+
			(h_inv[1]*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))))*fields.f_p->eval(0, 0, domain.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(std::complex<double>(0.0, 0.408248290463863)*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.0)))*((-h_inv[1]*fields.f_p->eval(0, 0, domain.voxelToWorld(vd+V2d(0.5, -0.5))))+
			(h_inv[1]*fields.f_p->eval(0, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))))));
	M_9(2, 2) = (std::pow(fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.0))), 2)+
			-(fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(1, 0, domain.voxelToWorld(vd+V2d(0.5, 0.0)))));
	Eigen::Matrix<double, 3, 3> M_9_real = (S*M_9*SInv).real();

	//M_10 ---
	Eigen::Matrix<std::complex<double>, 3, 3> M_10;
	M_10(0, 1) = (0.408248290464*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5))));
	M_10(0, 2) = -(0.408248290464*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5))));
	M_10(1, 0) = (0.408248290464*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(0, 0, domain.voxelToWorld(vd+V2d(0.0, 0.5))));
	M_10(2, 0) = -(0.408248290464*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(0, 0, domain.voxelToWorld(vd+V2d(0.5, 0.0))));
	Eigen::Matrix<double, 3, 3> M_10_real = (S*M_10*SInv).real();

	//M_11 ---
	Eigen::Matrix<std::complex<double>, 3, 3> M_11;
	M_11(0, 1) = -(std::complex<double>(0.0, 0.408248290463863)*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5))));
	M_11(0, 2) = -(std::complex<double>(0.0, 0.408248290463863)*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5))));
	M_11(1, 0) = (std::complex<double>(0.0, 0.408248290463863)*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(0, 0, domain.voxelToWorld(vd+V2d(0.0, 0.5))));
	M_11(2, 0) = (std::complex<double>(0.0, 0.408248290463863)*fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(0, 0, domain.voxelToWorld(vd+V2d(0.5, 0.0))));
	Eigen::Matrix<double, 3, 3> M_11_real = (S*M_11*SInv).real();

	//M_12 ---
	// all components vanish

	//b ---
	Eigen::Matrix<std::complex<double>, 3, 1> b;
	b(0, 0) = (-(0.408248290464*((-h_inv[0]*fields.q->eval(1, -1, domain.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(h_inv[0]*fields.q->eval(1, -1, domain.voxelToWorld(vd+V2d(1.0, 0.5))))))+
			(0.408248290464*((-h_inv[0]*fields.q->eval(1, 1, domain.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(h_inv[0]*fields.q->eval(1, 1, domain.voxelToWorld(vd+V2d(1.0, 0.5))))))+
			(std::complex<double>(0.0, 0.408248290463863)*((-h_inv[1]*fields.q->eval(1, -1, domain.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(h_inv[1]*fields.q->eval(1, -1, domain.voxelToWorld(vd+V2d(0.5, 1.0))))))+
			(std::complex<double>(0.0, 0.408248290463863)*((-h_inv[1]*fields.q->eval(1, 1, domain.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(h_inv[1]*fields.q->eval(1, 1, domain.voxelToWorld(vd+V2d(0.5, 1.0))))))+
			(fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.q->eval(0, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	b(1, 0) = (-(0.408248290464*((-h_inv[0]*fields.q->eval(0, 0, domain.voxelToWorld(vd+V2d(-0.5, 0.5))))+
			(h_inv[0]*fields.q->eval(0, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5))))))+
			-(std::complex<double>(0.0, 0.408248290463863)*((-h_inv[1]*fields.q->eval(0, 0, domain.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(h_inv[1]*fields.q->eval(0, 0, domain.voxelToWorld(vd+V2d(0.0, 1.0))))))+
			(fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.q->eval(1, -1, domain.voxelToWorld(vd+V2d(0.0, 0.5)))));
	b(2, 0) = ((0.408248290464*((-h_inv[0]*fields.q->eval(0, 0, domain.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(h_inv[0]*fields.q->eval(0, 0, domain.voxelToWorld(vd+V2d(1.0, 0.0))))))+
			-(std::complex<double>(0.0, 0.408248290463863)*((-h_inv[1]*fields.q->eval(0, 0, domain.voxelToWorld(vd+V2d(0.5, -0.5))))+
			(h_inv[1]*fields.q->eval(0, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5))))))+
			(fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.q->eval(1, 1, domain.voxelToWorld(vd+V2d(0.5, 0.0)))));
	Eigen::Matrix<double, 3, 1> b_real = (S*b).real();

	// Assembling global system =============
	ctx.coeff_A( 0, vi + V2i(-1,0), 0 ) += (h_inv[0]*h_inv[0]*-0.333333333333);
	ctx.coeff_A( 0, vi + V2i(0,0), 0 ) += -(h_inv[0]*h_inv[0]*-0.333333333333);
	ctx.coeff_A( 0, vi + V2i(0,0), 0 ) += -(h_inv[0]*h_inv[0]*-0.333333333333);
	ctx.coeff_A( 0, vi + V2i(1,0), 0 ) += (h_inv[0]*h_inv[0]*-0.333333333333);
	ctx.coeff_A( 1, vi + V2i(-1,0), 1 ) += (h_inv[0]*h_inv[0]*-0.6);
	ctx.coeff_A( 1, vi + V2i(0,0), 1 ) += -(h_inv[0]*h_inv[0]*-0.6);
	ctx.coeff_A( 1, vi + V2i(0,0), 1 ) += -(h_inv[0]*h_inv[0]*-0.6);
	ctx.coeff_A( 1, vi + V2i(1,0), 1 ) += (h_inv[0]*h_inv[0]*-0.6);
	ctx.coeff_A( 2, vi + V2i(-1,0), 2 ) += (h_inv[0]*h_inv[0]*-0.2);
	ctx.coeff_A( 2, vi + V2i(0,0), 2 ) += -(h_inv[0]*h_inv[0]*-0.2);
	ctx.coeff_A( 2, vi + V2i(0,0), 2 ) += -(h_inv[0]*h_inv[0]*-0.2);
	ctx.coeff_A( 2, vi + V2i(1,0), 2 ) += (h_inv[0]*h_inv[0]*-0.2);
	ctx.coeff_A( 2, vi + V2i(0,-1), 1 ) += (h_inv[1]*h_inv[0]*-0.2);
	ctx.coeff_A( 2, vi + V2i(0,0), 1 ) += -(h_inv[1]*h_inv[0]*-0.2);
	ctx.coeff_A( 2, vi + V2i(1,-1), 1 ) += -(h_inv[1]*h_inv[0]*-0.2);
	ctx.coeff_A( 2, vi + V2i(1,0), 1 ) += (h_inv[1]*h_inv[0]*-0.2);
	ctx.coeff_A( 1, vi + V2i(-1,0), 2 ) += (h_inv[1]*h_inv[0]*-0.2);
	ctx.coeff_A( 1, vi + V2i(-1,1), 2 ) += -(h_inv[1]*h_inv[0]*-0.2);
	ctx.coeff_A( 1, vi + V2i(0,0), 2 ) += -(h_inv[1]*h_inv[0]*-0.2);
	ctx.coeff_A( 1, vi + V2i(0,1), 2 ) += (h_inv[1]*h_inv[0]*-0.2);
	ctx.coeff_A( 2, vi + V2i(0,-1), 1 ) += (h_inv[0]*h_inv[1]*-0.2);
	ctx.coeff_A( 2, vi + V2i(1,-1), 1 ) += -(h_inv[0]*h_inv[1]*-0.2);
	ctx.coeff_A( 2, vi + V2i(0,0), 1 ) += -(h_inv[0]*h_inv[1]*-0.2);
	ctx.coeff_A( 2, vi + V2i(1,0), 1 ) += (h_inv[0]*h_inv[1]*-0.2);
	ctx.coeff_A( 1, vi + V2i(-1,0), 2 ) += (h_inv[0]*h_inv[1]*-0.2);
	ctx.coeff_A( 1, vi + V2i(0,0), 2 ) += -(h_inv[0]*h_inv[1]*-0.2);
	ctx.coeff_A( 1, vi + V2i(-1,1), 2 ) += -(h_inv[0]*h_inv[1]*-0.2);
	ctx.coeff_A( 1, vi + V2i(0,1), 2 ) += (h_inv[0]*h_inv[1]*-0.2);
	ctx.coeff_A( 0, vi + V2i(0,-1), 0 ) += (h_inv[1]*h_inv[1]*-0.333333333333);
	ctx.coeff_A( 0, vi + V2i(0,0), 0 ) += -(h_inv[1]*h_inv[1]*-0.333333333333);
	ctx.coeff_A( 0, vi + V2i(0,0), 0 ) += -(h_inv[1]*h_inv[1]*-0.333333333333);
	ctx.coeff_A( 0, vi + V2i(0,1), 0 ) += (h_inv[1]*h_inv[1]*-0.333333333333);
	ctx.coeff_A( 1, vi + V2i(0,-1), 1 ) += (h_inv[1]*h_inv[1]*-0.2);
	ctx.coeff_A( 1, vi + V2i(0,0), 1 ) += -(h_inv[1]*h_inv[1]*-0.2);
	ctx.coeff_A( 1, vi + V2i(0,0), 1 ) += -(h_inv[1]*h_inv[1]*-0.2);
	ctx.coeff_A( 1, vi + V2i(0,1), 1 ) += (h_inv[1]*h_inv[1]*-0.2);
	ctx.coeff_A( 2, vi + V2i(0,-1), 2 ) += (h_inv[1]*h_inv[1]*-0.6);
	ctx.coeff_A( 2, vi + V2i(0,0), 2 ) += -(h_inv[1]*h_inv[1]*-0.6);
	ctx.coeff_A( 2, vi + V2i(0,0), 2 ) += -(h_inv[1]*h_inv[1]*-0.6);
	ctx.coeff_A( 2, vi + V2i(0,1), 2 ) += (h_inv[1]*h_inv[1]*-0.6);
	ctx.coeff_A( 0, vi + V2i(0,0), 0 ) += M_9_real.coeffRef(0, 0);
	ctx.coeff_A( 1, vi + V2i(0,0), 0 ) += (0.5*M_9_real.coeffRef(1, 0));
	ctx.coeff_A( 1, vi + V2i(-1,0), 0 ) += (0.5*M_9_real.coeffRef(1, 0));
	ctx.coeff_A( 2, vi + V2i(0,0), 0 ) += (0.5*M_9_real.coeffRef(2, 0));
	ctx.coeff_A( 2, vi + V2i(0,-1), 0 ) += (0.5*M_9_real.coeffRef(2, 0));
	ctx.coeff_A( 0, vi + V2i(1,0), 1 ) += (0.5*M_9_real.coeffRef(0, 1));
	ctx.coeff_A( 0, vi + V2i(0,0), 1 ) += (0.5*M_9_real.coeffRef(0, 1));
	ctx.coeff_A( 1, vi + V2i(0,0), 1 ) += M_9_real.coeffRef(1, 1);
	ctx.coeff_A( 0, vi + V2i(0,1), 2 ) += (0.5*M_9_real.coeffRef(0, 2));
	ctx.coeff_A( 0, vi + V2i(0,0), 2 ) += (0.5*M_9_real.coeffRef(0, 2));
	ctx.coeff_A( 2, vi + V2i(0,0), 2 ) += M_9_real.coeffRef(2, 2);
	ctx.coeff_A( 1, vi + V2i(-1,0), 0 ) += -(h_inv[0]*M_10_real.coeffRef(1, 0));
	ctx.coeff_A( 1, vi + V2i(0,0), 0 ) += (h_inv[0]*M_10_real.coeffRef(1, 0));
	ctx.coeff_A( 2, vi + V2i(-1,-1), 0 ) += -(0.25*h_inv[0]*M_10_real.coeffRef(2, 0));
	ctx.coeff_A( 2, vi + V2i(-1,0), 0 ) += -(0.25*h_inv[0]*M_10_real.coeffRef(2, 0));
	ctx.coeff_A( 2, vi + V2i(0,-1), 0 ) += -(0.25*h_inv[0]*M_10_real.coeffRef(2, 0));
	ctx.coeff_A( 2, vi + V2i(0,0), 0 ) += -(0.25*h_inv[0]*M_10_real.coeffRef(2, 0));
	ctx.coeff_A( 2, vi + V2i(0,-1), 0 ) += (0.25*h_inv[0]*M_10_real.coeffRef(2, 0));
	ctx.coeff_A( 2, vi + V2i(0,0), 0 ) += (0.25*h_inv[0]*M_10_real.coeffRef(2, 0));
	ctx.coeff_A( 2, vi + V2i(1,-1), 0 ) += (0.25*h_inv[0]*M_10_real.coeffRef(2, 0));
	ctx.coeff_A( 2, vi + V2i(1,0), 0 ) += (0.25*h_inv[0]*M_10_real.coeffRef(2, 0));
	ctx.coeff_A( 0, vi + V2i(0,0), 1 ) += -(h_inv[0]*M_10_real.coeffRef(0, 1));
	ctx.coeff_A( 0, vi + V2i(1,0), 1 ) += (h_inv[0]*M_10_real.coeffRef(0, 1));
	ctx.coeff_A( 0, vi + V2i(-1,0), 2 ) += -(0.25*h_inv[0]*M_10_real.coeffRef(0, 2));
	ctx.coeff_A( 0, vi + V2i(-1,1), 2 ) += -(0.25*h_inv[0]*M_10_real.coeffRef(0, 2));
	ctx.coeff_A( 0, vi + V2i(0,0), 2 ) += -(0.25*h_inv[0]*M_10_real.coeffRef(0, 2));
	ctx.coeff_A( 0, vi + V2i(0,1), 2 ) += -(0.25*h_inv[0]*M_10_real.coeffRef(0, 2));
	ctx.coeff_A( 0, vi + V2i(0,0), 2 ) += (0.25*h_inv[0]*M_10_real.coeffRef(0, 2));
	ctx.coeff_A( 0, vi + V2i(0,1), 2 ) += (0.25*h_inv[0]*M_10_real.coeffRef(0, 2));
	ctx.coeff_A( 0, vi + V2i(1,0), 2 ) += (0.25*h_inv[0]*M_10_real.coeffRef(0, 2));
	ctx.coeff_A( 0, vi + V2i(1,1), 2 ) += (0.25*h_inv[0]*M_10_real.coeffRef(0, 2));
	ctx.coeff_A( 1, vi + V2i(-1,-1), 0 ) += -(0.25*h_inv[1]*M_11_real.coeffRef(1, 0));
	ctx.coeff_A( 1, vi + V2i(-1,0), 0 ) += -(0.25*h_inv[1]*M_11_real.coeffRef(1, 0));
	ctx.coeff_A( 1, vi + V2i(0,-1), 0 ) += -(0.25*h_inv[1]*M_11_real.coeffRef(1, 0));
	ctx.coeff_A( 1, vi + V2i(0,0), 0 ) += -(0.25*h_inv[1]*M_11_real.coeffRef(1, 0));
	ctx.coeff_A( 1, vi + V2i(-1,0), 0 ) += (0.25*h_inv[1]*M_11_real.coeffRef(1, 0));
	ctx.coeff_A( 1, vi + V2i(-1,1), 0 ) += (0.25*h_inv[1]*M_11_real.coeffRef(1, 0));
	ctx.coeff_A( 1, vi + V2i(0,0), 0 ) += (0.25*h_inv[1]*M_11_real.coeffRef(1, 0));
	ctx.coeff_A( 1, vi + V2i(0,1), 0 ) += (0.25*h_inv[1]*M_11_real.coeffRef(1, 0));
	ctx.coeff_A( 2, vi + V2i(0,-1), 0 ) += -(h_inv[1]*M_11_real.coeffRef(2, 0));
	ctx.coeff_A( 2, vi + V2i(0,0), 0 ) += (h_inv[1]*M_11_real.coeffRef(2, 0));
	ctx.coeff_A( 0, vi + V2i(0,-1), 1 ) += -(0.25*h_inv[1]*M_11_real.coeffRef(0, 1));
	ctx.coeff_A( 0, vi + V2i(0,0), 1 ) += -(0.25*h_inv[1]*M_11_real.coeffRef(0, 1));
	ctx.coeff_A( 0, vi + V2i(1,-1), 1 ) += -(0.25*h_inv[1]*M_11_real.coeffRef(0, 1));
	ctx.coeff_A( 0, vi + V2i(1,0), 1 ) += -(0.25*h_inv[1]*M_11_real.coeffRef(0, 1));
	ctx.coeff_A( 0, vi + V2i(0,0), 1 ) += (0.25*h_inv[1]*M_11_real.coeffRef(0, 1));
	ctx.coeff_A( 0, vi + V2i(0,1), 1 ) += (0.25*h_inv[1]*M_11_real.coeffRef(0, 1));
	ctx.coeff_A( 0, vi + V2i(1,0), 1 ) += (0.25*h_inv[1]*M_11_real.coeffRef(0, 1));
	ctx.coeff_A( 0, vi + V2i(1,1), 1 ) += (0.25*h_inv[1]*M_11_real.coeffRef(0, 1));
	ctx.coeff_A( 0, vi + V2i(0,0), 2 ) += -(h_inv[1]*M_11_real.coeffRef(0, 2));
	ctx.coeff_A( 0, vi + V2i(0,1), 2 ) += (h_inv[1]*M_11_real.coeffRef(0, 2));
	ctx.coeff_b( 0 ) += b_real.coeffRef(0, 0);
	ctx.coeff_b( 1 ) += b_real.coeffRef(1, 0);
	ctx.coeff_b( 2 ) += b_real.coeffRef(2, 0);
}
V2i stencil_sopn_p1_sg_get_offset(int coeff)
{
	switch(coeff)
	{
		case 0:return V2i(1, 1);break;
		case 1:return V2i(0, 1);break;
		case 2:return V2i(1, 0);break;
		default:throw std::runtime_error("unexpected coefficient index");break;
	};
}
REGISTER_STENCIL(stencil_sopn_p1_sg, 1, 1)
