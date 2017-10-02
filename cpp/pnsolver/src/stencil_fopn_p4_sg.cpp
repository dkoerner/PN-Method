// This file was generated by stencil.py

#include <PNSystem.h>

void stencil_fopn_p4_sg(PNSystem::Stencil::Context& ctx)
{
	V2i vi = ctx.getVoxel();
	V2d vd = vi.cast<double>();
	const Domain& domain = ctx.getDomain();
	const PNSystem::Fields& fields = ctx.getFields();
	V2d h_inv( 1.0/(1*domain.getVoxelSize()[0]), 1.0/(1*domain.getVoxelSize()[1]) );

	Eigen::Matrix<std::complex<double>, 15, 15> S;
	S.coeffRef(0, 0) = std::complex<double>(1.0, 0.0);
	S.coeffRef(1, 1) = std::complex<double>(0.7071067811865475, 0.0);
	S.coeffRef(1, 2) = std::complex<double>(-0.7071067811865475, 0.0);
	S.coeffRef(2, 1) = std::complex<double>(-0.0, -0.7071067811865475);
	S.coeffRef(2, 2) = std::complex<double>(-0.0, -0.7071067811865475);
	S.coeffRef(3, 3) = std::complex<double>(0.7071067811865475, 0.0);
	S.coeffRef(3, 5) = std::complex<double>(0.7071067811865475, 0.0);
	S.coeffRef(4, 3) = std::complex<double>(-0.0, -0.7071067811865475);
	S.coeffRef(4, 5) = std::complex<double>(0.0, 0.7071067811865475);
	S.coeffRef(5, 4) = std::complex<double>(1.0, 0.0);
	S.coeffRef(6, 6) = std::complex<double>(0.7071067811865475, 0.0);
	S.coeffRef(6, 9) = std::complex<double>(-0.7071067811865475, 0.0);
	S.coeffRef(7, 6) = std::complex<double>(-0.0, -0.7071067811865475);
	S.coeffRef(7, 9) = std::complex<double>(-0.0, -0.7071067811865475);
	S.coeffRef(8, 7) = std::complex<double>(0.7071067811865475, 0.0);
	S.coeffRef(8, 8) = std::complex<double>(-0.7071067811865475, 0.0);
	S.coeffRef(9, 7) = std::complex<double>(-0.0, -0.7071067811865475);
	S.coeffRef(9, 8) = std::complex<double>(-0.0, -0.7071067811865475);
	S.coeffRef(10, 10) = std::complex<double>(0.7071067811865475, 0.0);
	S.coeffRef(10, 14) = std::complex<double>(0.7071067811865475, 0.0);
	S.coeffRef(11, 10) = std::complex<double>(-0.0, -0.7071067811865475);
	S.coeffRef(11, 14) = std::complex<double>(0.0, 0.7071067811865475);
	S.coeffRef(12, 11) = std::complex<double>(0.7071067811865475, 0.0);
	S.coeffRef(12, 13) = std::complex<double>(0.7071067811865475, 0.0);
	S.coeffRef(13, 11) = std::complex<double>(-0.0, -0.7071067811865475);
	S.coeffRef(13, 13) = std::complex<double>(0.0, 0.7071067811865475);
	S.coeffRef(14, 12) = std::complex<double>(1.0, 0.0);
	Eigen::Matrix<std::complex<double>, 15, 15> SInv;
	SInv.coeffRef(0, 0) = std::complex<double>(1.0, 0.0);
	SInv.coeffRef(1, 1) = std::complex<double>(0.7071067811865476, 0.0);
	SInv.coeffRef(1, 2) = std::complex<double>(0.0, 0.7071067811865476);
	SInv.coeffRef(2, 1) = std::complex<double>(-0.7071067811865476, 0.0);
	SInv.coeffRef(2, 2) = std::complex<double>(-0.0, 0.7071067811865476);
	SInv.coeffRef(3, 3) = std::complex<double>(0.7071067811865476, 0.0);
	SInv.coeffRef(3, 4) = std::complex<double>(0.0, 0.7071067811865476);
	SInv.coeffRef(4, 5) = std::complex<double>(1.0, 0.0);
	SInv.coeffRef(5, 3) = std::complex<double>(0.7071067811865476, 0.0);
	SInv.coeffRef(5, 4) = std::complex<double>(0.0, -0.7071067811865476);
	SInv.coeffRef(6, 6) = std::complex<double>(0.7071067811865476, 0.0);
	SInv.coeffRef(6, 7) = std::complex<double>(0.0, 0.7071067811865476);
	SInv.coeffRef(7, 8) = std::complex<double>(0.7071067811865476, 0.0);
	SInv.coeffRef(7, 9) = std::complex<double>(0.0, 0.7071067811865476);
	SInv.coeffRef(8, 8) = std::complex<double>(-0.7071067811865476, 0.0);
	SInv.coeffRef(8, 9) = std::complex<double>(-0.0, 0.7071067811865476);
	SInv.coeffRef(9, 6) = std::complex<double>(-0.7071067811865476, 0.0);
	SInv.coeffRef(9, 7) = std::complex<double>(-0.0, 0.7071067811865476);
	SInv.coeffRef(10, 10) = std::complex<double>(0.7071067811865476, 0.0);
	SInv.coeffRef(10, 11) = std::complex<double>(0.0, 0.7071067811865476);
	SInv.coeffRef(11, 12) = std::complex<double>(0.7071067811865476, 0.0);
	SInv.coeffRef(11, 13) = std::complex<double>(0.0, 0.7071067811865476);
	SInv.coeffRef(12, 14) = std::complex<double>(1.0, 0.0);
	SInv.coeffRef(13, 12) = std::complex<double>(0.7071067811865476, 0.0);
	SInv.coeffRef(13, 13) = std::complex<double>(0.0, -0.7071067811865476);
	SInv.coeffRef(14, 10) = std::complex<double>(0.7071067811865476, 0.0);
	SInv.coeffRef(14, 11) = std::complex<double>(0.0, -0.7071067811865476);

	//Producing complex-valued matrices =============
	//M_0dxL + M_1dyL + M_2dzL + M_3L = b

	//M_0 ---
	// is constant

	//M_1 ---
	// is constant

	//M_2 ---
	// all components vanish

	//M_3 ---
	Eigen::Matrix<std::complex<double>, 15, 15> M_3;
	M_3(0, 0) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(0, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(1, 1) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.0, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(1, 0, domain.voxelToWorld(vd+V2d(0.0, 0.5)))));
	M_3(2, 2) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.0)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(1, 0, domain.voxelToWorld(vd+V2d(0.5, 0.0)))));
	M_3(3, 3) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(2, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(4, 4) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.0, 0.0)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.0, 0.0)))*fields.f_p->eval(2, 0, domain.voxelToWorld(vd+V2d(0.0, 0.0)))));
	M_3(5, 5) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(2, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(6, 6) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.0, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(3, 0, domain.voxelToWorld(vd+V2d(0.0, 0.5)))));
	M_3(7, 7) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.0)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(3, 0, domain.voxelToWorld(vd+V2d(0.5, 0.0)))));
	M_3(8, 8) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.0, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(3, 0, domain.voxelToWorld(vd+V2d(0.0, 0.5)))));
	M_3(9, 9) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.0)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(3, 0, domain.voxelToWorld(vd+V2d(0.5, 0.0)))));
	M_3(10, 10) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(4, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(11, 11) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.0, 0.0)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.0, 0.0)))*fields.f_p->eval(4, 0, domain.voxelToWorld(vd+V2d(0.0, 0.0)))));
	M_3(12, 12) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(4, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(13, 13) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.0, 0.0)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.0, 0.0)))*fields.f_p->eval(4, 0, domain.voxelToWorld(vd+V2d(0.0, 0.0)))));
	M_3(14, 14) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(4, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	Eigen::Matrix<double, 15, 15> M_3_real = (S*M_3*SInv).real();

	//b ---
	Eigen::Matrix<std::complex<double>, 15, 1> b;
	b(0, 0) = fields.q->eval(0, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(1, 0) = fields.q->eval(1, -1, domain.voxelToWorld(vd+V2d(0.0, 0.5)));
	b(2, 0) = fields.q->eval(1, 1, domain.voxelToWorld(vd+V2d(0.5, 0.0)));
	b(3, 0) = fields.q->eval(2, -2, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(4, 0) = fields.q->eval(2, 0, domain.voxelToWorld(vd+V2d(0.0, 0.0)));
	b(5, 0) = fields.q->eval(2, 2, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(6, 0) = fields.q->eval(3, -3, domain.voxelToWorld(vd+V2d(0.0, 0.5)));
	b(7, 0) = fields.q->eval(3, -1, domain.voxelToWorld(vd+V2d(0.5, 0.0)));
	b(8, 0) = fields.q->eval(3, 1, domain.voxelToWorld(vd+V2d(0.0, 0.5)));
	b(9, 0) = fields.q->eval(3, 3, domain.voxelToWorld(vd+V2d(0.5, 0.0)));
	b(10, 0) = fields.q->eval(4, -4, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(11, 0) = fields.q->eval(4, -2, domain.voxelToWorld(vd+V2d(0.0, 0.0)));
	b(12, 0) = fields.q->eval(4, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(13, 0) = fields.q->eval(4, 2, domain.voxelToWorld(vd+V2d(0.0, 0.0)));
	b(14, 0) = fields.q->eval(4, 4, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	Eigen::Matrix<double, 15, 1> b_real = (S*b).real();

	// Assembling global system =============
	ctx.coeff_A( 1, vi + V2i(-1,0), 0 ) += -(h_inv[0]*0.57735026919);
	ctx.coeff_A( 1, vi + V2i(0,0), 0 ) += (h_inv[0]*0.57735026919);
	ctx.coeff_A( 0, vi + V2i(0,0), 1 ) += -(h_inv[0]*0.57735026919);
	ctx.coeff_A( 0, vi + V2i(1,0), 1 ) += (h_inv[0]*0.57735026919);
	ctx.coeff_A( 3, vi + V2i(0,0), 1 ) += -(h_inv[0]*0.4472135955);
	ctx.coeff_A( 3, vi + V2i(1,0), 1 ) += (h_inv[0]*0.4472135955);
	ctx.coeff_A( 5, vi + V2i(0,0), 1 ) += -(h_inv[0]*-0.258198889747);
	ctx.coeff_A( 5, vi + V2i(1,0), 1 ) += (h_inv[0]*-0.258198889747);
	ctx.coeff_A( 4, vi + V2i(-1,0), 2 ) += -(h_inv[0]*0.4472135955);
	ctx.coeff_A( 4, vi + V2i(0,0), 2 ) += (h_inv[0]*0.4472135955);
	ctx.coeff_A( 1, vi + V2i(-1,0), 3 ) += -(h_inv[0]*0.4472135955);
	ctx.coeff_A( 1, vi + V2i(0,0), 3 ) += (h_inv[0]*0.4472135955);
	ctx.coeff_A( 6, vi + V2i(-1,0), 3 ) += -(h_inv[0]*0.462910049886);
	ctx.coeff_A( 6, vi + V2i(0,0), 3 ) += (h_inv[0]*0.462910049886);
	ctx.coeff_A( 8, vi + V2i(-1,0), 3 ) += -(h_inv[0]*-0.119522860933);
	ctx.coeff_A( 8, vi + V2i(0,0), 3 ) += (h_inv[0]*-0.119522860933);
	ctx.coeff_A( 2, vi + V2i(0,0), 4 ) += -(h_inv[0]*0.4472135955);
	ctx.coeff_A( 2, vi + V2i(1,0), 4 ) += (h_inv[0]*0.4472135955);
	ctx.coeff_A( 7, vi + V2i(0,0), 4 ) += -(h_inv[0]*0.462910049886);
	ctx.coeff_A( 7, vi + V2i(1,0), 4 ) += (h_inv[0]*0.462910049886);
	ctx.coeff_A( 9, vi + V2i(0,0), 4 ) += -(h_inv[0]*-0.119522860933);
	ctx.coeff_A( 9, vi + V2i(1,0), 4 ) += (h_inv[0]*-0.119522860933);
	ctx.coeff_A( 1, vi + V2i(-1,0), 5 ) += -(h_inv[0]*-0.258198889747);
	ctx.coeff_A( 1, vi + V2i(0,0), 5 ) += (h_inv[0]*-0.258198889747);
	ctx.coeff_A( 8, vi + V2i(-1,0), 5 ) += -(h_inv[0]*0.414039335605);
	ctx.coeff_A( 8, vi + V2i(0,0), 5 ) += (h_inv[0]*0.414039335605);
	ctx.coeff_A( 3, vi + V2i(0,0), 6 ) += -(h_inv[0]*0.462910049886);
	ctx.coeff_A( 3, vi + V2i(1,0), 6 ) += (h_inv[0]*0.462910049886);
	ctx.coeff_A( 10, vi + V2i(0,0), 6 ) += -(h_inv[0]*0.471404520791);
	ctx.coeff_A( 10, vi + V2i(1,0), 6 ) += (h_inv[0]*0.471404520791);
	ctx.coeff_A( 12, vi + V2i(0,0), 6 ) += -(h_inv[0]*-0.0890870806375);
	ctx.coeff_A( 12, vi + V2i(1,0), 6 ) += (h_inv[0]*-0.0890870806375);
	ctx.coeff_A( 4, vi + V2i(-1,0), 7 ) += -(h_inv[0]*0.462910049886);
	ctx.coeff_A( 4, vi + V2i(0,0), 7 ) += (h_inv[0]*0.462910049886);
	ctx.coeff_A( 11, vi + V2i(-1,0), 7 ) += -(h_inv[0]*0.471404520791);
	ctx.coeff_A( 11, vi + V2i(0,0), 7 ) += (h_inv[0]*0.471404520791);
	ctx.coeff_A( 13, vi + V2i(-1,0), 7 ) += -(h_inv[0]*-0.0890870806375);
	ctx.coeff_A( 13, vi + V2i(0,0), 7 ) += (h_inv[0]*-0.0890870806375);
	ctx.coeff_A( 3, vi + V2i(0,0), 8 ) += -(h_inv[0]*-0.119522860933);
	ctx.coeff_A( 3, vi + V2i(1,0), 8 ) += (h_inv[0]*-0.119522860933);
	ctx.coeff_A( 5, vi + V2i(0,0), 8 ) += -(h_inv[0]*0.414039335605);
	ctx.coeff_A( 5, vi + V2i(1,0), 8 ) += (h_inv[0]*0.414039335605);
	ctx.coeff_A( 12, vi + V2i(0,0), 8 ) += -(h_inv[0]*0.345032779671);
	ctx.coeff_A( 12, vi + V2i(1,0), 8 ) += (h_inv[0]*0.345032779671);
	ctx.coeff_A( 14, vi + V2i(0,0), 8 ) += -(h_inv[0]*-0.308606699924);
	ctx.coeff_A( 14, vi + V2i(1,0), 8 ) += (h_inv[0]*-0.308606699924);
	ctx.coeff_A( 4, vi + V2i(-1,0), 9 ) += -(h_inv[0]*-0.119522860933);
	ctx.coeff_A( 4, vi + V2i(0,0), 9 ) += (h_inv[0]*-0.119522860933);
	ctx.coeff_A( 13, vi + V2i(-1,0), 9 ) += -(h_inv[0]*0.345032779671);
	ctx.coeff_A( 13, vi + V2i(0,0), 9 ) += (h_inv[0]*0.345032779671);
	ctx.coeff_A( 6, vi + V2i(-1,0), 10 ) += -(h_inv[0]*0.471404520791);
	ctx.coeff_A( 6, vi + V2i(0,0), 10 ) += (h_inv[0]*0.471404520791);
	ctx.coeff_A( 7, vi + V2i(0,0), 11 ) += -(h_inv[0]*0.471404520791);
	ctx.coeff_A( 7, vi + V2i(1,0), 11 ) += (h_inv[0]*0.471404520791);
	ctx.coeff_A( 6, vi + V2i(-1,0), 12 ) += -(h_inv[0]*-0.0890870806375);
	ctx.coeff_A( 6, vi + V2i(0,0), 12 ) += (h_inv[0]*-0.0890870806375);
	ctx.coeff_A( 8, vi + V2i(-1,0), 12 ) += -(h_inv[0]*0.345032779671);
	ctx.coeff_A( 8, vi + V2i(0,0), 12 ) += (h_inv[0]*0.345032779671);
	ctx.coeff_A( 7, vi + V2i(0,0), 13 ) += -(h_inv[0]*-0.0890870806375);
	ctx.coeff_A( 7, vi + V2i(1,0), 13 ) += (h_inv[0]*-0.0890870806375);
	ctx.coeff_A( 9, vi + V2i(0,0), 13 ) += -(h_inv[0]*0.345032779671);
	ctx.coeff_A( 9, vi + V2i(1,0), 13 ) += (h_inv[0]*0.345032779671);
	ctx.coeff_A( 8, vi + V2i(-1,0), 14 ) += -(h_inv[0]*-0.308606699924);
	ctx.coeff_A( 8, vi + V2i(0,0), 14 ) += (h_inv[0]*-0.308606699924);
	ctx.coeff_A( 2, vi + V2i(0,-1), 0 ) += -(h_inv[1]*0.57735026919);
	ctx.coeff_A( 2, vi + V2i(0,0), 0 ) += (h_inv[1]*0.57735026919);
	ctx.coeff_A( 4, vi + V2i(0,-1), 1 ) += -(h_inv[1]*0.4472135955);
	ctx.coeff_A( 4, vi + V2i(0,0), 1 ) += (h_inv[1]*0.4472135955);
	ctx.coeff_A( 0, vi + V2i(0,0), 2 ) += -(h_inv[1]*0.57735026919);
	ctx.coeff_A( 0, vi + V2i(0,1), 2 ) += (h_inv[1]*0.57735026919);
	ctx.coeff_A( 3, vi + V2i(0,0), 2 ) += -(h_inv[1]*-0.4472135955);
	ctx.coeff_A( 3, vi + V2i(0,1), 2 ) += (h_inv[1]*-0.4472135955);
	ctx.coeff_A( 5, vi + V2i(0,0), 2 ) += -(h_inv[1]*-0.258198889747);
	ctx.coeff_A( 5, vi + V2i(0,1), 2 ) += (h_inv[1]*-0.258198889747);
	ctx.coeff_A( 2, vi + V2i(0,-1), 3 ) += -(h_inv[1]*-0.4472135955);
	ctx.coeff_A( 2, vi + V2i(0,0), 3 ) += (h_inv[1]*-0.4472135955);
	ctx.coeff_A( 7, vi + V2i(0,-1), 3 ) += -(h_inv[1]*0.462910049886);
	ctx.coeff_A( 7, vi + V2i(0,0), 3 ) += (h_inv[1]*0.462910049886);
	ctx.coeff_A( 9, vi + V2i(0,-1), 3 ) += -(h_inv[1]*0.119522860933);
	ctx.coeff_A( 9, vi + V2i(0,0), 3 ) += (h_inv[1]*0.119522860933);
	ctx.coeff_A( 1, vi + V2i(0,0), 4 ) += -(h_inv[1]*0.4472135955);
	ctx.coeff_A( 1, vi + V2i(0,1), 4 ) += (h_inv[1]*0.4472135955);
	ctx.coeff_A( 6, vi + V2i(0,0), 4 ) += -(h_inv[1]*-0.462910049886);
	ctx.coeff_A( 6, vi + V2i(0,1), 4 ) += (h_inv[1]*-0.462910049886);
	ctx.coeff_A( 8, vi + V2i(0,0), 4 ) += -(h_inv[1]*-0.119522860933);
	ctx.coeff_A( 8, vi + V2i(0,1), 4 ) += (h_inv[1]*-0.119522860933);
	ctx.coeff_A( 2, vi + V2i(0,-1), 5 ) += -(h_inv[1]*-0.258198889747);
	ctx.coeff_A( 2, vi + V2i(0,0), 5 ) += (h_inv[1]*-0.258198889747);
	ctx.coeff_A( 9, vi + V2i(0,-1), 5 ) += -(h_inv[1]*0.414039335605);
	ctx.coeff_A( 9, vi + V2i(0,0), 5 ) += (h_inv[1]*0.414039335605);
	ctx.coeff_A( 4, vi + V2i(0,-1), 6 ) += -(h_inv[1]*-0.462910049886);
	ctx.coeff_A( 4, vi + V2i(0,0), 6 ) += (h_inv[1]*-0.462910049886);
	ctx.coeff_A( 11, vi + V2i(0,-1), 6 ) += -(h_inv[1]*0.471404520791);
	ctx.coeff_A( 11, vi + V2i(0,0), 6 ) += (h_inv[1]*0.471404520791);
	ctx.coeff_A( 13, vi + V2i(0,-1), 6 ) += -(h_inv[1]*0.0890870806375);
	ctx.coeff_A( 13, vi + V2i(0,0), 6 ) += (h_inv[1]*0.0890870806375);
	ctx.coeff_A( 3, vi + V2i(0,0), 7 ) += -(h_inv[1]*0.462910049886);
	ctx.coeff_A( 3, vi + V2i(0,1), 7 ) += (h_inv[1]*0.462910049886);
	ctx.coeff_A( 10, vi + V2i(0,0), 7 ) += -(h_inv[1]*-0.471404520791);
	ctx.coeff_A( 10, vi + V2i(0,1), 7 ) += (h_inv[1]*-0.471404520791);
	ctx.coeff_A( 12, vi + V2i(0,0), 7 ) += -(h_inv[1]*-0.0890870806375);
	ctx.coeff_A( 12, vi + V2i(0,1), 7 ) += (h_inv[1]*-0.0890870806375);
	ctx.coeff_A( 4, vi + V2i(0,-1), 8 ) += -(h_inv[1]*-0.119522860933);
	ctx.coeff_A( 4, vi + V2i(0,0), 8 ) += (h_inv[1]*-0.119522860933);
	ctx.coeff_A( 13, vi + V2i(0,-1), 8 ) += -(h_inv[1]*0.345032779671);
	ctx.coeff_A( 13, vi + V2i(0,0), 8 ) += (h_inv[1]*0.345032779671);
	ctx.coeff_A( 3, vi + V2i(0,0), 9 ) += -(h_inv[1]*0.119522860933);
	ctx.coeff_A( 3, vi + V2i(0,1), 9 ) += (h_inv[1]*0.119522860933);
	ctx.coeff_A( 5, vi + V2i(0,0), 9 ) += -(h_inv[1]*0.414039335605);
	ctx.coeff_A( 5, vi + V2i(0,1), 9 ) += (h_inv[1]*0.414039335605);
	ctx.coeff_A( 12, vi + V2i(0,0), 9 ) += -(h_inv[1]*-0.345032779671);
	ctx.coeff_A( 12, vi + V2i(0,1), 9 ) += (h_inv[1]*-0.345032779671);
	ctx.coeff_A( 14, vi + V2i(0,0), 9 ) += -(h_inv[1]*-0.308606699924);
	ctx.coeff_A( 14, vi + V2i(0,1), 9 ) += (h_inv[1]*-0.308606699924);
	ctx.coeff_A( 7, vi + V2i(0,-1), 10 ) += -(h_inv[1]*-0.471404520791);
	ctx.coeff_A( 7, vi + V2i(0,0), 10 ) += (h_inv[1]*-0.471404520791);
	ctx.coeff_A( 6, vi + V2i(0,0), 11 ) += -(h_inv[1]*0.471404520791);
	ctx.coeff_A( 6, vi + V2i(0,1), 11 ) += (h_inv[1]*0.471404520791);
	ctx.coeff_A( 7, vi + V2i(0,-1), 12 ) += -(h_inv[1]*-0.0890870806375);
	ctx.coeff_A( 7, vi + V2i(0,0), 12 ) += (h_inv[1]*-0.0890870806375);
	ctx.coeff_A( 9, vi + V2i(0,-1), 12 ) += -(h_inv[1]*-0.345032779671);
	ctx.coeff_A( 9, vi + V2i(0,0), 12 ) += (h_inv[1]*-0.345032779671);
	ctx.coeff_A( 6, vi + V2i(0,0), 13 ) += -(h_inv[1]*0.0890870806375);
	ctx.coeff_A( 6, vi + V2i(0,1), 13 ) += (h_inv[1]*0.0890870806375);
	ctx.coeff_A( 8, vi + V2i(0,0), 13 ) += -(h_inv[1]*0.345032779671);
	ctx.coeff_A( 8, vi + V2i(0,1), 13 ) += (h_inv[1]*0.345032779671);
	ctx.coeff_A( 9, vi + V2i(0,-1), 14 ) += -(h_inv[1]*-0.308606699924);
	ctx.coeff_A( 9, vi + V2i(0,0), 14 ) += (h_inv[1]*-0.308606699924);
	ctx.coeff_A( 0, vi + V2i(0,0), 0 ) += M_3_real.coeffRef(0, 0);
	ctx.coeff_A( 1, vi + V2i(0,0), 1 ) += M_3_real.coeffRef(1, 1);
	ctx.coeff_A( 2, vi + V2i(0,0), 2 ) += M_3_real.coeffRef(2, 2);
	ctx.coeff_A( 3, vi + V2i(0,0), 3 ) += M_3_real.coeffRef(3, 3);
	ctx.coeff_A( 4, vi + V2i(0,0), 4 ) += M_3_real.coeffRef(4, 4);
	ctx.coeff_A( 5, vi + V2i(0,0), 5 ) += M_3_real.coeffRef(5, 5);
	ctx.coeff_A( 6, vi + V2i(0,0), 6 ) += M_3_real.coeffRef(6, 6);
	ctx.coeff_A( 7, vi + V2i(0,0), 7 ) += M_3_real.coeffRef(7, 7);
	ctx.coeff_A( 8, vi + V2i(0,0), 8 ) += M_3_real.coeffRef(8, 8);
	ctx.coeff_A( 9, vi + V2i(0,0), 9 ) += M_3_real.coeffRef(9, 9);
	ctx.coeff_A( 10, vi + V2i(0,0), 10 ) += M_3_real.coeffRef(10, 10);
	ctx.coeff_A( 11, vi + V2i(0,0), 11 ) += M_3_real.coeffRef(11, 11);
	ctx.coeff_A( 12, vi + V2i(0,0), 12 ) += M_3_real.coeffRef(12, 12);
	ctx.coeff_A( 13, vi + V2i(0,0), 13 ) += M_3_real.coeffRef(13, 13);
	ctx.coeff_A( 14, vi + V2i(0,0), 14 ) += M_3_real.coeffRef(14, 14);
	ctx.coeff_b( 0 ) += b_real.coeffRef(0, 0);
	ctx.coeff_b( 1 ) += b_real.coeffRef(1, 0);
	ctx.coeff_b( 2 ) += b_real.coeffRef(2, 0);
	ctx.coeff_b( 3 ) += b_real.coeffRef(3, 0);
	ctx.coeff_b( 4 ) += b_real.coeffRef(4, 0);
	ctx.coeff_b( 5 ) += b_real.coeffRef(5, 0);
	ctx.coeff_b( 6 ) += b_real.coeffRef(6, 0);
	ctx.coeff_b( 7 ) += b_real.coeffRef(7, 0);
	ctx.coeff_b( 8 ) += b_real.coeffRef(8, 0);
	ctx.coeff_b( 9 ) += b_real.coeffRef(9, 0);
	ctx.coeff_b( 10 ) += b_real.coeffRef(10, 0);
	ctx.coeff_b( 11 ) += b_real.coeffRef(11, 0);
	ctx.coeff_b( 12 ) += b_real.coeffRef(12, 0);
	ctx.coeff_b( 13 ) += b_real.coeffRef(13, 0);
	ctx.coeff_b( 14 ) += b_real.coeffRef(14, 0);
}
V2i stencil_fopn_p4_sg_get_offset(int coeff)
{
	switch(coeff)
	{
		case 0:return V2i(1, 1);break;
		case 1:return V2i(0, 1);break;
		case 2:return V2i(1, 0);break;
		case 3:return V2i(1, 1);break;
		case 4:return V2i(0, 0);break;
		case 5:return V2i(1, 1);break;
		case 6:return V2i(0, 1);break;
		case 7:return V2i(1, 0);break;
		case 8:return V2i(0, 1);break;
		case 9:return V2i(1, 0);break;
		case 10:return V2i(1, 1);break;
		case 11:return V2i(0, 0);break;
		case 12:return V2i(1, 1);break;
		case 13:return V2i(0, 0);break;
		case 14:return V2i(1, 1);break;
		default:throw std::runtime_error("unexpected coefficient index");break;
	};
}
REGISTER_STENCIL(stencil_fopn_p4_sg, 4, 1)
