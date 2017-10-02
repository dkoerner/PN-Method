// This file was generated by stencil.py

#include <PNSystem.h>

void stencil_fopn_p5_cg(PNSystem& sys,
					const V2i& voxel)
{
	V2i vi = voxel;
	V2d vd = vi.cast<double>();
	const Domain& domain = sys.getDomain();
	const PNSystem::Fields& fields = sys.getFields();
	V2d h_inv( 1.0/(2*domain.getVoxelSize()[0]), 1.0/(2*domain.getVoxelSize()[1]) );

	Eigen::Matrix<std::complex<double>, 21, 21> S;
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
	S.coeffRef(15, 15) = std::complex<double>(0.7071067811865475, 0.0);
	S.coeffRef(15, 20) = std::complex<double>(-0.7071067811865475, 0.0);
	S.coeffRef(16, 15) = std::complex<double>(-0.0, -0.7071067811865475);
	S.coeffRef(16, 20) = std::complex<double>(-0.0, -0.7071067811865475);
	S.coeffRef(17, 16) = std::complex<double>(0.7071067811865475, 0.0);
	S.coeffRef(17, 19) = std::complex<double>(-0.7071067811865475, 0.0);
	S.coeffRef(18, 16) = std::complex<double>(-0.0, -0.7071067811865475);
	S.coeffRef(18, 19) = std::complex<double>(-0.0, -0.7071067811865475);
	S.coeffRef(19, 17) = std::complex<double>(0.7071067811865475, 0.0);
	S.coeffRef(19, 18) = std::complex<double>(-0.7071067811865475, 0.0);
	S.coeffRef(20, 17) = std::complex<double>(-0.0, -0.7071067811865475);
	S.coeffRef(20, 18) = std::complex<double>(-0.0, -0.7071067811865475);
	Eigen::Matrix<std::complex<double>, 21, 21> SInv;
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
	SInv.coeffRef(15, 15) = std::complex<double>(0.7071067811865476, 0.0);
	SInv.coeffRef(15, 16) = std::complex<double>(0.0, 0.7071067811865476);
	SInv.coeffRef(16, 17) = std::complex<double>(0.7071067811865476, 0.0);
	SInv.coeffRef(16, 18) = std::complex<double>(0.0, 0.7071067811865476);
	SInv.coeffRef(17, 19) = std::complex<double>(0.7071067811865476, 0.0);
	SInv.coeffRef(17, 20) = std::complex<double>(0.0, 0.7071067811865476);
	SInv.coeffRef(18, 19) = std::complex<double>(-0.7071067811865476, 0.0);
	SInv.coeffRef(18, 20) = std::complex<double>(-0.0, 0.7071067811865476);
	SInv.coeffRef(19, 17) = std::complex<double>(-0.7071067811865476, 0.0);
	SInv.coeffRef(19, 18) = std::complex<double>(-0.0, 0.7071067811865476);
	SInv.coeffRef(20, 15) = std::complex<double>(-0.7071067811865476, 0.0);
	SInv.coeffRef(20, 16) = std::complex<double>(-0.0, 0.7071067811865476);

	//Producing complex-valued matrices =============
	//M_0dxL + M_1dyL + M_2dzL + M_3L = b

	//M_0 ---
	// is constant

	//M_1 ---
	// is constant

	//M_2 ---
	// all components vanish

	//M_3 ---
	Eigen::Matrix<std::complex<double>, 21, 21> M_3;
	M_3(0, 0) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(0, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(1, 1) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(2, 2) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(3, 3) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(2, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(4, 4) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(2, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(5, 5) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(2, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(6, 6) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(3, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(7, 7) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(3, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(8, 8) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(3, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(9, 9) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(3, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(10, 10) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(4, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(11, 11) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(4, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(12, 12) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(4, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(13, 13) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(4, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(14, 14) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(4, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(15, 15) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(5, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(16, 16) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(5, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(17, 17) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(5, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(18, 18) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(5, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(19, 19) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(5, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	M_3(20, 20) = (fields.sigma_t->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))+
			-(fields.sigma_s->eval(domain.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(5, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)))));
	Eigen::Matrix<double, 21, 21> M_3_real = (S*M_3*SInv).real();

	//b ---
	Eigen::Matrix<std::complex<double>, 21, 1> b;
	b(0, 0) = fields.q->eval(0, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(1, 0) = fields.q->eval(1, -1, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(2, 0) = fields.q->eval(1, 1, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(3, 0) = fields.q->eval(2, -2, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(4, 0) = fields.q->eval(2, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(5, 0) = fields.q->eval(2, 2, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(6, 0) = fields.q->eval(3, -3, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(7, 0) = fields.q->eval(3, -1, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(8, 0) = fields.q->eval(3, 1, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(9, 0) = fields.q->eval(3, 3, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(10, 0) = fields.q->eval(4, -4, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(11, 0) = fields.q->eval(4, -2, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(12, 0) = fields.q->eval(4, 0, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(13, 0) = fields.q->eval(4, 2, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(14, 0) = fields.q->eval(4, 4, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(15, 0) = fields.q->eval(5, -5, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(16, 0) = fields.q->eval(5, -3, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(17, 0) = fields.q->eval(5, -1, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(18, 0) = fields.q->eval(5, 1, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(19, 0) = fields.q->eval(5, 3, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	b(20, 0) = fields.q->eval(5, 5, domain.voxelToWorld(vd+V2d(0.5, 0.5)));
	Eigen::Matrix<double, 21, 1> b_real = (S*b).real();

	// Assembling global system =============
	sys.coeff_A( vi, 1, vi + V2i(-1,0), 0 ) += -(h_inv[0]*0.57735026919);
	sys.coeff_A( vi, 1, vi + V2i(1,0), 0 ) += (h_inv[0]*0.57735026919);
	sys.coeff_A( vi, 0, vi + V2i(-1,0), 1 ) += -(h_inv[0]*0.57735026919);
	sys.coeff_A( vi, 0, vi + V2i(1,0), 1 ) += (h_inv[0]*0.57735026919);
	sys.coeff_A( vi, 3, vi + V2i(-1,0), 1 ) += -(h_inv[0]*0.4472135955);
	sys.coeff_A( vi, 3, vi + V2i(1,0), 1 ) += (h_inv[0]*0.4472135955);
	sys.coeff_A( vi, 5, vi + V2i(-1,0), 1 ) += -(h_inv[0]*-0.258198889747);
	sys.coeff_A( vi, 5, vi + V2i(1,0), 1 ) += (h_inv[0]*-0.258198889747);
	sys.coeff_A( vi, 4, vi + V2i(-1,0), 2 ) += -(h_inv[0]*0.4472135955);
	sys.coeff_A( vi, 4, vi + V2i(1,0), 2 ) += (h_inv[0]*0.4472135955);
	sys.coeff_A( vi, 1, vi + V2i(-1,0), 3 ) += -(h_inv[0]*0.4472135955);
	sys.coeff_A( vi, 1, vi + V2i(1,0), 3 ) += (h_inv[0]*0.4472135955);
	sys.coeff_A( vi, 6, vi + V2i(-1,0), 3 ) += -(h_inv[0]*0.462910049886);
	sys.coeff_A( vi, 6, vi + V2i(1,0), 3 ) += (h_inv[0]*0.462910049886);
	sys.coeff_A( vi, 8, vi + V2i(-1,0), 3 ) += -(h_inv[0]*-0.119522860933);
	sys.coeff_A( vi, 8, vi + V2i(1,0), 3 ) += (h_inv[0]*-0.119522860933);
	sys.coeff_A( vi, 2, vi + V2i(-1,0), 4 ) += -(h_inv[0]*0.4472135955);
	sys.coeff_A( vi, 2, vi + V2i(1,0), 4 ) += (h_inv[0]*0.4472135955);
	sys.coeff_A( vi, 7, vi + V2i(-1,0), 4 ) += -(h_inv[0]*0.462910049886);
	sys.coeff_A( vi, 7, vi + V2i(1,0), 4 ) += (h_inv[0]*0.462910049886);
	sys.coeff_A( vi, 9, vi + V2i(-1,0), 4 ) += -(h_inv[0]*-0.119522860933);
	sys.coeff_A( vi, 9, vi + V2i(1,0), 4 ) += (h_inv[0]*-0.119522860933);
	sys.coeff_A( vi, 1, vi + V2i(-1,0), 5 ) += -(h_inv[0]*-0.258198889747);
	sys.coeff_A( vi, 1, vi + V2i(1,0), 5 ) += (h_inv[0]*-0.258198889747);
	sys.coeff_A( vi, 8, vi + V2i(-1,0), 5 ) += -(h_inv[0]*0.414039335605);
	sys.coeff_A( vi, 8, vi + V2i(1,0), 5 ) += (h_inv[0]*0.414039335605);
	sys.coeff_A( vi, 3, vi + V2i(-1,0), 6 ) += -(h_inv[0]*0.462910049886);
	sys.coeff_A( vi, 3, vi + V2i(1,0), 6 ) += (h_inv[0]*0.462910049886);
	sys.coeff_A( vi, 10, vi + V2i(-1,0), 6 ) += -(h_inv[0]*0.471404520791);
	sys.coeff_A( vi, 10, vi + V2i(1,0), 6 ) += (h_inv[0]*0.471404520791);
	sys.coeff_A( vi, 12, vi + V2i(-1,0), 6 ) += -(h_inv[0]*-0.0890870806375);
	sys.coeff_A( vi, 12, vi + V2i(1,0), 6 ) += (h_inv[0]*-0.0890870806375);
	sys.coeff_A( vi, 4, vi + V2i(-1,0), 7 ) += -(h_inv[0]*0.462910049886);
	sys.coeff_A( vi, 4, vi + V2i(1,0), 7 ) += (h_inv[0]*0.462910049886);
	sys.coeff_A( vi, 11, vi + V2i(-1,0), 7 ) += -(h_inv[0]*0.471404520791);
	sys.coeff_A( vi, 11, vi + V2i(1,0), 7 ) += (h_inv[0]*0.471404520791);
	sys.coeff_A( vi, 13, vi + V2i(-1,0), 7 ) += -(h_inv[0]*-0.0890870806375);
	sys.coeff_A( vi, 13, vi + V2i(1,0), 7 ) += (h_inv[0]*-0.0890870806375);
	sys.coeff_A( vi, 3, vi + V2i(-1,0), 8 ) += -(h_inv[0]*-0.119522860933);
	sys.coeff_A( vi, 3, vi + V2i(1,0), 8 ) += (h_inv[0]*-0.119522860933);
	sys.coeff_A( vi, 5, vi + V2i(-1,0), 8 ) += -(h_inv[0]*0.414039335605);
	sys.coeff_A( vi, 5, vi + V2i(1,0), 8 ) += (h_inv[0]*0.414039335605);
	sys.coeff_A( vi, 12, vi + V2i(-1,0), 8 ) += -(h_inv[0]*0.345032779671);
	sys.coeff_A( vi, 12, vi + V2i(1,0), 8 ) += (h_inv[0]*0.345032779671);
	sys.coeff_A( vi, 14, vi + V2i(-1,0), 8 ) += -(h_inv[0]*-0.308606699924);
	sys.coeff_A( vi, 14, vi + V2i(1,0), 8 ) += (h_inv[0]*-0.308606699924);
	sys.coeff_A( vi, 4, vi + V2i(-1,0), 9 ) += -(h_inv[0]*-0.119522860933);
	sys.coeff_A( vi, 4, vi + V2i(1,0), 9 ) += (h_inv[0]*-0.119522860933);
	sys.coeff_A( vi, 13, vi + V2i(-1,0), 9 ) += -(h_inv[0]*0.345032779671);
	sys.coeff_A( vi, 13, vi + V2i(1,0), 9 ) += (h_inv[0]*0.345032779671);
	sys.coeff_A( vi, 6, vi + V2i(-1,0), 10 ) += -(h_inv[0]*0.471404520791);
	sys.coeff_A( vi, 6, vi + V2i(1,0), 10 ) += (h_inv[0]*0.471404520791);
	sys.coeff_A( vi, 15, vi + V2i(-1,0), 10 ) += -(h_inv[0]*0.476731294623);
	sys.coeff_A( vi, 15, vi + V2i(1,0), 10 ) += (h_inv[0]*0.476731294623);
	sys.coeff_A( vi, 17, vi + V2i(-1,0), 10 ) += -(h_inv[0]*-0.0710669054519);
	sys.coeff_A( vi, 17, vi + V2i(1,0), 10 ) += (h_inv[0]*-0.0710669054519);
	sys.coeff_A( vi, 7, vi + V2i(-1,0), 11 ) += -(h_inv[0]*0.471404520791);
	sys.coeff_A( vi, 7, vi + V2i(1,0), 11 ) += (h_inv[0]*0.471404520791);
	sys.coeff_A( vi, 16, vi + V2i(-1,0), 11 ) += -(h_inv[0]*0.476731294623);
	sys.coeff_A( vi, 16, vi + V2i(1,0), 11 ) += (h_inv[0]*0.476731294623);
	sys.coeff_A( vi, 18, vi + V2i(-1,0), 11 ) += -(h_inv[0]*-0.0710669054519);
	sys.coeff_A( vi, 18, vi + V2i(1,0), 11 ) += (h_inv[0]*-0.0710669054519);
	sys.coeff_A( vi, 6, vi + V2i(-1,0), 12 ) += -(h_inv[0]*-0.0890870806375);
	sys.coeff_A( vi, 6, vi + V2i(1,0), 12 ) += (h_inv[0]*-0.0890870806375);
	sys.coeff_A( vi, 8, vi + V2i(-1,0), 12 ) += -(h_inv[0]*0.345032779671);
	sys.coeff_A( vi, 8, vi + V2i(1,0), 12 ) += (h_inv[0]*0.345032779671);
	sys.coeff_A( vi, 17, vi + V2i(-1,0), 12 ) += -(h_inv[0]*0.376050716545);
	sys.coeff_A( vi, 17, vi + V2i(1,0), 12 ) += (h_inv[0]*0.376050716545);
	sys.coeff_A( vi, 19, vi + V2i(-1,0), 12 ) += -(h_inv[0]*-0.174077655956);
	sys.coeff_A( vi, 19, vi + V2i(1,0), 12 ) += (h_inv[0]*-0.174077655956);
	sys.coeff_A( vi, 7, vi + V2i(-1,0), 13 ) += -(h_inv[0]*-0.0890870806375);
	sys.coeff_A( vi, 7, vi + V2i(1,0), 13 ) += (h_inv[0]*-0.0890870806375);
	sys.coeff_A( vi, 9, vi + V2i(-1,0), 13 ) += -(h_inv[0]*0.345032779671);
	sys.coeff_A( vi, 9, vi + V2i(1,0), 13 ) += (h_inv[0]*0.345032779671);
	sys.coeff_A( vi, 18, vi + V2i(-1,0), 13 ) += -(h_inv[0]*0.376050716545);
	sys.coeff_A( vi, 18, vi + V2i(1,0), 13 ) += (h_inv[0]*0.376050716545);
	sys.coeff_A( vi, 20, vi + V2i(-1,0), 13 ) += -(h_inv[0]*-0.174077655956);
	sys.coeff_A( vi, 20, vi + V2i(1,0), 13 ) += (h_inv[0]*-0.174077655956);
	sys.coeff_A( vi, 8, vi + V2i(-1,0), 14 ) += -(h_inv[0]*-0.308606699924);
	sys.coeff_A( vi, 8, vi + V2i(1,0), 14 ) += (h_inv[0]*-0.308606699924);
	sys.coeff_A( vi, 19, vi + V2i(-1,0), 14 ) += -(h_inv[0]*0.389249472081);
	sys.coeff_A( vi, 19, vi + V2i(1,0), 14 ) += (h_inv[0]*0.389249472081);
	sys.coeff_A( vi, 10, vi + V2i(-1,0), 15 ) += -(h_inv[0]*0.476731294623);
	sys.coeff_A( vi, 10, vi + V2i(1,0), 15 ) += (h_inv[0]*0.476731294623);
	sys.coeff_A( vi, 11, vi + V2i(-1,0), 16 ) += -(h_inv[0]*0.476731294623);
	sys.coeff_A( vi, 11, vi + V2i(1,0), 16 ) += (h_inv[0]*0.476731294623);
	sys.coeff_A( vi, 10, vi + V2i(-1,0), 17 ) += -(h_inv[0]*-0.0710669054519);
	sys.coeff_A( vi, 10, vi + V2i(1,0), 17 ) += (h_inv[0]*-0.0710669054519);
	sys.coeff_A( vi, 12, vi + V2i(-1,0), 17 ) += -(h_inv[0]*0.376050716545);
	sys.coeff_A( vi, 12, vi + V2i(1,0), 17 ) += (h_inv[0]*0.376050716545);
	sys.coeff_A( vi, 11, vi + V2i(-1,0), 18 ) += -(h_inv[0]*-0.0710669054519);
	sys.coeff_A( vi, 11, vi + V2i(1,0), 18 ) += (h_inv[0]*-0.0710669054519);
	sys.coeff_A( vi, 13, vi + V2i(-1,0), 18 ) += -(h_inv[0]*0.376050716545);
	sys.coeff_A( vi, 13, vi + V2i(1,0), 18 ) += (h_inv[0]*0.376050716545);
	sys.coeff_A( vi, 12, vi + V2i(-1,0), 19 ) += -(h_inv[0]*-0.174077655956);
	sys.coeff_A( vi, 12, vi + V2i(1,0), 19 ) += (h_inv[0]*-0.174077655956);
	sys.coeff_A( vi, 14, vi + V2i(-1,0), 19 ) += -(h_inv[0]*0.389249472081);
	sys.coeff_A( vi, 14, vi + V2i(1,0), 19 ) += (h_inv[0]*0.389249472081);
	sys.coeff_A( vi, 13, vi + V2i(-1,0), 20 ) += -(h_inv[0]*-0.174077655956);
	sys.coeff_A( vi, 13, vi + V2i(1,0), 20 ) += (h_inv[0]*-0.174077655956);
	sys.coeff_A( vi, 2, vi + V2i(0,-1), 0 ) += -(h_inv[1]*0.57735026919);
	sys.coeff_A( vi, 2, vi + V2i(0,1), 0 ) += (h_inv[1]*0.57735026919);
	sys.coeff_A( vi, 4, vi + V2i(0,-1), 1 ) += -(h_inv[1]*0.4472135955);
	sys.coeff_A( vi, 4, vi + V2i(0,1), 1 ) += (h_inv[1]*0.4472135955);
	sys.coeff_A( vi, 0, vi + V2i(0,-1), 2 ) += -(h_inv[1]*0.57735026919);
	sys.coeff_A( vi, 0, vi + V2i(0,1), 2 ) += (h_inv[1]*0.57735026919);
	sys.coeff_A( vi, 3, vi + V2i(0,-1), 2 ) += -(h_inv[1]*-0.4472135955);
	sys.coeff_A( vi, 3, vi + V2i(0,1), 2 ) += (h_inv[1]*-0.4472135955);
	sys.coeff_A( vi, 5, vi + V2i(0,-1), 2 ) += -(h_inv[1]*-0.258198889747);
	sys.coeff_A( vi, 5, vi + V2i(0,1), 2 ) += (h_inv[1]*-0.258198889747);
	sys.coeff_A( vi, 2, vi + V2i(0,-1), 3 ) += -(h_inv[1]*-0.4472135955);
	sys.coeff_A( vi, 2, vi + V2i(0,1), 3 ) += (h_inv[1]*-0.4472135955);
	sys.coeff_A( vi, 7, vi + V2i(0,-1), 3 ) += -(h_inv[1]*0.462910049886);
	sys.coeff_A( vi, 7, vi + V2i(0,1), 3 ) += (h_inv[1]*0.462910049886);
	sys.coeff_A( vi, 9, vi + V2i(0,-1), 3 ) += -(h_inv[1]*0.119522860933);
	sys.coeff_A( vi, 9, vi + V2i(0,1), 3 ) += (h_inv[1]*0.119522860933);
	sys.coeff_A( vi, 1, vi + V2i(0,-1), 4 ) += -(h_inv[1]*0.4472135955);
	sys.coeff_A( vi, 1, vi + V2i(0,1), 4 ) += (h_inv[1]*0.4472135955);
	sys.coeff_A( vi, 6, vi + V2i(0,-1), 4 ) += -(h_inv[1]*-0.462910049886);
	sys.coeff_A( vi, 6, vi + V2i(0,1), 4 ) += (h_inv[1]*-0.462910049886);
	sys.coeff_A( vi, 8, vi + V2i(0,-1), 4 ) += -(h_inv[1]*-0.119522860933);
	sys.coeff_A( vi, 8, vi + V2i(0,1), 4 ) += (h_inv[1]*-0.119522860933);
	sys.coeff_A( vi, 2, vi + V2i(0,-1), 5 ) += -(h_inv[1]*-0.258198889747);
	sys.coeff_A( vi, 2, vi + V2i(0,1), 5 ) += (h_inv[1]*-0.258198889747);
	sys.coeff_A( vi, 9, vi + V2i(0,-1), 5 ) += -(h_inv[1]*0.414039335605);
	sys.coeff_A( vi, 9, vi + V2i(0,1), 5 ) += (h_inv[1]*0.414039335605);
	sys.coeff_A( vi, 4, vi + V2i(0,-1), 6 ) += -(h_inv[1]*-0.462910049886);
	sys.coeff_A( vi, 4, vi + V2i(0,1), 6 ) += (h_inv[1]*-0.462910049886);
	sys.coeff_A( vi, 11, vi + V2i(0,-1), 6 ) += -(h_inv[1]*0.471404520791);
	sys.coeff_A( vi, 11, vi + V2i(0,1), 6 ) += (h_inv[1]*0.471404520791);
	sys.coeff_A( vi, 13, vi + V2i(0,-1), 6 ) += -(h_inv[1]*0.0890870806375);
	sys.coeff_A( vi, 13, vi + V2i(0,1), 6 ) += (h_inv[1]*0.0890870806375);
	sys.coeff_A( vi, 3, vi + V2i(0,-1), 7 ) += -(h_inv[1]*0.462910049886);
	sys.coeff_A( vi, 3, vi + V2i(0,1), 7 ) += (h_inv[1]*0.462910049886);
	sys.coeff_A( vi, 10, vi + V2i(0,-1), 7 ) += -(h_inv[1]*-0.471404520791);
	sys.coeff_A( vi, 10, vi + V2i(0,1), 7 ) += (h_inv[1]*-0.471404520791);
	sys.coeff_A( vi, 12, vi + V2i(0,-1), 7 ) += -(h_inv[1]*-0.0890870806375);
	sys.coeff_A( vi, 12, vi + V2i(0,1), 7 ) += (h_inv[1]*-0.0890870806375);
	sys.coeff_A( vi, 4, vi + V2i(0,-1), 8 ) += -(h_inv[1]*-0.119522860933);
	sys.coeff_A( vi, 4, vi + V2i(0,1), 8 ) += (h_inv[1]*-0.119522860933);
	sys.coeff_A( vi, 13, vi + V2i(0,-1), 8 ) += -(h_inv[1]*0.345032779671);
	sys.coeff_A( vi, 13, vi + V2i(0,1), 8 ) += (h_inv[1]*0.345032779671);
	sys.coeff_A( vi, 3, vi + V2i(0,-1), 9 ) += -(h_inv[1]*0.119522860933);
	sys.coeff_A( vi, 3, vi + V2i(0,1), 9 ) += (h_inv[1]*0.119522860933);
	sys.coeff_A( vi, 5, vi + V2i(0,-1), 9 ) += -(h_inv[1]*0.414039335605);
	sys.coeff_A( vi, 5, vi + V2i(0,1), 9 ) += (h_inv[1]*0.414039335605);
	sys.coeff_A( vi, 12, vi + V2i(0,-1), 9 ) += -(h_inv[1]*-0.345032779671);
	sys.coeff_A( vi, 12, vi + V2i(0,1), 9 ) += (h_inv[1]*-0.345032779671);
	sys.coeff_A( vi, 14, vi + V2i(0,-1), 9 ) += -(h_inv[1]*-0.308606699924);
	sys.coeff_A( vi, 14, vi + V2i(0,1), 9 ) += (h_inv[1]*-0.308606699924);
	sys.coeff_A( vi, 7, vi + V2i(0,-1), 10 ) += -(h_inv[1]*-0.471404520791);
	sys.coeff_A( vi, 7, vi + V2i(0,1), 10 ) += (h_inv[1]*-0.471404520791);
	sys.coeff_A( vi, 16, vi + V2i(0,-1), 10 ) += -(h_inv[1]*0.476731294623);
	sys.coeff_A( vi, 16, vi + V2i(0,1), 10 ) += (h_inv[1]*0.476731294623);
	sys.coeff_A( vi, 18, vi + V2i(0,-1), 10 ) += -(h_inv[1]*0.0710669054519);
	sys.coeff_A( vi, 18, vi + V2i(0,1), 10 ) += (h_inv[1]*0.0710669054519);
	sys.coeff_A( vi, 6, vi + V2i(0,-1), 11 ) += -(h_inv[1]*0.471404520791);
	sys.coeff_A( vi, 6, vi + V2i(0,1), 11 ) += (h_inv[1]*0.471404520791);
	sys.coeff_A( vi, 15, vi + V2i(0,-1), 11 ) += -(h_inv[1]*-0.476731294623);
	sys.coeff_A( vi, 15, vi + V2i(0,1), 11 ) += (h_inv[1]*-0.476731294623);
	sys.coeff_A( vi, 17, vi + V2i(0,-1), 11 ) += -(h_inv[1]*-0.0710669054519);
	sys.coeff_A( vi, 17, vi + V2i(0,1), 11 ) += (h_inv[1]*-0.0710669054519);
	sys.coeff_A( vi, 7, vi + V2i(0,-1), 12 ) += -(h_inv[1]*-0.0890870806375);
	sys.coeff_A( vi, 7, vi + V2i(0,1), 12 ) += (h_inv[1]*-0.0890870806375);
	sys.coeff_A( vi, 9, vi + V2i(0,-1), 12 ) += -(h_inv[1]*-0.345032779671);
	sys.coeff_A( vi, 9, vi + V2i(0,1), 12 ) += (h_inv[1]*-0.345032779671);
	sys.coeff_A( vi, 18, vi + V2i(0,-1), 12 ) += -(h_inv[1]*0.376050716545);
	sys.coeff_A( vi, 18, vi + V2i(0,1), 12 ) += (h_inv[1]*0.376050716545);
	sys.coeff_A( vi, 20, vi + V2i(0,-1), 12 ) += -(h_inv[1]*0.174077655956);
	sys.coeff_A( vi, 20, vi + V2i(0,1), 12 ) += (h_inv[1]*0.174077655956);
	sys.coeff_A( vi, 6, vi + V2i(0,-1), 13 ) += -(h_inv[1]*0.0890870806375);
	sys.coeff_A( vi, 6, vi + V2i(0,1), 13 ) += (h_inv[1]*0.0890870806375);
	sys.coeff_A( vi, 8, vi + V2i(0,-1), 13 ) += -(h_inv[1]*0.345032779671);
	sys.coeff_A( vi, 8, vi + V2i(0,1), 13 ) += (h_inv[1]*0.345032779671);
	sys.coeff_A( vi, 17, vi + V2i(0,-1), 13 ) += -(h_inv[1]*-0.376050716545);
	sys.coeff_A( vi, 17, vi + V2i(0,1), 13 ) += (h_inv[1]*-0.376050716545);
	sys.coeff_A( vi, 19, vi + V2i(0,-1), 13 ) += -(h_inv[1]*-0.174077655956);
	sys.coeff_A( vi, 19, vi + V2i(0,1), 13 ) += (h_inv[1]*-0.174077655956);
	sys.coeff_A( vi, 9, vi + V2i(0,-1), 14 ) += -(h_inv[1]*-0.308606699924);
	sys.coeff_A( vi, 9, vi + V2i(0,1), 14 ) += (h_inv[1]*-0.308606699924);
	sys.coeff_A( vi, 20, vi + V2i(0,-1), 14 ) += -(h_inv[1]*0.389249472081);
	sys.coeff_A( vi, 20, vi + V2i(0,1), 14 ) += (h_inv[1]*0.389249472081);
	sys.coeff_A( vi, 11, vi + V2i(0,-1), 15 ) += -(h_inv[1]*-0.476731294623);
	sys.coeff_A( vi, 11, vi + V2i(0,1), 15 ) += (h_inv[1]*-0.476731294623);
	sys.coeff_A( vi, 10, vi + V2i(0,-1), 16 ) += -(h_inv[1]*0.476731294623);
	sys.coeff_A( vi, 10, vi + V2i(0,1), 16 ) += (h_inv[1]*0.476731294623);
	sys.coeff_A( vi, 11, vi + V2i(0,-1), 17 ) += -(h_inv[1]*-0.0710669054519);
	sys.coeff_A( vi, 11, vi + V2i(0,1), 17 ) += (h_inv[1]*-0.0710669054519);
	sys.coeff_A( vi, 13, vi + V2i(0,-1), 17 ) += -(h_inv[1]*-0.376050716545);
	sys.coeff_A( vi, 13, vi + V2i(0,1), 17 ) += (h_inv[1]*-0.376050716545);
	sys.coeff_A( vi, 10, vi + V2i(0,-1), 18 ) += -(h_inv[1]*0.0710669054519);
	sys.coeff_A( vi, 10, vi + V2i(0,1), 18 ) += (h_inv[1]*0.0710669054519);
	sys.coeff_A( vi, 12, vi + V2i(0,-1), 18 ) += -(h_inv[1]*0.376050716545);
	sys.coeff_A( vi, 12, vi + V2i(0,1), 18 ) += (h_inv[1]*0.376050716545);
	sys.coeff_A( vi, 13, vi + V2i(0,-1), 19 ) += -(h_inv[1]*-0.174077655956);
	sys.coeff_A( vi, 13, vi + V2i(0,1), 19 ) += (h_inv[1]*-0.174077655956);
	sys.coeff_A( vi, 12, vi + V2i(0,-1), 20 ) += -(h_inv[1]*0.174077655956);
	sys.coeff_A( vi, 12, vi + V2i(0,1), 20 ) += (h_inv[1]*0.174077655956);
	sys.coeff_A( vi, 14, vi + V2i(0,-1), 20 ) += -(h_inv[1]*0.389249472081);
	sys.coeff_A( vi, 14, vi + V2i(0,1), 20 ) += (h_inv[1]*0.389249472081);
	sys.coeff_A( vi, 0, vi + V2i(0,0), 0 ) += M_3_real.coeffRef(0, 0);
	sys.coeff_A( vi, 1, vi + V2i(0,0), 1 ) += M_3_real.coeffRef(1, 1);
	sys.coeff_A( vi, 2, vi + V2i(0,0), 2 ) += M_3_real.coeffRef(2, 2);
	sys.coeff_A( vi, 3, vi + V2i(0,0), 3 ) += M_3_real.coeffRef(3, 3);
	sys.coeff_A( vi, 4, vi + V2i(0,0), 4 ) += M_3_real.coeffRef(4, 4);
	sys.coeff_A( vi, 5, vi + V2i(0,0), 5 ) += M_3_real.coeffRef(5, 5);
	sys.coeff_A( vi, 6, vi + V2i(0,0), 6 ) += M_3_real.coeffRef(6, 6);
	sys.coeff_A( vi, 7, vi + V2i(0,0), 7 ) += M_3_real.coeffRef(7, 7);
	sys.coeff_A( vi, 8, vi + V2i(0,0), 8 ) += M_3_real.coeffRef(8, 8);
	sys.coeff_A( vi, 9, vi + V2i(0,0), 9 ) += M_3_real.coeffRef(9, 9);
	sys.coeff_A( vi, 10, vi + V2i(0,0), 10 ) += M_3_real.coeffRef(10, 10);
	sys.coeff_A( vi, 11, vi + V2i(0,0), 11 ) += M_3_real.coeffRef(11, 11);
	sys.coeff_A( vi, 12, vi + V2i(0,0), 12 ) += M_3_real.coeffRef(12, 12);
	sys.coeff_A( vi, 13, vi + V2i(0,0), 13 ) += M_3_real.coeffRef(13, 13);
	sys.coeff_A( vi, 14, vi + V2i(0,0), 14 ) += M_3_real.coeffRef(14, 14);
	sys.coeff_A( vi, 15, vi + V2i(0,0), 15 ) += M_3_real.coeffRef(15, 15);
	sys.coeff_A( vi, 16, vi + V2i(0,0), 16 ) += M_3_real.coeffRef(16, 16);
	sys.coeff_A( vi, 17, vi + V2i(0,0), 17 ) += M_3_real.coeffRef(17, 17);
	sys.coeff_A( vi, 18, vi + V2i(0,0), 18 ) += M_3_real.coeffRef(18, 18);
	sys.coeff_A( vi, 19, vi + V2i(0,0), 19 ) += M_3_real.coeffRef(19, 19);
	sys.coeff_A( vi, 20, vi + V2i(0,0), 20 ) += M_3_real.coeffRef(20, 20);
	sys.coeff_b( vi, 0 ) += b_real.coeffRef(0, 0);
	sys.coeff_b( vi, 1 ) += b_real.coeffRef(1, 0);
	sys.coeff_b( vi, 2 ) += b_real.coeffRef(2, 0);
	sys.coeff_b( vi, 3 ) += b_real.coeffRef(3, 0);
	sys.coeff_b( vi, 4 ) += b_real.coeffRef(4, 0);
	sys.coeff_b( vi, 5 ) += b_real.coeffRef(5, 0);
	sys.coeff_b( vi, 6 ) += b_real.coeffRef(6, 0);
	sys.coeff_b( vi, 7 ) += b_real.coeffRef(7, 0);
	sys.coeff_b( vi, 8 ) += b_real.coeffRef(8, 0);
	sys.coeff_b( vi, 9 ) += b_real.coeffRef(9, 0);
	sys.coeff_b( vi, 10 ) += b_real.coeffRef(10, 0);
	sys.coeff_b( vi, 11 ) += b_real.coeffRef(11, 0);
	sys.coeff_b( vi, 12 ) += b_real.coeffRef(12, 0);
	sys.coeff_b( vi, 13 ) += b_real.coeffRef(13, 0);
	sys.coeff_b( vi, 14 ) += b_real.coeffRef(14, 0);
	sys.coeff_b( vi, 15 ) += b_real.coeffRef(15, 0);
	sys.coeff_b( vi, 16 ) += b_real.coeffRef(16, 0);
	sys.coeff_b( vi, 17 ) += b_real.coeffRef(17, 0);
	sys.coeff_b( vi, 18 ) += b_real.coeffRef(18, 0);
	sys.coeff_b( vi, 19 ) += b_real.coeffRef(19, 0);
	sys.coeff_b( vi, 20 ) += b_real.coeffRef(20, 0);
}
V2i stencil_fopn_p5_cg_get_offset(int coeff)
{
	switch(coeff)
	{
		case 0:return V2i(1, 1);break;
		case 1:return V2i(1, 1);break;
		case 2:return V2i(1, 1);break;
		case 3:return V2i(1, 1);break;
		case 4:return V2i(1, 1);break;
		case 5:return V2i(1, 1);break;
		case 6:return V2i(1, 1);break;
		case 7:return V2i(1, 1);break;
		case 8:return V2i(1, 1);break;
		case 9:return V2i(1, 1);break;
		case 10:return V2i(1, 1);break;
		case 11:return V2i(1, 1);break;
		case 12:return V2i(1, 1);break;
		case 13:return V2i(1, 1);break;
		case 14:return V2i(1, 1);break;
		case 15:return V2i(1, 1);break;
		case 16:return V2i(1, 1);break;
		case 17:return V2i(1, 1);break;
		case 18:return V2i(1, 1);break;
		case 19:return V2i(1, 1);break;
		case 20:return V2i(1, 1);break;
		default:throw std::runtime_error("unexpected coefficient index");break;
	};
}
REGISTER_STENCIL(stencil_fopn_p5_cg, 5, 1)
