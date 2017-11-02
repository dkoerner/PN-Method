
#include <PNSystem.h>

void stencil_cda(PNSystem::Stencil::Context& ctx)
{
	V3i vi = ctx.getVoxelCoord();
	V3d vd = vi.cast<double>();
	const Domain& domain = ctx.getDomain();
	const PNSystem::Fields& fields = ctx.getFields();
	V3d vs = domain.getVoxelSize();
	V3d h_inv( 1.0/vs[0], 1.0/vs[1], 1.0/vs[2] );
	V3d h_inv2( 1.0/(vs[0]*vs[0]), 1.0/(vs[1]*vs[1]), 1.0/(vs[2]*vs[2]) );

	V3d voxel_center_VS = vd + V3d(0.5, 0.5, 0.5);
	P3d voxel_center_WS = domain.voxelToWorld(voxel_center_VS);


	double D_xph = 1.0/(3.0*fields.sigma_t->eval(voxel_center_WS+V3d(vs[0], 0.0, 0.0)).real());
	double D_xmh = 1.0/(3.0*fields.sigma_t->eval(voxel_center_WS-V3d(vs[0], 0.0, 0.0)).real());
	double D_yph = 1.0/(3.0*fields.sigma_t->eval(voxel_center_WS+V3d(0.0, vs[1], 0.0)).real());
	double D_ymh = 1.0/(3.0*fields.sigma_t->eval(voxel_center_WS-V3d(0.0, vs[1], 0.0)).real());
	double D_zph = 1.0/(3.0*fields.sigma_t->eval(voxel_center_WS+V3d(0.0, 0.0, vs[2])).real());
	double D_zmh = 1.0/(3.0*fields.sigma_t->eval(voxel_center_WS-V3d(0.0, 0.0, vs[2])).real());

	// compute diffusion coefficient

	// Assembling global system =============
	double aii = 0.0;
	aii += -h_inv2[0]*(D_xph+D_xmh);
	aii += -h_inv2[1]*(D_yph+D_ymh);
	aii += -h_inv2[2]*(D_zph+D_zmh);
	aii += -fields.sigma_a->eval(voxel_center_WS).real();
	ctx.coeff_A( 0, vi + V3i(0,0,0), 0 ) += aii;

	ctx.coeff_A( 0, vi + V3i(1,0,0), 0 ) += h_inv2[0]*D_xph;
	ctx.coeff_A( 0, vi + V3i(-1,0,0), 0 ) += h_inv2[0]*D_xmh;
	ctx.coeff_A( 0, vi + V3i(0,1,0), 0 ) += h_inv2[1]*D_yph;
	ctx.coeff_A( 0, vi + V3i(0,-1,0), 0 ) += h_inv2[1]*D_ymh;
	ctx.coeff_A( 0, vi + V3i(0,0,1), 0 ) += h_inv2[2]*D_zph;
	ctx.coeff_A( 0, vi + V3i(0,0,-1), 0 ) += h_inv2[2]*D_zmh;

	ctx.coeff_b( 0 ) += -fields.q->eval(0,0,voxel_center_WS).real();
}
V3i stencil_cda_get_offset(int coeff)
{
	switch(coeff)
	{
		case 0:return V3i(1, 1, 1);break;
		default:throw std::runtime_error("unexpected coefficient index");break;
	};
}
REGISTER_STENCIL(stencil_cda, 0, 1, 1)
