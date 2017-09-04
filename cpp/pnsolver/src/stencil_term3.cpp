#include <PNSystem.h>

void set_system_row(PNSystem::Voxel& sys,
					PNSystem::Fields& fields)
{
	V2i vi = sys.getVoxel();
	V2d vd = sys.getVoxel().cast<double>();

	// term 0 ----------------
	sys.A( 2, vi + V2i(0,0), 0 ) += -(0.204124145232*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(-0.5, 0.5))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));
	sys.A( 2, vi + V2i(-1,0), 0 ) += -(0.204124145232*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(-0.5, 0.5))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));

	// term 1 ----------------
	sys.A( 2, vi + V2i(0,0), 0 ) += -(0.204124145232*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*(-(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(-0.5, 0.5))))+
			(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))))));
	sys.A( 2, vi + V2i(-1,0), 0 ) += -(0.204124145232*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*(-(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(-0.5, 0.5))))+
			(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))))));

	// term 2 ----------------
	sys.A( 2, vi + V2i(-1,0), 0 ) += (4.08248290464*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));
	sys.A( 2, vi + V2i(0,0), 0 ) += -(4.08248290464*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));

	// term 3 ----------------
	sys.A( 0, vi + V2i(0,1), 1 ) += (0.204124145232*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(1.0, 0.5)))))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(0,0), 1 ) += (0.204124145232*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(1.0, 0.5)))))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));

	// term 4 ----------------
	sys.A( 0, vi + V2i(0,1), 1 ) += (0.204124145232*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*(-(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(1.0, 0.5))))));
	sys.A( 0, vi + V2i(0,0), 1 ) += (0.204124145232*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*(-(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(1.0, 0.5))))));

	// term 5 ----------------
	sys.A( 0, vi + V2i(-1,0), 1 ) += -(1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(-1,1), 1 ) += -(1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(0,0), 1 ) += -(1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(0,1), 1 ) += -(1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(0,0), 1 ) += (1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(0,1), 1 ) += (1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(1,0), 1 ) += (1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(1,1), 1 ) += (1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));

	// term 6 ----------------
	sys.A( 1, vi + V2i(0,0), 0 ) += (0.204124145232*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(1.0, 0.0)))))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));
	sys.A( 1, vi + V2i(0,-1), 0 ) += (0.204124145232*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(1.0, 0.0)))))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));

	// term 7 ----------------
	sys.A( 1, vi + V2i(0,0), 0 ) += (0.204124145232*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*(-(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(1.0, 0.0))))));
	sys.A( 1, vi + V2i(0,-1), 0 ) += (0.204124145232*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*(-(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(1.0, 0.0))))));

	// term 8 ----------------
	sys.A( 1, vi + V2i(-1,-1), 0 ) += -(1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));
	sys.A( 1, vi + V2i(-1,0), 0 ) += -(1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));
	sys.A( 1, vi + V2i(0,-1), 0 ) += -(1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));
	sys.A( 1, vi + V2i(0,0), 0 ) += -(1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));
	sys.A( 1, vi + V2i(0,-1), 0 ) += (1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));
	sys.A( 1, vi + V2i(0,0), 0 ) += (1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));
	sys.A( 1, vi + V2i(1,-1), 0 ) += (1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));
	sys.A( 1, vi + V2i(1,0), 0 ) += (1.02062072616*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));

	// term 9 ----------------
	sys.A( 0, vi + V2i(1,0), 2 ) += -(0.204124145232*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(1.0, 0.5)))))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(0,0), 2 ) += -(0.204124145232*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(1.0, 0.5)))))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));

	// term 10 ----------------
	sys.A( 0, vi + V2i(1,0), 2 ) += -(0.204124145232*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*(-(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(1.0, 0.5))))));
	sys.A( 0, vi + V2i(0,0), 2 ) += -(0.204124145232*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*(-(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))))+
			(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(1.0, 0.5))))));

	// term 11 ----------------
	sys.A( 0, vi + V2i(0,0), 2 ) += (4.08248290464*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(1,0), 2 ) += -(4.08248290464*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));

	// term 12 ----------------
	sys.A( 2, vi + V2i(0,0), 0 ) += (std::complex<double>(0.0, 0.2041241452319315)*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 1.0)))))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));
	sys.A( 2, vi + V2i(-1,0), 0 ) += (std::complex<double>(0.0, 0.2041241452319315)*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 1.0)))))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));

	// term 13 ----------------
	sys.A( 2, vi + V2i(0,0), 0 ) += (std::complex<double>(0.0, 0.2041241452319315)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*(-(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 1.0))))));
	sys.A( 2, vi + V2i(-1,0), 0 ) += (std::complex<double>(0.0, 0.2041241452319315)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*(-(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.0))))+
			(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 1.0))))));

	// term 14 ----------------
	sys.A( 2, vi + V2i(-1,-1), 0 ) += -(std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));
	sys.A( 2, vi + V2i(-1,0), 0 ) += -(std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));
	sys.A( 2, vi + V2i(0,-1), 0 ) += -(std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));
	sys.A( 2, vi + V2i(0,0), 0 ) += -(std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));
	sys.A( 2, vi + V2i(-1,0), 0 ) += (std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));
	sys.A( 2, vi + V2i(-1,1), 0 ) += (std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));
	sys.A( 2, vi + V2i(0,0), 0 ) += (std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));
	sys.A( 2, vi + V2i(0,1), 0 ) += (std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.0, 0.5)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.0, 0.5))));

	// term 15 ----------------
	sys.A( 0, vi + V2i(0,1), 1 ) += -(std::complex<double>(0.0, 0.2041241452319315)*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 1.0)))))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(0,0), 1 ) += -(std::complex<double>(0.0, 0.2041241452319315)*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 1.0)))))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));

	// term 16 ----------------
	sys.A( 0, vi + V2i(0,1), 1 ) += -(std::complex<double>(0.0, 0.2041241452319315)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*(-(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 1.0))))));
	sys.A( 0, vi + V2i(0,0), 1 ) += -(std::complex<double>(0.0, 0.2041241452319315)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*(-(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 1.0))))));

	// term 17 ----------------
	sys.A( 0, vi + V2i(0,0), 1 ) += (std::complex<double>(0.0, 4.08248290463863)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(0,1), 1 ) += -(std::complex<double>(0.0, 4.08248290463863)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));

	// term 18 ----------------
	sys.A( 1, vi + V2i(0,0), 0 ) += (std::complex<double>(0.0, 0.2041241452319315)*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, -0.5))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));
	sys.A( 1, vi + V2i(0,-1), 0 ) += (std::complex<double>(0.0, 0.2041241452319315)*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, -0.5))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));

	// term 19 ----------------
	sys.A( 1, vi + V2i(0,0), 0 ) += (std::complex<double>(0.0, 0.2041241452319315)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*(-(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, -0.5))))+
			(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))))));
	sys.A( 1, vi + V2i(0,-1), 0 ) += (std::complex<double>(0.0, 0.2041241452319315)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*(-(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, -0.5))))+
			(10.0*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))))));

	// term 20 ----------------
	sys.A( 1, vi + V2i(0,-1), 0 ) += -(std::complex<double>(0.0, 4.08248290463863)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));
	sys.A( 1, vi + V2i(0,0), 0 ) += (std::complex<double>(0.0, 4.08248290463863)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0)))*fields.f_p->eval(0, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))));

	// term 21 ----------------
	sys.A( 0, vi + V2i(1,0), 2 ) += -(std::complex<double>(0.0, 0.2041241452319315)*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 1.0)))))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(0,0), 2 ) += -(std::complex<double>(0.0, 0.2041241452319315)*(-(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(10.0*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 1.0)))))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));

	// term 22 ----------------
	sys.A( 0, vi + V2i(1,0), 2 ) += -(std::complex<double>(0.0, 0.2041241452319315)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*(-(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 1.0))))));
	sys.A( 0, vi + V2i(0,0), 2 ) += -(std::complex<double>(0.0, 0.2041241452319315)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*(-(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.0))))+
			(10.0*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 1.0))))));

	// term 23 ----------------
	sys.A( 0, vi + V2i(0,-1), 2 ) += (std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(0,0), 2 ) += (std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(1,-1), 2 ) += (std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(1,0), 2 ) += (std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(0,0), 2 ) += -(std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(0,1), 2 ) += -(std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(1,0), 2 ) += -(std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));
	sys.A( 0, vi + V2i(1,1), 2 ) += -(std::complex<double>(0.0, 1.0206207261596576)*fields.sigma_s->eval(sys.voxelToWorld(vd+V2d(0.5, 0.5)))*fields.f_p->eval(1, 0, sys.voxelToWorld(vd+V2d(0.5, 0.5))));

	// term 24 ----------------

	// term 25 ----------------

	// term 26 ----------------

	// term 27 ----------------

	// term 28 ----------------

	// term 29 ----------------
}
