// This file was generated by stencil2.py on Sunday, 11. February 2018 04:10PM

#include <PNSystem.h>

void stencil_p5_2d_cg(PNSystem::Stencil::Context& ctx)
{
	V3i vi = ctx.getVoxelCoord();
	V3d vd = vi.cast<double>();
	const Domain& domain = ctx.getDomain();
	const PNVolume& problem = ctx.getProblem();
	V3d h_inv( 1.0/(2*domain.getVoxelSize()[0]), 1.0/(2*domain.getVoxelSize()[1]), 1.0/(2*domain.getVoxelSize()[2]) );
	int color_channel = 0;

	// row=0 l=0 m=0 --------------------------
	ctx.coeff_A( 0, vi+V3i(-1, 0, 0), 2 ) += (0.57735026919*h_inv[0]);
	ctx.coeff_A( 0, vi+V3i(1, 0, 0), 2 ) += -(0.57735026919*h_inv[0]);
	ctx.coeff_A( 0, vi+V3i(0, -1, 0), 1 ) += (0.57735026919*h_inv[1]);
	ctx.coeff_A( 0, vi+V3i(0, 1, 0), 1 ) += -(0.57735026919*h_inv[1]);
	{
		double c = 0.0;
		c+=ctx.evalExtinction(0, 0, 0)[color_channel];
		c+=-(3.54490770181*ctx.evalScattering(0, 0, 0)[color_channel]*ctx.evalPhase(0, 0, 0, 0, 0)[color_channel]);
		ctx.coeff_A( 0, vi+V3i(0, 0, 0), 0 ) += c;
	}
	{
		double c = 0.0;
		c+=ctx.evalEmission(0, 0, 0, 0, 0)[color_channel];
		ctx.coeff_b( 0 ) += c;
	}
	// row=1 l=1 m=-1 --------------------------
	ctx.coeff_A( 1, vi+V3i(0, -1, 0), 5 ) += -(0.4472135955*h_inv[1]);
	ctx.coeff_A( 1, vi+V3i(0, 1, 0), 5 ) += (0.4472135955*h_inv[1]);
	ctx.coeff_A( 1, vi+V3i(0, -1, 0), 0 ) += (0.57735026919*h_inv[1]);
	ctx.coeff_A( 1, vi+V3i(0, 1, 0), 0 ) += -(0.57735026919*h_inv[1]);
	ctx.coeff_A( 1, vi+V3i(0, -1, 0), 4 ) += -(0.258198889747*h_inv[1]);
	ctx.coeff_A( 1, vi+V3i(0, 1, 0), 4 ) += (0.258198889747*h_inv[1]);
	ctx.coeff_A( 1, vi+V3i(-1, 0, 0), 3 ) += (0.4472135955*h_inv[0]);
	ctx.coeff_A( 1, vi+V3i(1, 0, 0), 3 ) += -(0.4472135955*h_inv[0]);
	{
		double c = 0.0;
		c+=ctx.evalExtinction(0, 0, 0)[color_channel];
		c+=-(2.04665341589*ctx.evalScattering(0, 0, 0)[color_channel]*ctx.evalPhase(1, 0, 0, 0, 0)[color_channel]);
		ctx.coeff_A( 1, vi+V3i(0, 0, 0), 1 ) += c;
	}
	{
		double c = 0.0;
		c+=ctx.evalEmission(1, -1, 0, 0, 0)[color_channel];
		ctx.coeff_b( 1 ) += c;
	}
	// row=2 l=1 m=1 --------------------------
	ctx.coeff_A( 2, vi+V3i(-1, 0, 0), 5 ) += (0.4472135955*h_inv[0]);
	ctx.coeff_A( 2, vi+V3i(1, 0, 0), 5 ) += -(0.4472135955*h_inv[0]);
	ctx.coeff_A( 2, vi+V3i(-1, 0, 0), 0 ) += (0.57735026919*h_inv[0]);
	ctx.coeff_A( 2, vi+V3i(1, 0, 0), 0 ) += -(0.57735026919*h_inv[0]);
	ctx.coeff_A( 2, vi+V3i(-1, 0, 0), 4 ) += -(0.258198889747*h_inv[0]);
	ctx.coeff_A( 2, vi+V3i(1, 0, 0), 4 ) += (0.258198889747*h_inv[0]);
	ctx.coeff_A( 2, vi+V3i(0, -1, 0), 3 ) += (0.4472135955*h_inv[1]);
	ctx.coeff_A( 2, vi+V3i(0, 1, 0), 3 ) += -(0.4472135955*h_inv[1]);
	{
		double c = 0.0;
		c+=ctx.evalExtinction(0, 0, 0)[color_channel];
		c+=-(2.04665341589*ctx.evalScattering(0, 0, 0)[color_channel]*ctx.evalPhase(1, 0, 0, 0, 0)[color_channel]);
		ctx.coeff_A( 2, vi+V3i(0, 0, 0), 2 ) += c;
	}
	{
		double c = 0.0;
		c+=ctx.evalEmission(1, 1, 0, 0, 0)[color_channel];
		ctx.coeff_b( 2 ) += c;
	}
	// row=3 l=2 m=-2 --------------------------
	ctx.coeff_A( 3, vi+V3i(0, -1, 0), 9 ) += -(0.462910049886*h_inv[1]);
	ctx.coeff_A( 3, vi+V3i(0, 1, 0), 9 ) += (0.462910049886*h_inv[1]);
	ctx.coeff_A( 3, vi+V3i(0, -1, 0), 2 ) += (0.4472135955*h_inv[1]);
	ctx.coeff_A( 3, vi+V3i(0, 1, 0), 2 ) += -(0.4472135955*h_inv[1]);
	ctx.coeff_A( 3, vi+V3i(0, -1, 0), 8 ) += -(0.119522860933*h_inv[1]);
	ctx.coeff_A( 3, vi+V3i(0, 1, 0), 8 ) += (0.119522860933*h_inv[1]);
	ctx.coeff_A( 3, vi+V3i(-1, 0, 0), 1 ) += (0.4472135955*h_inv[0]);
	ctx.coeff_A( 3, vi+V3i(1, 0, 0), 1 ) += -(0.4472135955*h_inv[0]);
	ctx.coeff_A( 3, vi+V3i(-1, 0, 0), 7 ) += -(0.119522860933*h_inv[0]);
	ctx.coeff_A( 3, vi+V3i(1, 0, 0), 7 ) += (0.119522860933*h_inv[0]);
	ctx.coeff_A( 3, vi+V3i(-1, 0, 0), 6 ) += (0.462910049886*h_inv[0]);
	ctx.coeff_A( 3, vi+V3i(1, 0, 0), 6 ) += -(0.462910049886*h_inv[0]);
	{
		double c = 0.0;
		c+=ctx.evalExtinction(0, 0, 0)[color_channel];
		c+=-(1.58533091904*ctx.evalScattering(0, 0, 0)[color_channel]*ctx.evalPhase(2, 0, 0, 0, 0)[color_channel]);
		ctx.coeff_A( 3, vi+V3i(0, 0, 0), 3 ) += c;
	}
	{
		double c = 0.0;
		c+=ctx.evalEmission(2, -2, 0, 0, 0)[color_channel];
		ctx.coeff_b( 3 ) += c;
	}
	// row=4 l=2 m=0 --------------------------
	ctx.coeff_A( 4, vi+V3i(-1, 0, 0), 2 ) += -(0.258198889747*h_inv[0]);
	ctx.coeff_A( 4, vi+V3i(1, 0, 0), 2 ) += (0.258198889747*h_inv[0]);
	ctx.coeff_A( 4, vi+V3i(-1, 0, 0), 8 ) += (0.414039335605*h_inv[0]);
	ctx.coeff_A( 4, vi+V3i(1, 0, 0), 8 ) += -(0.414039335605*h_inv[0]);
	ctx.coeff_A( 4, vi+V3i(0, -1, 0), 1 ) += -(0.258198889747*h_inv[1]);
	ctx.coeff_A( 4, vi+V3i(0, 1, 0), 1 ) += (0.258198889747*h_inv[1]);
	ctx.coeff_A( 4, vi+V3i(0, -1, 0), 7 ) += (0.414039335605*h_inv[1]);
	ctx.coeff_A( 4, vi+V3i(0, 1, 0), 7 ) += -(0.414039335605*h_inv[1]);
	{
		double c = 0.0;
		c+=ctx.evalExtinction(0, 0, 0)[color_channel];
		c+=-(1.58533091904*ctx.evalScattering(0, 0, 0)[color_channel]*ctx.evalPhase(2, 0, 0, 0, 0)[color_channel]);
		ctx.coeff_A( 4, vi+V3i(0, 0, 0), 4 ) += c;
	}
	{
		double c = 0.0;
		c+=ctx.evalEmission(2, 0, 0, 0, 0)[color_channel];
		ctx.coeff_b( 4 ) += c;
	}
	// row=5 l=2 m=2 --------------------------
	ctx.coeff_A( 5, vi+V3i(-1, 0, 0), 9 ) += (0.462910049886*h_inv[0]);
	ctx.coeff_A( 5, vi+V3i(1, 0, 0), 9 ) += -(0.462910049886*h_inv[0]);
	ctx.coeff_A( 5, vi+V3i(-1, 0, 0), 2 ) += (0.4472135955*h_inv[0]);
	ctx.coeff_A( 5, vi+V3i(1, 0, 0), 2 ) += -(0.4472135955*h_inv[0]);
	ctx.coeff_A( 5, vi+V3i(-1, 0, 0), 8 ) += -(0.119522860933*h_inv[0]);
	ctx.coeff_A( 5, vi+V3i(1, 0, 0), 8 ) += (0.119522860933*h_inv[0]);
	ctx.coeff_A( 5, vi+V3i(0, -1, 0), 6 ) += (0.462910049886*h_inv[1]);
	ctx.coeff_A( 5, vi+V3i(0, 1, 0), 6 ) += -(0.462910049886*h_inv[1]);
	ctx.coeff_A( 5, vi+V3i(0, -1, 0), 1 ) += -(0.4472135955*h_inv[1]);
	ctx.coeff_A( 5, vi+V3i(0, 1, 0), 1 ) += (0.4472135955*h_inv[1]);
	ctx.coeff_A( 5, vi+V3i(0, -1, 0), 7 ) += (0.119522860933*h_inv[1]);
	ctx.coeff_A( 5, vi+V3i(0, 1, 0), 7 ) += -(0.119522860933*h_inv[1]);
	{
		double c = 0.0;
		c+=ctx.evalExtinction(0, 0, 0)[color_channel];
		c+=-(1.58533091904*ctx.evalScattering(0, 0, 0)[color_channel]*ctx.evalPhase(2, 0, 0, 0, 0)[color_channel]);
		ctx.coeff_A( 5, vi+V3i(0, 0, 0), 5 ) += c;
	}
	{
		double c = 0.0;
		c+=ctx.evalEmission(2, 2, 0, 0, 0)[color_channel];
		ctx.coeff_b( 5 ) += c;
	}
	// row=6 l=3 m=-3 --------------------------
	ctx.coeff_A( 6, vi+V3i(0, -1, 0), 14 ) += -(0.471404520791*h_inv[1]);
	ctx.coeff_A( 6, vi+V3i(0, 1, 0), 14 ) += (0.471404520791*h_inv[1]);
	ctx.coeff_A( 6, vi+V3i(0, -1, 0), 5 ) += (0.462910049886*h_inv[1]);
	ctx.coeff_A( 6, vi+V3i(0, 1, 0), 5 ) += -(0.462910049886*h_inv[1]);
	ctx.coeff_A( 6, vi+V3i(0, -1, 0), 13 ) += -(0.0890870806375*h_inv[1]);
	ctx.coeff_A( 6, vi+V3i(0, 1, 0), 13 ) += (0.0890870806375*h_inv[1]);
	ctx.coeff_A( 6, vi+V3i(-1, 0, 0), 3 ) += (0.462910049886*h_inv[0]);
	ctx.coeff_A( 6, vi+V3i(1, 0, 0), 3 ) += -(0.462910049886*h_inv[0]);
	ctx.coeff_A( 6, vi+V3i(-1, 0, 0), 11 ) += -(0.0890870806375*h_inv[0]);
	ctx.coeff_A( 6, vi+V3i(1, 0, 0), 11 ) += (0.0890870806375*h_inv[0]);
	ctx.coeff_A( 6, vi+V3i(-1, 0, 0), 10 ) += (0.471404520791*h_inv[0]);
	ctx.coeff_A( 6, vi+V3i(1, 0, 0), 10 ) += -(0.471404520791*h_inv[0]);
	{
		double c = 0.0;
		c+=ctx.evalExtinction(0, 0, 0)[color_channel];
		c+=-(1.33984917138*ctx.evalScattering(0, 0, 0)[color_channel]*ctx.evalPhase(3, 0, 0, 0, 0)[color_channel]);
		ctx.coeff_A( 6, vi+V3i(0, 0, 0), 6 ) += c;
	}
	{
		double c = 0.0;
		c+=ctx.evalEmission(3, -3, 0, 0, 0)[color_channel];
		ctx.coeff_b( 6 ) += c;
	}
	// row=7 l=3 m=-1 --------------------------
	ctx.coeff_A( 7, vi+V3i(0, -1, 0), 5 ) += (0.119522860933*h_inv[1]);
	ctx.coeff_A( 7, vi+V3i(0, 1, 0), 5 ) += -(0.119522860933*h_inv[1]);
	ctx.coeff_A( 7, vi+V3i(0, -1, 0), 13 ) += -(0.345032779671*h_inv[1]);
	ctx.coeff_A( 7, vi+V3i(0, 1, 0), 13 ) += (0.345032779671*h_inv[1]);
	ctx.coeff_A( 7, vi+V3i(0, -1, 0), 4 ) += (0.414039335605*h_inv[1]);
	ctx.coeff_A( 7, vi+V3i(0, 1, 0), 4 ) += -(0.414039335605*h_inv[1]);
	ctx.coeff_A( 7, vi+V3i(0, -1, 0), 12 ) += -(0.308606699924*h_inv[1]);
	ctx.coeff_A( 7, vi+V3i(0, 1, 0), 12 ) += (0.308606699924*h_inv[1]);
	ctx.coeff_A( 7, vi+V3i(-1, 0, 0), 3 ) += -(0.119522860933*h_inv[0]);
	ctx.coeff_A( 7, vi+V3i(1, 0, 0), 3 ) += (0.119522860933*h_inv[0]);
	ctx.coeff_A( 7, vi+V3i(-1, 0, 0), 11 ) += (0.345032779671*h_inv[0]);
	ctx.coeff_A( 7, vi+V3i(1, 0, 0), 11 ) += -(0.345032779671*h_inv[0]);
	{
		double c = 0.0;
		c+=ctx.evalExtinction(0, 0, 0)[color_channel];
		c+=-(1.33984917138*ctx.evalScattering(0, 0, 0)[color_channel]*ctx.evalPhase(3, 0, 0, 0, 0)[color_channel]);
		ctx.coeff_A( 7, vi+V3i(0, 0, 0), 7 ) += c;
	}
	{
		double c = 0.0;
		c+=ctx.evalEmission(3, -1, 0, 0, 0)[color_channel];
		ctx.coeff_b( 7 ) += c;
	}
	// row=8 l=3 m=1 --------------------------
	ctx.coeff_A( 8, vi+V3i(-1, 0, 0), 5 ) += -(0.119522860933*h_inv[0]);
	ctx.coeff_A( 8, vi+V3i(1, 0, 0), 5 ) += (0.119522860933*h_inv[0]);
	ctx.coeff_A( 8, vi+V3i(-1, 0, 0), 13 ) += (0.345032779671*h_inv[0]);
	ctx.coeff_A( 8, vi+V3i(1, 0, 0), 13 ) += -(0.345032779671*h_inv[0]);
	ctx.coeff_A( 8, vi+V3i(-1, 0, 0), 4 ) += (0.414039335605*h_inv[0]);
	ctx.coeff_A( 8, vi+V3i(1, 0, 0), 4 ) += -(0.414039335605*h_inv[0]);
	ctx.coeff_A( 8, vi+V3i(-1, 0, 0), 12 ) += -(0.308606699924*h_inv[0]);
	ctx.coeff_A( 8, vi+V3i(1, 0, 0), 12 ) += (0.308606699924*h_inv[0]);
	ctx.coeff_A( 8, vi+V3i(0, -1, 0), 3 ) += -(0.119522860933*h_inv[1]);
	ctx.coeff_A( 8, vi+V3i(0, 1, 0), 3 ) += (0.119522860933*h_inv[1]);
	ctx.coeff_A( 8, vi+V3i(0, -1, 0), 11 ) += (0.345032779671*h_inv[1]);
	ctx.coeff_A( 8, vi+V3i(0, 1, 0), 11 ) += -(0.345032779671*h_inv[1]);
	{
		double c = 0.0;
		c+=ctx.evalExtinction(0, 0, 0)[color_channel];
		c+=-(1.33984917138*ctx.evalScattering(0, 0, 0)[color_channel]*ctx.evalPhase(3, 0, 0, 0, 0)[color_channel]);
		ctx.coeff_A( 8, vi+V3i(0, 0, 0), 8 ) += c;
	}
	{
		double c = 0.0;
		c+=ctx.evalEmission(3, 1, 0, 0, 0)[color_channel];
		ctx.coeff_b( 8 ) += c;
	}
	// row=9 l=3 m=3 --------------------------
	ctx.coeff_A( 9, vi+V3i(-1, 0, 0), 14 ) += (0.471404520791*h_inv[0]);
	ctx.coeff_A( 9, vi+V3i(1, 0, 0), 14 ) += -(0.471404520791*h_inv[0]);
	ctx.coeff_A( 9, vi+V3i(-1, 0, 0), 5 ) += (0.462910049886*h_inv[0]);
	ctx.coeff_A( 9, vi+V3i(1, 0, 0), 5 ) += -(0.462910049886*h_inv[0]);
	ctx.coeff_A( 9, vi+V3i(-1, 0, 0), 13 ) += -(0.0890870806375*h_inv[0]);
	ctx.coeff_A( 9, vi+V3i(1, 0, 0), 13 ) += (0.0890870806375*h_inv[0]);
	ctx.coeff_A( 9, vi+V3i(0, -1, 0), 10 ) += (0.471404520791*h_inv[1]);
	ctx.coeff_A( 9, vi+V3i(0, 1, 0), 10 ) += -(0.471404520791*h_inv[1]);
	ctx.coeff_A( 9, vi+V3i(0, -1, 0), 3 ) += -(0.462910049886*h_inv[1]);
	ctx.coeff_A( 9, vi+V3i(0, 1, 0), 3 ) += (0.462910049886*h_inv[1]);
	ctx.coeff_A( 9, vi+V3i(0, -1, 0), 11 ) += (0.0890870806375*h_inv[1]);
	ctx.coeff_A( 9, vi+V3i(0, 1, 0), 11 ) += -(0.0890870806375*h_inv[1]);
	{
		double c = 0.0;
		c+=ctx.evalExtinction(0, 0, 0)[color_channel];
		c+=-(1.33984917138*ctx.evalScattering(0, 0, 0)[color_channel]*ctx.evalPhase(3, 0, 0, 0, 0)[color_channel]);
		ctx.coeff_A( 9, vi+V3i(0, 0, 0), 9 ) += c;
	}
	{
		double c = 0.0;
		c+=ctx.evalEmission(3, 3, 0, 0, 0)[color_channel];
		ctx.coeff_b( 9 ) += c;
	}
	// row=10 l=4 m=-4 --------------------------
	ctx.coeff_A( 10, vi+V3i(0, -1, 0), 20 ) += -(0.476731294623*h_inv[1]);
	ctx.coeff_A( 10, vi+V3i(0, 1, 0), 20 ) += (0.476731294623*h_inv[1]);
	ctx.coeff_A( 10, vi+V3i(0, -1, 0), 9 ) += (0.471404520791*h_inv[1]);
	ctx.coeff_A( 10, vi+V3i(0, 1, 0), 9 ) += -(0.471404520791*h_inv[1]);
	ctx.coeff_A( 10, vi+V3i(0, -1, 0), 19 ) += -(0.0710669054519*h_inv[1]);
	ctx.coeff_A( 10, vi+V3i(0, 1, 0), 19 ) += (0.0710669054519*h_inv[1]);
	ctx.coeff_A( 10, vi+V3i(-1, 0, 0), 6 ) += (0.471404520791*h_inv[0]);
	ctx.coeff_A( 10, vi+V3i(1, 0, 0), 6 ) += -(0.471404520791*h_inv[0]);
	ctx.coeff_A( 10, vi+V3i(-1, 0, 0), 16 ) += -(0.0710669054519*h_inv[0]);
	ctx.coeff_A( 10, vi+V3i(1, 0, 0), 16 ) += (0.0710669054519*h_inv[0]);
	ctx.coeff_A( 10, vi+V3i(-1, 0, 0), 15 ) += (0.476731294623*h_inv[0]);
	ctx.coeff_A( 10, vi+V3i(1, 0, 0), 15 ) += -(0.476731294623*h_inv[0]);
	{
		double c = 0.0;
		c+=ctx.evalExtinction(0, 0, 0)[color_channel];
		c+=-(1.1816359006*ctx.evalScattering(0, 0, 0)[color_channel]*ctx.evalPhase(4, 0, 0, 0, 0)[color_channel]);
		ctx.coeff_A( 10, vi+V3i(0, 0, 0), 10 ) += c;
	}
	{
		double c = 0.0;
		c+=ctx.evalEmission(4, -4, 0, 0, 0)[color_channel];
		ctx.coeff_b( 10 ) += c;
	}
	// row=11 l=4 m=-2 --------------------------
	ctx.coeff_A( 11, vi+V3i(0, -1, 0), 9 ) += (0.0890870806375*h_inv[1]);
	ctx.coeff_A( 11, vi+V3i(0, 1, 0), 9 ) += -(0.0890870806375*h_inv[1]);
	ctx.coeff_A( 11, vi+V3i(0, -1, 0), 19 ) += -(0.376050716545*h_inv[1]);
	ctx.coeff_A( 11, vi+V3i(0, 1, 0), 19 ) += (0.376050716545*h_inv[1]);
	ctx.coeff_A( 11, vi+V3i(0, -1, 0), 8 ) += (0.345032779671*h_inv[1]);
	ctx.coeff_A( 11, vi+V3i(0, 1, 0), 8 ) += -(0.345032779671*h_inv[1]);
	ctx.coeff_A( 11, vi+V3i(0, -1, 0), 18 ) += -(0.174077655956*h_inv[1]);
	ctx.coeff_A( 11, vi+V3i(0, 1, 0), 18 ) += (0.174077655956*h_inv[1]);
	ctx.coeff_A( 11, vi+V3i(-1, 0, 0), 6 ) += -(0.0890870806375*h_inv[0]);
	ctx.coeff_A( 11, vi+V3i(1, 0, 0), 6 ) += (0.0890870806375*h_inv[0]);
	ctx.coeff_A( 11, vi+V3i(-1, 0, 0), 7 ) += (0.345032779671*h_inv[0]);
	ctx.coeff_A( 11, vi+V3i(1, 0, 0), 7 ) += -(0.345032779671*h_inv[0]);
	ctx.coeff_A( 11, vi+V3i(-1, 0, 0), 17 ) += -(0.174077655956*h_inv[0]);
	ctx.coeff_A( 11, vi+V3i(1, 0, 0), 17 ) += (0.174077655956*h_inv[0]);
	ctx.coeff_A( 11, vi+V3i(-1, 0, 0), 16 ) += (0.376050716545*h_inv[0]);
	ctx.coeff_A( 11, vi+V3i(1, 0, 0), 16 ) += -(0.376050716545*h_inv[0]);
	{
		double c = 0.0;
		c+=ctx.evalExtinction(0, 0, 0)[color_channel];
		c+=-(1.1816359006*ctx.evalScattering(0, 0, 0)[color_channel]*ctx.evalPhase(4, 0, 0, 0, 0)[color_channel]);
		ctx.coeff_A( 11, vi+V3i(0, 0, 0), 11 ) += c;
	}
	{
		double c = 0.0;
		c+=ctx.evalEmission(4, -2, 0, 0, 0)[color_channel];
		ctx.coeff_b( 11 ) += c;
	}
	// row=12 l=4 m=0 --------------------------
	ctx.coeff_A( 12, vi+V3i(-1, 0, 0), 8 ) += -(0.308606699924*h_inv[0]);
	ctx.coeff_A( 12, vi+V3i(1, 0, 0), 8 ) += (0.308606699924*h_inv[0]);
	ctx.coeff_A( 12, vi+V3i(-1, 0, 0), 18 ) += (0.389249472081*h_inv[0]);
	ctx.coeff_A( 12, vi+V3i(1, 0, 0), 18 ) += -(0.389249472081*h_inv[0]);
	ctx.coeff_A( 12, vi+V3i(0, -1, 0), 7 ) += -(0.308606699924*h_inv[1]);
	ctx.coeff_A( 12, vi+V3i(0, 1, 0), 7 ) += (0.308606699924*h_inv[1]);
	ctx.coeff_A( 12, vi+V3i(0, -1, 0), 17 ) += (0.389249472081*h_inv[1]);
	ctx.coeff_A( 12, vi+V3i(0, 1, 0), 17 ) += -(0.389249472081*h_inv[1]);
	{
		double c = 0.0;
		c+=ctx.evalExtinction(0, 0, 0)[color_channel];
		c+=-(1.1816359006*ctx.evalScattering(0, 0, 0)[color_channel]*ctx.evalPhase(4, 0, 0, 0, 0)[color_channel]);
		ctx.coeff_A( 12, vi+V3i(0, 0, 0), 12 ) += c;
	}
	{
		double c = 0.0;
		c+=ctx.evalEmission(4, 0, 0, 0, 0)[color_channel];
		ctx.coeff_b( 12 ) += c;
	}
	// row=13 l=4 m=2 --------------------------
	ctx.coeff_A( 13, vi+V3i(-1, 0, 0), 9 ) += -(0.0890870806375*h_inv[0]);
	ctx.coeff_A( 13, vi+V3i(1, 0, 0), 9 ) += (0.0890870806375*h_inv[0]);
	ctx.coeff_A( 13, vi+V3i(-1, 0, 0), 19 ) += (0.376050716545*h_inv[0]);
	ctx.coeff_A( 13, vi+V3i(1, 0, 0), 19 ) += -(0.376050716545*h_inv[0]);
	ctx.coeff_A( 13, vi+V3i(-1, 0, 0), 8 ) += (0.345032779671*h_inv[0]);
	ctx.coeff_A( 13, vi+V3i(1, 0, 0), 8 ) += -(0.345032779671*h_inv[0]);
	ctx.coeff_A( 13, vi+V3i(-1, 0, 0), 18 ) += -(0.174077655956*h_inv[0]);
	ctx.coeff_A( 13, vi+V3i(1, 0, 0), 18 ) += (0.174077655956*h_inv[0]);
	ctx.coeff_A( 13, vi+V3i(0, -1, 0), 6 ) += -(0.0890870806375*h_inv[1]);
	ctx.coeff_A( 13, vi+V3i(0, 1, 0), 6 ) += (0.0890870806375*h_inv[1]);
	ctx.coeff_A( 13, vi+V3i(0, -1, 0), 16 ) += (0.376050716545*h_inv[1]);
	ctx.coeff_A( 13, vi+V3i(0, 1, 0), 16 ) += -(0.376050716545*h_inv[1]);
	ctx.coeff_A( 13, vi+V3i(0, -1, 0), 7 ) += -(0.345032779671*h_inv[1]);
	ctx.coeff_A( 13, vi+V3i(0, 1, 0), 7 ) += (0.345032779671*h_inv[1]);
	ctx.coeff_A( 13, vi+V3i(0, -1, 0), 17 ) += (0.174077655956*h_inv[1]);
	ctx.coeff_A( 13, vi+V3i(0, 1, 0), 17 ) += -(0.174077655956*h_inv[1]);
	{
		double c = 0.0;
		c+=ctx.evalExtinction(0, 0, 0)[color_channel];
		c+=-(1.1816359006*ctx.evalScattering(0, 0, 0)[color_channel]*ctx.evalPhase(4, 0, 0, 0, 0)[color_channel]);
		ctx.coeff_A( 13, vi+V3i(0, 0, 0), 13 ) += c;
	}
	{
		double c = 0.0;
		c+=ctx.evalEmission(4, 2, 0, 0, 0)[color_channel];
		ctx.coeff_b( 13 ) += c;
	}
	// row=14 l=4 m=4 --------------------------
	ctx.coeff_A( 14, vi+V3i(-1, 0, 0), 20 ) += (0.476731294623*h_inv[0]);
	ctx.coeff_A( 14, vi+V3i(1, 0, 0), 20 ) += -(0.476731294623*h_inv[0]);
	ctx.coeff_A( 14, vi+V3i(-1, 0, 0), 9 ) += (0.471404520791*h_inv[0]);
	ctx.coeff_A( 14, vi+V3i(1, 0, 0), 9 ) += -(0.471404520791*h_inv[0]);
	ctx.coeff_A( 14, vi+V3i(-1, 0, 0), 19 ) += -(0.0710669054519*h_inv[0]);
	ctx.coeff_A( 14, vi+V3i(1, 0, 0), 19 ) += (0.0710669054519*h_inv[0]);
	ctx.coeff_A( 14, vi+V3i(0, -1, 0), 15 ) += (0.476731294623*h_inv[1]);
	ctx.coeff_A( 14, vi+V3i(0, 1, 0), 15 ) += -(0.476731294623*h_inv[1]);
	ctx.coeff_A( 14, vi+V3i(0, -1, 0), 6 ) += -(0.471404520791*h_inv[1]);
	ctx.coeff_A( 14, vi+V3i(0, 1, 0), 6 ) += (0.471404520791*h_inv[1]);
	ctx.coeff_A( 14, vi+V3i(0, -1, 0), 16 ) += (0.0710669054519*h_inv[1]);
	ctx.coeff_A( 14, vi+V3i(0, 1, 0), 16 ) += -(0.0710669054519*h_inv[1]);
	{
		double c = 0.0;
		c+=ctx.evalExtinction(0, 0, 0)[color_channel];
		c+=-(1.1816359006*ctx.evalScattering(0, 0, 0)[color_channel]*ctx.evalPhase(4, 0, 0, 0, 0)[color_channel]);
		ctx.coeff_A( 14, vi+V3i(0, 0, 0), 14 ) += c;
	}
	{
		double c = 0.0;
		c+=ctx.evalEmission(4, 4, 0, 0, 0)[color_channel];
		ctx.coeff_b( 14 ) += c;
	}
	// row=15 l=5 m=-5 --------------------------
	ctx.coeff_A( 15, vi+V3i(0, -1, 0), 14 ) += (0.476731294623*h_inv[1]);
	ctx.coeff_A( 15, vi+V3i(0, 1, 0), 14 ) += -(0.476731294623*h_inv[1]);
	ctx.coeff_A( 15, vi+V3i(-1, 0, 0), 10 ) += (0.476731294623*h_inv[0]);
	ctx.coeff_A( 15, vi+V3i(1, 0, 0), 10 ) += -(0.476731294623*h_inv[0]);
	{
		double c = 0.0;
		c+=ctx.evalExtinction(0, 0, 0)[color_channel];
		c+=-(1.06882988758*ctx.evalScattering(0, 0, 0)[color_channel]*ctx.evalPhase(5, 0, 0, 0, 0)[color_channel]);
		ctx.coeff_A( 15, vi+V3i(0, 0, 0), 15 ) += c;
	}
	{
		double c = 0.0;
		c+=ctx.evalEmission(5, -5, 0, 0, 0)[color_channel];
		ctx.coeff_b( 15 ) += c;
	}
	// row=16 l=5 m=-3 --------------------------
	ctx.coeff_A( 16, vi+V3i(0, -1, 0), 14 ) += (0.0710669054519*h_inv[1]);
	ctx.coeff_A( 16, vi+V3i(0, 1, 0), 14 ) += -(0.0710669054519*h_inv[1]);
	ctx.coeff_A( 16, vi+V3i(0, -1, 0), 13 ) += (0.376050716545*h_inv[1]);
	ctx.coeff_A( 16, vi+V3i(0, 1, 0), 13 ) += -(0.376050716545*h_inv[1]);
	ctx.coeff_A( 16, vi+V3i(-1, 0, 0), 10 ) += -(0.0710669054519*h_inv[0]);
	ctx.coeff_A( 16, vi+V3i(1, 0, 0), 10 ) += (0.0710669054519*h_inv[0]);
	ctx.coeff_A( 16, vi+V3i(-1, 0, 0), 11 ) += (0.376050716545*h_inv[0]);
	ctx.coeff_A( 16, vi+V3i(1, 0, 0), 11 ) += -(0.376050716545*h_inv[0]);
	{
		double c = 0.0;
		c+=ctx.evalExtinction(0, 0, 0)[color_channel];
		c+=-(1.06882988758*ctx.evalScattering(0, 0, 0)[color_channel]*ctx.evalPhase(5, 0, 0, 0, 0)[color_channel]);
		ctx.coeff_A( 16, vi+V3i(0, 0, 0), 16 ) += c;
	}
	{
		double c = 0.0;
		c+=ctx.evalEmission(5, -3, 0, 0, 0)[color_channel];
		ctx.coeff_b( 16 ) += c;
	}
	// row=17 l=5 m=-1 --------------------------
	ctx.coeff_A( 17, vi+V3i(0, -1, 0), 13 ) += (0.174077655956*h_inv[1]);
	ctx.coeff_A( 17, vi+V3i(0, 1, 0), 13 ) += -(0.174077655956*h_inv[1]);
	ctx.coeff_A( 17, vi+V3i(0, -1, 0), 12 ) += (0.389249472081*h_inv[1]);
	ctx.coeff_A( 17, vi+V3i(0, 1, 0), 12 ) += -(0.389249472081*h_inv[1]);
	ctx.coeff_A( 17, vi+V3i(-1, 0, 0), 11 ) += -(0.174077655956*h_inv[0]);
	ctx.coeff_A( 17, vi+V3i(1, 0, 0), 11 ) += (0.174077655956*h_inv[0]);
	{
		double c = 0.0;
		c+=ctx.evalExtinction(0, 0, 0)[color_channel];
		c+=-(1.06882988758*ctx.evalScattering(0, 0, 0)[color_channel]*ctx.evalPhase(5, 0, 0, 0, 0)[color_channel]);
		ctx.coeff_A( 17, vi+V3i(0, 0, 0), 17 ) += c;
	}
	{
		double c = 0.0;
		c+=ctx.evalEmission(5, -1, 0, 0, 0)[color_channel];
		ctx.coeff_b( 17 ) += c;
	}
	// row=18 l=5 m=1 --------------------------
	ctx.coeff_A( 18, vi+V3i(-1, 0, 0), 13 ) += -(0.174077655956*h_inv[0]);
	ctx.coeff_A( 18, vi+V3i(1, 0, 0), 13 ) += (0.174077655956*h_inv[0]);
	ctx.coeff_A( 18, vi+V3i(-1, 0, 0), 12 ) += (0.389249472081*h_inv[0]);
	ctx.coeff_A( 18, vi+V3i(1, 0, 0), 12 ) += -(0.389249472081*h_inv[0]);
	ctx.coeff_A( 18, vi+V3i(0, -1, 0), 11 ) += -(0.174077655956*h_inv[1]);
	ctx.coeff_A( 18, vi+V3i(0, 1, 0), 11 ) += (0.174077655956*h_inv[1]);
	{
		double c = 0.0;
		c+=ctx.evalExtinction(0, 0, 0)[color_channel];
		c+=-(1.06882988758*ctx.evalScattering(0, 0, 0)[color_channel]*ctx.evalPhase(5, 0, 0, 0, 0)[color_channel]);
		ctx.coeff_A( 18, vi+V3i(0, 0, 0), 18 ) += c;
	}
	{
		double c = 0.0;
		c+=ctx.evalEmission(5, 1, 0, 0, 0)[color_channel];
		ctx.coeff_b( 18 ) += c;
	}
	// row=19 l=5 m=3 --------------------------
	ctx.coeff_A( 19, vi+V3i(-1, 0, 0), 14 ) += -(0.0710669054519*h_inv[0]);
	ctx.coeff_A( 19, vi+V3i(1, 0, 0), 14 ) += (0.0710669054519*h_inv[0]);
	ctx.coeff_A( 19, vi+V3i(-1, 0, 0), 13 ) += (0.376050716545*h_inv[0]);
	ctx.coeff_A( 19, vi+V3i(1, 0, 0), 13 ) += -(0.376050716545*h_inv[0]);
	ctx.coeff_A( 19, vi+V3i(0, -1, 0), 10 ) += -(0.0710669054519*h_inv[1]);
	ctx.coeff_A( 19, vi+V3i(0, 1, 0), 10 ) += (0.0710669054519*h_inv[1]);
	ctx.coeff_A( 19, vi+V3i(0, -1, 0), 11 ) += -(0.376050716545*h_inv[1]);
	ctx.coeff_A( 19, vi+V3i(0, 1, 0), 11 ) += (0.376050716545*h_inv[1]);
	{
		double c = 0.0;
		c+=ctx.evalExtinction(0, 0, 0)[color_channel];
		c+=-(1.06882988758*ctx.evalScattering(0, 0, 0)[color_channel]*ctx.evalPhase(5, 0, 0, 0, 0)[color_channel]);
		ctx.coeff_A( 19, vi+V3i(0, 0, 0), 19 ) += c;
	}
	{
		double c = 0.0;
		c+=ctx.evalEmission(5, 3, 0, 0, 0)[color_channel];
		ctx.coeff_b( 19 ) += c;
	}
	// row=20 l=5 m=5 --------------------------
	ctx.coeff_A( 20, vi+V3i(-1, 0, 0), 14 ) += (0.476731294623*h_inv[0]);
	ctx.coeff_A( 20, vi+V3i(1, 0, 0), 14 ) += -(0.476731294623*h_inv[0]);
	ctx.coeff_A( 20, vi+V3i(0, -1, 0), 10 ) += -(0.476731294623*h_inv[1]);
	ctx.coeff_A( 20, vi+V3i(0, 1, 0), 10 ) += (0.476731294623*h_inv[1]);
	{
		double c = 0.0;
		c+=ctx.evalExtinction(0, 0, 0)[color_channel];
		c+=-(1.06882988758*ctx.evalScattering(0, 0, 0)[color_channel]*ctx.evalPhase(5, 0, 0, 0, 0)[color_channel]);
		ctx.coeff_A( 20, vi+V3i(0, 0, 0), 20 ) += c;
	}
	{
		double c = 0.0;
		c+=ctx.evalEmission(5, 5, 0, 0, 0)[color_channel];
		ctx.coeff_b( 20 ) += c;
	}
}
V3i stencil_p5_2d_cg_get_offset(int coeff)
{
	switch(coeff)
	{
		case 0:return V3i(1, 1, 1);break;
		case 1:return V3i(1, 1, 1);break;
		case 2:return V3i(1, 1, 1);break;
		case 3:return V3i(1, 1, 1);break;
		case 4:return V3i(1, 1, 1);break;
		case 5:return V3i(1, 1, 1);break;
		case 6:return V3i(1, 1, 1);break;
		case 7:return V3i(1, 1, 1);break;
		case 8:return V3i(1, 1, 1);break;
		case 9:return V3i(1, 1, 1);break;
		case 10:return V3i(1, 1, 1);break;
		case 11:return V3i(1, 1, 1);break;
		case 12:return V3i(1, 1, 1);break;
		case 13:return V3i(1, 1, 1);break;
		case 14:return V3i(1, 1, 1);break;
		case 15:return V3i(1, 1, 1);break;
		case 16:return V3i(1, 1, 1);break;
		case 17:return V3i(1, 1, 1);break;
		case 18:return V3i(1, 1, 1);break;
		case 19:return V3i(1, 1, 1);break;
		case 20:return V3i(1, 1, 1);break;
		default:throw std::runtime_error("unexpected coefficient index");break;
	};
}
REGISTER_STENCIL(stencil_p5_2d_cg, 5, 21, 1)
