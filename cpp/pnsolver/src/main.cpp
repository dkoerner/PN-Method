#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <complex>
#include <iostream>

#include <math/common.h>
#include <math/vector.h>
#include <math/ray.h>

#include <util/timer.h>
#include <util/mem.h>

#include<common/Domain.h>
#include<field/VoxelGridField.h>
#include<field/Constant.h>
#include<field/SHEXP.h>
#include <PNSystem.h>

#include <solver.h>


namespace py = pybind11;



void setup_solver( Solver& mg, PNSystem& sys, int numLevels = 1 )
{
	mg = Solver(numLevels);

	PNSystem::Stencil& stencil = sys.getStencil();
	int boundaryConditions = sys.getBoundaryConditions();
	Domain domain = sys.getDomain();
	PNSystem::Fields fields = sys.getFields();

	for( int i=0;i<numLevels;++i )
	{
		PNSystem::RealMatrix A;
		Eigen::VectorXd x;
		PNSystem::RealMatrix downsample;
		PNSystem::RealMatrix upsample;

		PNSystem sys_level(stencil, domain, boundaryConditions);
		sys_level.setFields( fields );
		sys_level.build();
		buildUpAndDownsamplingMatrices(sys_level, downsample, upsample);


		A = sys_level.get_A_real().transpose()*sys_level.get_A_real();
		Eigen::VectorXd b = sys_level.get_A_real().transpose()*sys_level.get_b_real();
		//std::cout << "warning: not using ATA form\n";
		//A = sys_level.get_A_real();
		//Eigen::VectorXd b = sys_level.get_b_real();

		x = Eigen::VectorXd(sys_level.getVoxelManager().getNumUnknowns());
		x.fill(0.0);

		mg.setMultigridLevel(i, A, x, b, downsample, upsample);

		// used for blockgs test
		mg.m_levels[i].m_voxelManager = sys_level.getVoxelManager();


		// downsample to next level
		domain = domain.downsample();
		fields = fields.createRestricted();
	}
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_multigrid(PNSystem& sys, int numLevels)
{
	Solver mg;
	setup_solver(mg, sys, numLevels);
	int maxIterations = 1000;
	auto result = mg.solve(maxIterations);

	//sys.debug_x = std::get<0>(result);
	//sys.debug_x_downsampled = sys.debug_downsample*sys.debug_x;
	//sys.debug_x_up_sampled_downsampled = sys.debug_upsample*sys.debug_x_downsampled;


	return std::make_tuple(sys.stripBoundary(std::get<0>(result)), std::get<1>(result), std::get<2>(result));
}


#include <cstdint>

constexpr int64_t ipow(int64_t base, int exp, int64_t result = 1)
{
	if(base==0)
		return 1;
	return exp < 1 ? result : ipow(base*base, exp/2, (exp % 2) ? result*base : result);
}


std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_multigrid2(PNSystem& sys, int numLevels)
{
	Solver mg;
	setup_solver(mg, sys, numLevels);

	/*
	for( int i=0;i<mg.m_levels.size();++i )
	{
		mg.m_levels[i].numPreSmoothingIterations = ipow(2,i);
		std::cout << "level=" << i <<" #iter=" << mg.m_levels[i].numPreSmoothingIterations << std::endl;
	}


	int maxIterations = 1000;
	auto result = mg.solve2(maxIterations);

	return std::make_tuple(sys.stripBoundary(std::get<0>(result)), std::get<1>(result), std::get<2>(result));
	*/

	/*
	Solver::System& lgs = mg.m_levels[0];

	Timer timer;
	std::vector<double> convergence;
	std::vector<double> convergence_time;

	convergence.push_back(lgs.computeRMSE());
	convergence_time.push_back(0.0);

	{
		Solver::System& lgs = mg.m_levels[1];
		mg.run_cg_iterations( lgs.A, lgs.x, lgs.b, lgs.r, 1000, 1.0e-10 );
		mg.m_levels[0].x = mg.m_levels[0].upsample*lgs.x;
	}

	timer.start();

	//std::vector<double> tolerances = {1.0e-10, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8, 1.0e-9, 1.0e-10};

//	for( int i=0;i<numLevels-1;++i )
//	{
//		int level = numLevels-1-i;
//		Solver::System& lgs = mg.m_levels[level];
//		double tol = 1.0e-10;
//		//double tol = tolerances[level];
//		mg.run_cg_iterations( lgs.A, lgs.x, lgs.b, lgs.r, 1000, tol );
//		mg.m_levels[level-1].x = mg.m_levels[level-1].upsample*lgs.x;
//	}
	mg.run_cg_iterations( lgs.A, lgs.x, lgs.b, lgs.r, 1000, 1.0e-10 );


	timer.stop();

	convergence.push_back(lgs.computeRMSE());
	convergence_time.push_back(timer.elapsedSeconds());



	return std::make_tuple(sys.stripBoundary(lgs.x), to_vector(convergence), to_vector(convergence_time));
	*/


	int maxIterations = 1000;
	auto result = mg.solve2(maxIterations);
	return std::make_tuple(sys.stripBoundary(std::get<0>(result)), std::get<1>(result), std::get<2>(result));
}


std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_gs(PNSystem& sys)
{
	Solver solver;
	setup_solver(solver, sys);
	auto result = solver.solve_gs();
	return std::make_tuple(sys.stripBoundary(std::get<0>(result)), std::get<1>(result), std::get<2>(result));
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_blockgs(PNSystem& sys)
{
	Solver solver;
	setup_solver(solver, sys);
	auto result = solver.solve_blockgs();
	return std::make_tuple(sys.stripBoundary(std::get<0>(result)), std::get<1>(result), std::get<2>(result));
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_cg(PNSystem& sys)
{
	Solver solver;
	setup_solver(solver, sys);
	auto result = solver.solve_cg();
	return std::make_tuple(sys.stripBoundary(std::get<0>(result)), std::get<1>(result), std::get<2>(result));
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_cg_eigen(PNSystem& sys)
{
	Solver solver;
	setup_solver(solver, sys);
	auto result = solver.solve_cg_eigen();
	return std::make_tuple(sys.stripBoundary(std::get<0>(result)), std::get<1>(result), std::get<2>(result));
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_sparseLU(PNSystem& sys)
{
	Solver solver;
	setup_solver(solver, sys);
	auto result = solver.solve_sparseLU();
	return std::make_tuple(sys.stripBoundary(std::get<0>(result)), std::get<1>(result), std::get<2>(result));
}


PYBIND11_MODULE(pnsolver, m)
{
	m.def( "solve_multigrid", &solve_multigrid);
	m.def( "solve_multigrid2", &solve_multigrid2);
	m.def( "solve_cg", &solve_cg);
	m.def( "solve_cg_eigen", &solve_cg_eigen);
	m.def( "solve_gs", &solve_gs);
	m.def( "solve_blockgs", &solve_blockgs);
	m.def( "solve_sparseLU", &solve_sparseLU);



	// Solver ================================
	py::class_<Solver> class_solver(m, "Solver");
	class_solver
	.def("__init__",
	[](Solver &m, int numLevels)
	{
		new (&m) Solver(numLevels);
	})
	.def("setMultigridLevel",&Solver::setMultigridLevel)
	.def("solve", &Solver::solve)
	.def("solve_gs", &Solver::solve_gs)
	.def("solve_cg", &Solver::solve_cg)
	;



	// PNSystem ==============================
	py::class_<PNSystem> class_pnsystem(m, "PNSystem");
	class_pnsystem
	.def("__init__",
	[](PNSystem &m, const std::string& stencil_name, Domain &domain, bool neumannBC)
	{
		PNSystem::Stencil stencil;

		if( stencil_name == "noop" )
			stencil = PNSystem::Stencil::noop();
		else
			stencil = PNSystem::findStencil(stencil_name);
		new (&m) PNSystem(stencil, domain, neumannBC);
	})
	.def("getNumCoefficients", &PNSystem::getNumCoefficients )
	.def("getResolution", &PNSystem::getResolution )
	.def("getNumVoxels", &PNSystem::getNumVoxels )
	.def("getOrder", &PNSystem::getOrder )
	.def("build", &PNSystem::build )
	.def("setField", &PNSystem::setField )
	.def("get_A_real", &PNSystem::get_A_real )
	.def("get_b_real", &PNSystem::get_b_real )
	.def("get_A_real_test", &PNSystem::get_A_real_test )
	.def("get_b_real_test", &PNSystem::get_b_real_test )
	.def("setDebugVoxel", &PNSystem::setDebugVoxel )
	.def("get_debug", &PNSystem::get_debug )

	//.def("computeGroundtruth", &PNSystem::computeGroundtruth )
	.def("getVoxelInfo", &PNSystem::getVoxelInfo )
	;

	// Domain ============================================================
	py::class_<Domain> class_domain(m, "Domain");
	class_domain
	.def("__init__",
	[](Domain &m, const Eigen::Matrix<double, 3, 1>& size, const Eigen::Matrix<int, 3, 1>& resolution, const Eigen::Matrix<double, 3, 1>& offset)
	{
		new (&m) Domain( size,
						 resolution,
						 offset);
	})
	.def("getResolution",
		[](Domain &m)
		{
			return m.getResolution();
		})
	.def("setResolution",
		[](Domain &m, const Eigen::Matrix<int, 3, 1>& resolution)
		{
			m.setResolution(resolution);
		})
	.def("getVoxelSize",
		[](Domain &m)
		{
			return m.getVoxelSize();
		})
	.def("numVoxels",
		[](Domain &m)
		{
			return m.numVoxels();
		})
	.def("worldToVoxel",
		[](Domain &m, const Eigen::Matrix<double, 3, 1>& pWS)
		{
			return m.worldToVoxel(pWS);
		})
	.def("voxelToWorld",
		[](Domain &m, const Eigen::Matrix<double, 3, 1>& pVS)
		{
			return m.voxelToWorld(pVS);
		})
	.def("getBoundMax",
		[](Domain &m)
		{
			return m.getBoundMax();
		})
	.def("getBoundMin",
		[](Domain &m)
		{
			return m.getBoundMin();
		})
	;


	// Field ============================================================
	py::class_<Field, Field::Ptr> class_field(m, "Field");
	class_field
	.def("__call__",
	[](Field &m, const Eigen::Matrix<double, 3, 1>& pWS)
	{
		return m.eval(pWS);
	})
	;

	// VoxelGrid ============================================================
	py::class_<VoxelGridField, VoxelGridField::Ptr> class_VoxelGridField(m, "VoxelGridField", class_field);
	class_VoxelGridField
	.def("__init__",
	[](VoxelGridField &m, py::array b, Domain& domain, const Eigen::Matrix<double, 3, 1>& offset)
	{
		typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> Strides;

		// Request a buffer descriptor from Python
		py::buffer_info info = b.request();

		// Some sanity checks ...
		if (info.format != py::format_descriptor<std::complex<double>>::format())
			throw std::runtime_error("Incompatible format: expected a complex array!");

		if (info.ndim != 3)
			throw std::runtime_error("Incompatible buffer dimension!");

		int res_x = int(info.shape[0]);
		int res_y = int(info.shape[1]);
		int res_z = int(info.shape[2]);

		auto data = static_cast<std::complex<double> *>(info.ptr);
		new (&m) VoxelGridField( data, domain, offset );
	})
	.def("test", &VoxelGridField::test)
	.def("getSlice", &VoxelGridField::getSlice)
	;

	// Constant ============================================================
	py::class_<Constant, Constant::Ptr> class_constant(m, "Constant", class_field);
	class_constant
	.def("__init__",
	[](Constant &m, std::complex<double> value)
	{
		new (&m) Constant(value);
	});

	// RadianceField ============================================================
	py::class_<RadianceField, RadianceField::Ptr> class_radiancefield(m, "RadianceField");
	class_radiancefield
	.def("__call__",
	[](RadianceField &m, const Eigen::Matrix<double, 3, 1>& pWS, const Eigen::Matrix<double, 3, 1>& omega)
	{
		return m.eval(pWS, omega);
	})
	;

	// SHEXP ============================================================
	py::class_<SHEXP, SHEXP::Ptr> class_shexp(m, "SHEXP", class_radiancefield);
	class_shexp
	.def("__init__",
	[](SHEXP &m, int order)
	{
		new (&m) SHEXP(order);
	})
	.def("setCoefficientField", &SHEXP::setCoefficientField)
	.def("getCoefficientField", &SHEXP::getCoefficientField)
	;
}
