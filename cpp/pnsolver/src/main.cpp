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

#include<common/Domain.h>
#include<field/VoxelGridField.h>
#include<field/Constant.h>
#include<field/SHEXP.h>
#include <PNSystem.h>

namespace py = pybind11;


Eigen::VectorXd to_vector( const std::vector<double>& values )
{
	Eigen::RowVectorXd x = Eigen::VectorXd( values.size() );
	for( int i=0;i<values.size();++i )
		x(i) = values[i];
	return x;
}


// returns the square root of the residual
double run_cg_iterations( PNSystem::RealMatrix& A, Eigen::VectorXd& b, Eigen::VectorXd& x, Eigen::VectorXd& r, int numIterations, double tol )
{
	r = b-A*x;
	Eigen::VectorXd p = r;
	Eigen::VectorXd Ap;
	double rsold = r.squaredNorm();

	for( int i=0;i<numIterations;++i )
	{
		Ap = A*p;
		double alpha = rsold/p.dot(Ap);
		x = x + alpha*p;
		r = r - alpha*Ap;
		double rsnew = r.squaredNorm();
		if( std::sqrt(rsnew) < tol )
			return std::sqrt(rsnew);
		p = r + (rsnew/rsold)*p;
		rsold = rsnew;
	}

	return std::sqrt(rsold);
}



// returns the solution vector x, and the convergence plot per iteration
std::tuple<Eigen::VectorXd, Eigen::VectorXd> solve_multigrid( PNSystem& sys_fine )
{
	Timer timer;
	timer.start();

	std::vector<double> solve_convergence;

	V3i res_fine = sys_fine.getDomain().getResolution();

	bool is2D = res_fine[2] == 1;

	// for multigrid, we require the resolution to be even,
	// so that we can do restriction and interpolation straigh forwardly
	if( (res_fine[0]%2!=0)||(res_fine[1]%2!=0)||(!is2D && (res_fine[2]%2!=0)))
		throw std::runtime_error("solve_multigrid multigrid currently requires even resolution");

	V3i res_coarse( res_fine[0]/2, res_fine[1]/2, is2D ? 1:res_fine[2]/2 );

	Domain domain_coarse( sys_fine.getDomain().getBound().getExtents(),
						  res_coarse,
						  sys_fine.getDomain().getBound().min );

	// build coarse level using restriction step on fields
	PNSystem sys_coarse(sys_fine.getStencil(), domain_coarse, sys_fine.getBoundaryConditions());
	sys_coarse.setFields( sys_fine.getFields().createRestricted() );
	sys_coarse.build();

	std::cout << "solve_multigrid: time for setup: " << timer.elapsedSeconds() << "s\n";

	//the following two lines solve the system on the fine grid and return the restricted result
	//NB: to visualize the solution, you have to use the correct coarse solution...
	//Eigen::VectorXd x_restricted = restricted( sys_coarse, solve_cg( get_b_real_test() ));
	//return sys_coarse.stripBoundary(x_restricted);

	//the following two lines solve the system on the coarse grid and return the interpolated result
	//Eigen::VectorXd x_interp = interpolate( sys_coarse, sys_coarse.solve_cg(sys_coarse.get_b_real_test()) );
	//return stripBoundary(x_interp);

	double tol = 1.0e-10; // convergence error tolerance
	int maxNumCycles = 1000; // maximum number of V-cycles
	int numSteps_fine = 10; // number of CG iterations on the fine level
	int numSteps_coarse = 10; // number of CG iterations on the coarse level

	PNSystem::RealMatrix A_coarse = sys_coarse.get_A_real().transpose()*sys_coarse.get_A_real();
	Eigen::VectorXd b_coarse( sys_coarse.getVoxelManager().getNumUnknowns() );
	Eigen::VectorXd x_coarse( sys_coarse.getVoxelManager().getNumUnknowns() );
	x_coarse.fill(0.0);
	Eigen::VectorXd r_coarse( sys_coarse.getVoxelManager().getNumUnknowns() );

	PNSystem::RealMatrix A_fine = sys_fine.get_A_real().transpose()*sys_fine.get_A_real();
	Eigen::VectorXd b_fine = sys_fine.get_A_real().transpose()*sys_fine.get_b_real();
	Eigen::VectorXd x_fine = b_fine; // initial guess for first fine level solve
	Eigen::VectorXd r_fine( sys_fine.getVoxelManager().getNumUnknowns() );


	// v-cycle
	solve_convergence.clear();
	for( int i=0;i<maxNumCycles;++i )
	{
		// solve on fine level
		double rmse = run_cg_iterations( A_fine, b_fine, x_fine, r_fine, numSteps_fine, tol );

		// check convergence
		solve_convergence.push_back(rmse);
		if(rmse < tol)
			// done
			break;

		// get residual as right hand side on the coarse level
		b_coarse = sys_fine.downsample(sys_coarse, r_fine);

		// solve correction equation on course level
		run_cg_iterations( A_coarse, b_coarse, x_coarse, r_coarse, numSteps_coarse, tol );

		// x_course is the error correction on the course level.
		// We interpolate it and apply it to the current final level solution.
		x_fine = x_fine + sys_fine.upsample(sys_coarse, x_coarse);
	}

	timer.stop();
	std::cout << "solve_multigrid: " << timer.elapsedSeconds() << "s\n";
	return std::make_tuple( sys_fine.stripBoundary(x_fine), to_vector(solve_convergence) );
}



std::tuple<Eigen::VectorXd, Eigen::VectorXd> solve_cg( PNSystem& sys )
{
	Timer timer;
	timer.start();

	std::vector<double> solve_convergence;

	PNSystem::RealMatrix A = sys.get_A_real().transpose()*sys.get_A_real();
	Eigen::VectorXd b = sys.get_A_real().transpose()*sys.get_b_real();
	Eigen::VectorXd x = b; // initial guess
	Eigen::VectorXd r( sys.getVoxelManager().getNumUnknowns() );



	// cg solver
	double tol = 1.0e-10; // convergence error tolerance
	int numIterations = b.rows();
	{
		r = b-A*x;
		Eigen::VectorXd p = r;
		Eigen::VectorXd Ap;
		double rsold = r.squaredNorm();
		for( int i=0;i<numIterations;++i )
		{
			Ap = A*p;
			double alpha = rsold/p.dot(Ap);
			x = x + alpha*p;
			r = r - alpha*Ap;
			double rsnew = r.squaredNorm();
			solve_convergence.push_back(std::sqrt(rsnew));
			if( std::sqrt(rsnew) < tol )
				break;
			p = r + (rsnew/rsold)*p;
			rsold = rsnew;
		}
	} // end of cg solver

	timer.stop();
	std::cout << "solve_cg: " << timer.elapsedSeconds() << "s\n";
	return std::make_tuple( sys.stripBoundary(x), to_vector(solve_convergence) );
}

PYBIND11_MODULE(pnsolver, m)
{
	m.def( "solve_multigrid", &solve_multigrid);
	m.def( "solve_cg", &solve_cg);


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
	.def("getNumVoxels", &PNSystem::getNumVoxels )
	.def("getOrder", &PNSystem::getOrder )
	.def("build", &PNSystem::build )
	.def("solve", &PNSystem::solve )
	.def("solve_cg", &PNSystem::solve_cg )
	.def("setField", &PNSystem::setField )
	.def("get_A_real", &PNSystem::get_A_real )
	.def("get_b_real", &PNSystem::get_b_real )
	.def("get_A_real_test", &PNSystem::get_A_real_test )
	.def("get_b_real_test", &PNSystem::get_b_real_test )
	.def("setDebugVoxel", &PNSystem::setDebugVoxel )
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
