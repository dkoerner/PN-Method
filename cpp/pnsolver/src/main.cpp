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


void buildUpAndDownsamplingMatrices( PNSystem& sys_fine, PNSystem::RealMatrix& downsampleMatrix, PNSystem::RealMatrix& upsampleMatrix )
{
	PNSystem::VoxelManager& vm_fine = sys_fine.getVoxelManager();
	PNSystem::VoxelManager vm_coarse = vm_fine.downsample();


	PNSystem::MatrixBuilderd downsampleMatrixBuilder;
	PNSystem::MatrixBuilderd upsampleMatrixBuilder;

	// downsample --------------------------


	// offset and weights for coefficients at grid with offset=1,1
	std::vector<std::pair<V3i, double>> stamp_11;
	stamp_11.push_back( std::make_pair(V3i(0,0,0), 0.75*0.75) );
	stamp_11.push_back( std::make_pair(V3i(0,1,0), 0.75*0.75) );
	stamp_11.push_back( std::make_pair(V3i(1,1,0), 0.75*0.75) );
	stamp_11.push_back( std::make_pair(V3i(1,0,0), 0.75*0.75) );
	stamp_11.push_back( std::make_pair(V3i(0,2,0), 0.25*0.75) );
	stamp_11.push_back( std::make_pair(V3i(1,2,0), 0.25*0.75) );
	stamp_11.push_back( std::make_pair(V3i(0,-1,0), 0.25*0.75) );
	stamp_11.push_back( std::make_pair(V3i(1,-1,0), 0.25*0.75) );
	stamp_11.push_back( std::make_pair(V3i(2,0,0), 0.25*0.75) );
	stamp_11.push_back( std::make_pair(V3i(2,1,0), 0.25*0.75) );
	stamp_11.push_back( std::make_pair(V3i(-1,0,0), 0.25*0.75) );
	stamp_11.push_back( std::make_pair(V3i(-1,1,0), 0.25*0.75) );
	stamp_11.push_back( std::make_pair(V3i(2,2,0), 0.25*0.25) );
	stamp_11.push_back( std::make_pair(V3i(-1,2,0), 0.25*0.25) );
	stamp_11.push_back( std::make_pair(V3i(-1,-1,0), 0.25*0.25) );
	stamp_11.push_back( std::make_pair(V3i(2,-1,0), 0.25*0.25) );

	// offset and weights for coefficients at grid with offset=0,1
	std::vector<std::pair<V3i, double>> stamp_01;
	stamp_01.push_back( std::make_pair(V3i(-1,-1,0), 0.5*0.25) );
	stamp_01.push_back( std::make_pair(V3i(-1,0,0), 0.5*0.75) );
	stamp_01.push_back( std::make_pair(V3i(-1,1,0), 0.5*0.75) );
	stamp_01.push_back( std::make_pair(V3i(-1,2,0), 0.5*0.25) );
	stamp_01.push_back( std::make_pair(V3i(0,-1,0), 0.25) );
	stamp_01.push_back( std::make_pair(V3i(0,0,0), 0.75) );
	stamp_01.push_back( std::make_pair(V3i(0,1,0), 0.75) );
	stamp_01.push_back( std::make_pair(V3i(0,2,0), 0.25) );
	stamp_01.push_back( std::make_pair(V3i(1,-1,0), 0.5*0.25) );
	stamp_01.push_back( std::make_pair(V3i(1,0,0), 0.5*0.75) );
	stamp_01.push_back( std::make_pair(V3i(1,1,0), 0.5*0.75) );
	stamp_01.push_back( std::make_pair(V3i(1,2,0), 0.5*0.25) );

	// offset and weights for coefficients at grid with offset=1,0
	std::vector<std::pair<V3i, double>> stamp_10;
	stamp_10.push_back( std::make_pair(V3i(-1,-1,0), 0.5*0.25) );
	stamp_10.push_back( std::make_pair(V3i(0,-1,0), 0.5*0.75) );
	stamp_10.push_back( std::make_pair(V3i(1,-1,0), 0.5*0.75) );
	stamp_10.push_back( std::make_pair(V3i(2,-1,0), 0.5*0.25) );
	stamp_10.push_back( std::make_pair(V3i(-1,0,0), 0.25) );
	stamp_10.push_back( std::make_pair(V3i(0,0,0), 0.75) );
	stamp_10.push_back( std::make_pair(V3i(1,0,0), 0.75) );
	stamp_10.push_back( std::make_pair(V3i(2,0,0), 0.25) );
	stamp_10.push_back( std::make_pair(V3i(-1,1,0), 0.5*0.25) );
	stamp_10.push_back( std::make_pair(V3i(0,1,0), 0.5*0.75) );
	stamp_10.push_back( std::make_pair(V3i(1,1,0), 0.5*0.75) );
	stamp_10.push_back( std::make_pair(V3i(2,1,0), 0.5*0.25) );

	// offset and weights for coefficients at grid with offset=0,0
	std::vector<std::pair<V3i, double>> stamp_00;
	stamp_00.push_back( std::make_pair(V3i(0,0,0), 1.0) );
	stamp_00.push_back( std::make_pair(V3i(0,-1,0), 0.5) );
	stamp_00.push_back( std::make_pair(V3i(0,1,0), 0.5) );
	stamp_00.push_back( std::make_pair(V3i(-1,0,0), 0.5) );
	stamp_00.push_back( std::make_pair(V3i(1,0,0), 0.5) );
	stamp_00.push_back( std::make_pair(V3i(1,-1,0), 0.5*0.5) );
	stamp_00.push_back( std::make_pair(V3i(1,1,0), 0.5*0.5) );
	stamp_00.push_back( std::make_pair(V3i(-1,1,0), 0.5*0.5) );
	stamp_00.push_back( std::make_pair(V3i(-1,-1,0), 0.5*0.5) );




	std::vector<std::vector<std::pair<V3i, double>>*> stamps(8, 0);
	stamps[0] = &stamp_00;
	stamps[1] = &stamp_10;
	stamps[2] = &stamp_11;
	stamps[3] = &stamp_01;



	// 2d --------------------
	// iterate all coarse voxels
	std::vector<PNSystem::Voxel>& voxels_coarse = vm_coarse.getVoxels();
	for( auto&v_coarse:voxels_coarse )
	{
		// coordinate of the corresponding fine voxel
		V3i coord_fine = v_coarse.coord*2;

		// Now compute the value of the coarse voxels coefficients from gathering the values from
		// the fine voxels which are affecting the coarse voxels values. These values will be weighted,
		// as fine voxels are affected by multiple coarse voxels.

		for( int coeff_index=0;coeff_index<sys_fine.getNumCoefficients();++coeff_index )
		{
			//int coeff_index = 0;

			// get the global index of the coarse voxel/coeff
			int gi_coarse = vm_coarse.getGlobalIndex( v_coarse, coeff_index );

			// initialize the value of the current coefficient at current coarse voxel
			// this value is computed from all affected fine voxels
			double weight_sum = 0.0;

			// now iterate over all affected fine voxels. This depends on the grid location
			int grid_index = sys_fine.getStencil().getGridIndexFromCoefficient(coeff_index);
			auto stamp = stamps[grid_index];
			std::vector<std::pair<int, double>> row_entries; // column index and value
			for( auto& tt:*stamp )
			{
				V3i v_fine_coord = coord_fine + std::get<0>(tt);
				double weight = std::get<1>(tt);

				// using the stamps offset on boundary voxels may create invalid voxel coordinates
				if( !vm_fine.voxelIsValid(v_fine_coord) )
					continue;

				PNSystem::Voxel& v_fine = vm_fine.getVoxel(v_fine_coord);

				// retrieve global index of fine voxel/coeff
				int gi_fine = vm_fine.getGlobalIndex(v_fine, coeff_index);
				if( gi_fine >= 0 )
				{
					// either the fine voxel/coeff is not a boundary voxel/coeff,
					// or we have a boundary voxel with neumann BC, pointing to an
					// interior voxel
					row_entries.push_back( std::make_pair(gi_fine, weight) );
				}
				// only set weights on voxels which are actually active
				if( !vm_fine.isBoundaryCoefficient(v_fine, coeff_index)&&(gi_coarse>=0) )
					upsampleMatrixBuilder.coeff(gi_fine, gi_coarse) += weight;

				// else: value_fine remains zero, which means we have a boundary voxel with dirichlet BC
				weight_sum += weight;
			} // for each stamp entry

			// build row
			if( !vm_coarse.isBoundaryCoefficient(v_coarse, coeff_index) )
				for( auto& it:row_entries )
				{
					int gi_fine = it.first;
					double weight = it.second;
					downsampleMatrixBuilder.coeff(gi_coarse, it.first) += weight/weight_sum;
				}
		} // for each coefficient
	} // for each coarse voxel



	//return x_coarse;
	downsampleMatrix = downsampleMatrixBuilder.build( vm_coarse.getNumUnknowns(), vm_fine.getNumUnknowns() );
	upsampleMatrix = upsampleMatrixBuilder.build( vm_fine.getNumUnknowns(), vm_coarse.getNumUnknowns() );
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
std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_multigrid( PNSystem& sys_fine )
{
	Timer timer;
	PNSystem::RealMatrix downsample;
	PNSystem::RealMatrix upsample;

	buildUpAndDownsamplingMatrices(sys_fine, downsample, upsample);

	std::vector<double> solve_convergence;
	std::vector<double> solve_convergence_timestamps;



	// build coarse level using restriction step on fields
	PNSystem sys_coarse(sys_fine.getStencil(), sys_fine.getDomain().downsample(), sys_fine.getBoundaryConditions());
	sys_coarse.setFields( sys_fine.getFields().createRestricted() );
	sys_coarse.build();

	std::cout << "solve_multigrid: time for setup: " << timer.elapsedSeconds() << "s\n";

	timer.start();
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
		solve_convergence_timestamps.push_back(timer.elapsedSeconds());
		//std::cout << "i=" << i << " rmse=" << rmse << std::endl;
		if(rmse < tol)
			// done
			break;

		// get residual as right hand side on the coarse level
		//b_coarse = sys_fine.downsample(sys_coarse.getVoxelManager(), r_fine);
		b_coarse = downsample*r_fine;

		// solve correction equation on course level
		run_cg_iterations( A_coarse, b_coarse, x_coarse, r_coarse, numSteps_coarse, 0.0 );

		// x_coarse is the error correction on the course level.
		// We interpolate it and apply it to the current final level solution.
		//x_fine = x_fine + sys_fine.upsample(sys_coarse.getVoxelManager(), x_coarse);
		x_fine = x_fine + upsample*x_coarse;
	}

	timer.stop();
	std::cout << "solve_multigrid: " << timer.elapsedSeconds() << "s numIterations=" << solve_convergence.size() <<  "\n";
	return std::make_tuple( sys_fine.stripBoundary(x_fine),
							to_vector(solve_convergence),
							to_vector(solve_convergence_timestamps));
}



std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_cg( PNSystem& sys )
{
	Timer timer;
	timer.start();

	std::vector<double> solve_convergence;
	std::vector<double> solve_convergence_timestamps;

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
			solve_convergence_timestamps.push_back(timer.elapsedSeconds());
			if( std::sqrt(rsnew) < tol )
				break;
			p = r + (rsnew/rsold)*p;
			rsold = rsnew;
		}
	} // end of cg solver

	timer.stop();
	std::cout << "solve_cg: " << timer.elapsedSeconds() << "s\n";
	return std::make_tuple( sys.stripBoundary(x),
							to_vector(solve_convergence),
							to_vector(solve_convergence_timestamps));
}


std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_sparseLU( PNSystem& sys )
{
	std::vector<double> solve_convergence;
	std::vector<double> solve_convergence_timestamps;


	Timer timer;
	timer.start();





	//Eigen::ConjugateGradient<RealMatrix> solver;
	//Eigen::BiCGSTAB<RealMatrix> solver;
	Eigen::SparseLU<PNSystem::RealMatrix> solver;
	solver.compute(sys.get_A_real());
	if(solver.info()!=Eigen::Success)
	{
		throw std::runtime_error("solve_sparseLU decomposition failed");
	}
	Eigen::VectorXd x = solver.solve(sys.get_b_real());
	if(solver.info()!=Eigen::Success)
	{
		throw std::runtime_error("solve_sparseLU solve failed");
	}


	timer.stop();
	std::cout << "solve_sparseLU: " << timer.elapsedSeconds() << "s\n";

	return std::make_tuple( sys.stripBoundary(x),
							to_vector(solve_convergence),
							to_vector(solve_convergence_timestamps));
}

PYBIND11_MODULE(pnsolver, m)
{
	m.def( "solve_multigrid", &solve_multigrid);
	m.def( "solve_cg", &solve_cg);
	m.def( "solve_sparseLU", &solve_sparseLU);


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
