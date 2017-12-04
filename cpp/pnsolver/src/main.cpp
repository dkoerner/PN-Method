#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <complex>
#include <iostream>

#include <math/common.h>
#include <math/vector.h>
#include <math/ray.h>
#include <math/sph.h>

#include <util/timer.h>
#include <util/mem.h>

#include<common/Domain.h>

#include<fields/VoxelGridField.h>

#include <PNSystem.h>
#include <PNSolution.h>
#include <PNVolume.h>

#include <solver.h>
#include <MultigridSolver.h>


namespace py = pybind11;

void save_solution( const std::string& filename, PNSystem& sys, const Eigen::VectorXd& x );


void setup_solver( Solver& mg, PNSystem& sys, int numLevels = 1 )
{
	mg = Solver(numLevels);

	PNSystem::Stencil& stencil = sys.getStencil();
	int boundaryConditions = sys.getBoundaryConditions();
	Domain domain = sys.getDomain();
	//PNSystem::Fields fields = sys.getFields();

	for( int i=0;i<numLevels;++i )
	{
		/*
		Eigen::VectorXd x;
		PNSystem::RealMatrix downsample;
		PNSystem::RealMatrix upsample;

		PNSystem sys_level(stencil, domain, boundaryConditions);
		//sys_level.setFields( fields );
		sys_level.build();

		std::cout << "fetching A...\n";std::flush(std::cout);
		PNSystem::RealMatrix& A_org = sys_level.get_A_real();
		std::cout << "computing A^T...\n";std::flush(std::cout);
		PNSystem::RealMatrix AT = A_org.transpose();
		std::cout << "computing A^T*A...\n";std::flush(std::cout);
		mg.m_levels[i].A = AT*A_org;


		//std::cout << "computing b...\n";std::flush(std::cout);
		mg.m_levels[i].b = AT*sys_level.get_b_real();
		//std::cout << "warning: not using ATA form\n";
		//A = sys_level.get_A_real();
		//Eigen::VectorXd b = sys_level.get_b_real();

		//std::cout << "computing x...\n";std::flush(std::cout);
		mg.m_levels[i].x = Eigen::VectorXd(sys_level.getVoxelManager().getNumUnknowns());
		mg.m_levels[i].x.fill(0.0);

		// used for blockgs test
		//mg.m_levels[i].m_voxelManager = sys_level.getVoxelManager();

		// downsample to next level
		if(i<numLevels-1)
		{
			buildUpAndDownsamplingMatrices(sys_level, mg.m_levels[i].downsample, mg.m_levels[i].upsample);
			domain = domain.downsample();
			//fields = fields.createRestricted();
		}

		//std::cout << "setting multigrid level...\n";std::flush(std::cout);
		//mg.setMultigridLevel(i, A, x, b, downsample, upsample);
		*/
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


	return std::make_tuple(std::get<0>(result), std::get<1>(result), std::get<2>(result));
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

	return std::make_tuple(std::get<0>(result), std::get<1>(result), std::get<2>(result));
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



	return std::make_tuple(lgs.x, to_vector(convergence), to_vector(convergence_time));
	*/


	int maxIterations = 1000;
	auto result = mg.solve2(maxIterations);
	return std::make_tuple(std::get<0>(result), std::get<1>(result), std::get<2>(result));
}

/*
std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_gs(PNSystem& sys)
{
	Solver solver;
	setup_solver(solver, sys);
	auto result = solver.solve_gs();
	return std::make_tuple(std::get<0>(result), std::get<1>(result), std::get<2>(result));
}
*/

std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_blockgs(PNSystem& sys)
{
	Solver solver;
	setup_solver(solver, sys);
	auto result = solver.solve_blockgs();
	return std::make_tuple(std::get<0>(result), std::get<1>(result), std::get<2>(result));
}

// this will solve the normal form of Ax=b by explicitly computing A^TA and A^Tb and using CG to solve it
std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_ls_cg(PNSystem& sys, double tol )
{
	sys.build();

	std::vector<double> list_rmse;
	std::vector<double> list_time;

	Timer timer;

	std::cout << "computing A^T...\n";std::flush(std::cout);
	PNSystem::RealMatrix AT = sys.get_A().transpose();


	// we solve using ATA
	std::cout << "computing A^T*A...\n";std::flush(std::cout);
	PNSystem::RealMatrix A = AT*sys.get_A();
	Eigen::VectorXd b = AT*sys.get_b();

	Eigen::VectorXd x = Eigen::VectorXd(b.rows()); // initial guess
	x.fill(0.0);
	Eigen::VectorXd r( x.rows() );


	std::cout << "solving...\n";std::flush(std::cout);

	//double tol = 1.0e-10; // convergence error tolerance
	int iteration = 0;
	int numIterations = b.rows();

	timer.start();
	// cg solver
	{
		r = b-A*x;
		Eigen::VectorXd p = r;
		Eigen::VectorXd Ap;
		double rsold = r.squaredNorm();
		for( iteration=0;iteration<numIterations;++iteration )
		{
			Ap = A*p;
			double alpha = rsold/p.dot(Ap);
			x = x + alpha*p;
			r = r - alpha*Ap;
			double rsnew = r.squaredNorm();
			//std::cout << "rmse=" << std::sqrt(rsnew) << std::endl;
			list_rmse.push_back(std::sqrt(rsnew));
			list_time.push_back(timer.elapsedSeconds());

			if( std::sqrt(rsnew) < tol )
				break;
			p = r + (rsnew/rsold)*p;
			rsold = rsnew;
		}
	} // end of cg solver
	timer.stop();
	std::cout << "solve_ls_cg: " << timer.elapsedSeconds() << "s #iterations=" << iteration << " rmse=" << (b - A*x).norm() <<  "\n";

	return std::make_tuple( x,
							to_vector(list_rmse),
							to_vector(list_time));
}


// this will solve Ax=b using CG. Note that for PN-problems this will not work as A is not symmetric
// However, it will work for problems which produce a symmetric, positive, definite matrix, such as classical diffusion
std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_cg(PNSystem& sys, double tol )
{
	sys.build();

	Timer timer;

	PNSystem::RealMatrix& A = sys.get_A();
	Eigen::VectorXd b = sys.get_b();

	Eigen::VectorXd x = Eigen::VectorXd(b.rows()); // initial guess
	x.fill(0.0);
	Eigen::VectorXd r( x.rows() );


	std::cout << "solving...\n";std::flush(std::cout);

	//double tol = 1.0e-10; // convergence error tolerance
	int iteration = 0;
	int numIterations = b.rows();

	timer.start();
	// cg solver
	{
		r = b-A*x;
		Eigen::VectorXd p = r;
		Eigen::VectorXd Ap;
		double rsold = r.squaredNorm();
		for( iteration=0;iteration<numIterations;++iteration )
		{
			Ap = A*p;
			double alpha = rsold/p.dot(Ap);
			x = x + alpha*p;
			r = r - alpha*Ap;
			double rsnew = r.squaredNorm();
			std::cout << "rmse=" << std::sqrt(rsnew) << std::endl;
			if( std::sqrt(rsnew) < tol )
				break;
			p = r + (rsnew/rsold)*p;
			rsold = rsnew;
		}
	} // end of cg solver
	timer.stop();
	std::cout << "solve_cg: " << timer.elapsedSeconds() << "s #iterations=" << iteration << "\n";

	return std::make_tuple( x,
							Eigen::VectorXd(),
							Eigen::VectorXd());
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_cg_eigen( PNSystem& sys, double tol )
{
	sys.build();

	std::cout << "solving cg eigen\n";std::flush(std::cout);

	std::vector<double> solve_convergence;
	std::vector<double> solve_convergence_timestamps;

	//double rmse = (sys.get_b() - sys.get_A()*x).norm();
	//solve_convergence.push_back(rmse);
	//solve_convergence_timestamps.push_back(0.0);


	Timer timer;
	timer.start();

	//Eigen::SparseLU<PNSystem::RealMatrix> solver;
	//typedef Eigen::ConjugateGradient<PNSystem::RealMatrix,Eigen::Lower, Eigen::IncompleteCholesky<double>> ICCG;
	//ICCG solver;

	Eigen::ConjugateGradient<PNSystem::RealMatrix> solver;
	solver.setTolerance(tol);
	solver.compute(sys.get_A());
	if(solver.info()!=Eigen::Success)
	{
		throw std::runtime_error("solve_cg_eigen decomposition failed");
	}
	Eigen::VectorXd x = solver.solve(sys.get_b());
	if(solver.info()!=Eigen::Success)
	{
		throw std::runtime_error("solve_cg_eigen solve failed");
	}

	timer.stop();

	double rmse = (sys.get_b() - sys.get_A()*x).norm();
	solve_convergence.push_back(rmse);
	solve_convergence_timestamps.push_back(timer.elapsedSeconds());


	std::cout << "solve_cg_eigen: " << timer.elapsedSeconds() << "s #iterations=" <<  solver.iterations() << "\n";
	return std::make_tuple( x,
							to_vector(solve_convergence),
							to_vector(solve_convergence_timestamps));
}

// solves the normal form A^TAx = A^Tb using eigens CG solver
std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_ls_cg_eigen( PNSystem& sys, double tol )
{
	sys.build();

	std::cout << "computing A^T...\n";std::flush(std::cout);
	PNSystem::RealMatrix AT = sys.get_A().transpose();


	// we solve using ATA
	std::cout << "computing A^T*A...\n";std::flush(std::cout);
	PNSystem::RealMatrix A = AT*sys.get_A();
	Eigen::VectorXd b = AT*sys.get_b();


	std::cout << "solving ls cg eigen\n";std::flush(std::cout);

	std::vector<double> solve_convergence;
	std::vector<double> solve_convergence_timestamps;

	//double rmse = (sys.get_b() - sys.get_A()*x).norm();
	//solve_convergence.push_back(rmse);
	//solve_convergence_timestamps.push_back(0.0);


	Timer timer;
	timer.start();

	//Eigen::SparseLU<PNSystem::RealMatrix> solver;
	//typedef Eigen::ConjugateGradient<PNSystem::RealMatrix,Eigen::Lower, Eigen::IncompleteCholesky<double>> ICCG;
	//ICCG solver;

	Eigen::ConjugateGradient<PNSystem::RealMatrix> solver;
	solver.setTolerance(tol);
	solver.compute(A);
	if(solver.info()!=Eigen::Success)
	{
		throw std::runtime_error("solve_cg_eigen decomposition failed");
	}
	Eigen::VectorXd x = solver.solve(b);
	if(solver.info()!=Eigen::Success)
	{
		throw std::runtime_error("solve_cg_eigen solve failed");
	}

	timer.stop();

	double rmse = (sys.get_b() - sys.get_A()*x).norm();
	solve_convergence.push_back(rmse);
	solve_convergence_timestamps.push_back(timer.elapsedSeconds());


	std::cout << "solve_ls_cg_eigen: " << timer.elapsedSeconds() << "s #iterations=" <<  solver.iterations() << "\n";
	return std::make_tuple( x,
							to_vector(solve_convergence),
							to_vector(solve_convergence_timestamps));
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_gs(PNSystem& sys, double tol, int maxIterations )
{
	sys.build();

	Timer timer;

	// here we convert the sparse matrix A to a row major A matrix. This allows us to iterate
	// sparse elements in a single row efficiently via eigen
	Eigen::SparseMatrix<double, Eigen::RowMajor> A = sys.get_A();
	Eigen::VectorXd b = sys.get_b();

	Eigen::VectorXd x = Eigen::VectorXd(b.rows()); // initial guess
	x.fill(0.0);
	Eigen::VectorXd r( x.rows() );


	std::cout << "solving...\n";std::flush(std::cout);

	int iteration = 0;

	timer.start();

	for( iteration=0;iteration<maxIterations;++iteration )
	{
		double rmse = (b - A*x).norm();
		//solve_convergence.push_back(rmse);
		//solve_convergence_timestamps.push_back(timer.elapsedSeconds());
		if( rmse < tol )
			break;

		//std::cout << "rmse=" << rmse << std::endl;

		// iterate all rows
		for( int i=0;i<A.rows();++i )
		{
			double aii = 0.0;
			double sum = b(i);
			// iterate all non-zero column elements (this is why we need RowMajor storage order for A)
			for( Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A, i);it;++it )
			{
				int j = it.col();
				if( i != j )
					sum -= it.value()*x.coeffRef(j);
				else
					aii = it.value();
			}
			x.coeffRef(i) = sum/aii;
		}
	}

	timer.stop();
	std::cout << "solve_gs: " << timer.elapsedSeconds() << "s #iterations=" << iteration << " rmse=" << (b - A*x).norm() <<  "\n";

	return std::make_tuple( x,
							Eigen::VectorXd(),
							Eigen::VectorXd());
}



// this solves the normal form for Ax=b by using a modified CG-solver which doesnt require an explicit
// computation of A^TA and A^Tb
std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_ls_lscg(PNSystem& sys, double tol)
{
	sys.build();

	Timer timer;


	PNSystem::RealMatrix& A = sys.get_A();
	Eigen::VectorXd b = sys.get_b();

	// we solve using ATA
	PNSystem::RealMatrix AT = A.transpose();

	Eigen::VectorXd x = Eigen::VectorXd(b.rows()); // initial guess
	x.fill(0.0);
	Eigen::VectorXd r( x.rows() );


	std::cout << "solving...\n";std::flush(std::cout);

	//double tol = 1.0e-10; // convergence error tolerance
	int iteration = 0;
	int numIterations = b.rows();

	timer.start();
	// cg solver
	{
		r = AT*(b-A*x);
		Eigen::VectorXd p = r;
		Eigen::VectorXd Ap;
		double rsold = r.squaredNorm();
		for( iteration=0;iteration<numIterations;++iteration )
		{
			Ap = AT*(A*p);
			double alpha = rsold/p.dot(Ap);
			x = x + alpha*p;
			r = r - alpha*Ap;
			double rsnew = r.squaredNorm();
			if( iteration % 10 == 0 )
				std::cout << "rmse=" << std::sqrt(rsnew) << std::endl;
			if( std::sqrt(rsnew) < tol )
				break;
			p = r + (rsnew/rsold)*p;
			rsold = rsnew;
		}
	} // end of cg solver
	timer.stop();
	std::cout << "solve_ls_lscg: " << timer.elapsedSeconds() << "s #iterations=" << iteration << "\n";

	return std::make_tuple( x,
							Eigen::VectorXd(),
							Eigen::VectorXd());
}



// solves Ax=b using a multigrid solver with gauss-seidel smoothing steps
std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_mg(PNSystem& sys, double tol, int maxIterations )
{
	MultigridSolver mgs;

	PNSystem::Stencil& stencil = sys.getStencil();
	PNVolume::Ptr problem = sys.getProblem();
	int boundaryConditions = sys.getBoundaryConditions();

	// first we setup the levels of the multigrid solver ---
	//int numLevels = 2;
	//for( int i=0;i<numLevels;++i )
	while(problem)
	{
		MultigridSolver::Level* level = mgs.addLevel();

		PNSystem sys_level(stencil, problem, boundaryConditions);
		sys_level.build();

		// debug
		level->sys_level = new PNSystem(stencil, problem, boundaryConditions);
		level->sys_level->build();


		//std::cout << "computing b...\n";std::flush(std::cout);
		//mg.m_levels[i].b = AT*sys_level.get_b_real();
		//std::cout << "warning: not using ATA form\n";
		//A = sys_level.get_A_real();
		//Eigen::VectorXd b = sys_level.get_b_real();

		level->A = sys_level.get_A();
		//level->b = MultigridSolver::RealVector(level->A.cols());
		//level->b.fill(0.0);
		level->b = sys_level.get_b();
		level->x = MultigridSolver::RealVector(level->A.cols());
		level->x.fill(0.0);

		//if(level->index == 0)
		//	level->b = sys_level.get_b();

		/*
		if(i==1)
		{
			MultigridSolver::Level* level = mgs.getLevel(i);
			Eigen::SparseLU<MultigridSolver::RealMatrix> solver;
			solver.compute(level->A);
			if(solver.info()!=Eigen::Success)
			{
				throw std::runtime_error("decomposition failed");
			}
			level->x = solver.solve(level->b);
			if(solver.info()!=Eigen::Success)
			{
				throw std::runtime_error("solve failed");
			}

			//std::cout << "rmse=" << level->computeRMSE() << std::endl;
			//save_solution( "test.pns", sys_level, level->x );

			std::cout << level->x.cols() << " " << level->x.rows() << std::endl;

			std::cout << "rmse(before)=" << mgs.getLevel(i-1)->computeRMSE() << std::endl;
			mgs.getLevel(i-1)->x = mgs.getLevel(i-1)->upsample*level->x;
			std::cout << mgs.getLevel(i-1)->x.cols() << " " << mgs.getLevel(i-1)->x.rows() << std::endl;
			std::cout << "rmse(after)=" << mgs.getLevel(i-1)->computeRMSE() << std::endl;
			//save_solution( "test.pns", *mgs.getLevel(i-1)->sys_level, mgs.getLevel(i-1)->x );
			save_solution( "test.pns", *mgs.getLevel(i)->sys_level, mgs.getLevel(i-1)->downsample*mgs.getLevel(i-1)->x );
		}
		*/



		//level->b = sys_level.get_b();

		//level->x = MultigridSolver::RealVector(sys_level.getVoxelManager().getNumUnknowns());
		//level->x.fill(0.0);


		// downsample to next coarser level
		problem = problem->downsample();
		if(problem)
		{
			// we have a coarser level, so create matrices for up- and downsampling the solution vector x
			buildUpAndDownsamplingMatrices2(sys_level, level->downsample, level->upsample);
		}
	}


	/*
	{
		MultigridSolver::Level* level = mgs.getLevel(0);
		Eigen::SparseLU<MultigridSolver::RealMatrix> solver;
		solver.compute(level->A);
		if(solver.info()!=Eigen::Success)
		{
			throw std::runtime_error("decomposition failed");
		}
		level->x = solver.solve(level->b);
		if(solver.info()!=Eigen::Success)
		{
			throw std::runtime_error("solve failed");
		}

		std::cout << "rmse=" << level->computeRMSE() << std::endl;

		return std::make_tuple( level->x,
								Eigen::VectorXd(),
								Eigen::VectorXd());
	}
	*/


	std::cout << "MultigridSolver numLevels=" << mgs.getNumLevels() << std::endl;


	auto result = mgs.solve(maxIterations);
	//save_solution( "test.pns", *mgs.getLevel(0)->sys_level, mgs.getLevel(0)->x );
	return result;

	return std::make_tuple( Eigen::VectorXd(),
							Eigen::VectorXd(),
							Eigen::VectorXd());
}




/*
std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_cg_eigen(PNSystem& sys)
{
	Solver solver;
	setup_solver(solver, sys);
	auto result = solver.solve_cg_eigen();
	return std::make_tuple(std::get<0>(result), std::get<1>(result), std::get<2>(result));
}
*/

std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_sparseLU(PNSystem& sys)
{
	Solver solver;
	setup_solver(solver, sys);
	auto result = solver.solve_sparseLU();
	return std::make_tuple(std::get<0>(result), std::get<1>(result), std::get<2>(result));
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_lscg(PNSystem& sys)
{
	sys.build();

	std::cout << "solve_lscg\n";

	std::vector<double> solve_convergence;
	std::vector<double> solve_convergence_timestamps;

	solve_convergence.push_back(sys.get_b().norm());
	solve_convergence_timestamps.push_back(0.0);


	Timer timer;
	timer.start();

	Eigen::LeastSquaresConjugateGradient<PNSystem::RealMatrix> solver;
	//solver.setMaxIterations(1);
	solver.compute(sys.get_A());
	if(solver.info()!=Eigen::Success)
	{
		throw std::runtime_error("solve_lscg decomposition failed");
	}
	Eigen::VectorXd x = solver.solve(sys.get_b());
	if(solver.info()!=Eigen::Success)
	{
		throw std::runtime_error("solve_lscg solve failed");
	}

	timer.stop();
	std::cout << "solve_lscg: " << timer.elapsedSeconds() << "s #iterations=" <<  solver.iterations() << "\n";

	solve_convergence.push_back((sys.get_b() - sys.get_A()*x).norm());
	solve_convergence_timestamps.push_back(timer.elapsedSeconds());


	return std::make_tuple(x, to_vector(solve_convergence), to_vector(solve_convergence_timestamps));
}





PNSystem::MatrixBuildercd::Matrix buildBlockDiagonalMatrix( const Eigen::MatrixXcd& M, int numBlocks )
{
	if( M.rows() != M.cols() )
		throw std::runtime_error("buildBlockDiagonalMatrix expected square matrix");


	int numElementsPerBlock = M.rows();
	PNSystem::MatrixBuildercd::Matrix M_sparse = M.sparseView();
	PNSystem::MatrixBuildercd builder;
	for(int k=0; k<M_sparse.outerSize(); ++k)
	{
		// Iterate over inside
		for(PNSystem::MatrixBuildercd::Matrix::InnerIterator it (M_sparse,k); it; ++it)
		{
			for( int block_i=0;block_i<numBlocks;++block_i )
			{
				int global_i = block_i*numElementsPerBlock + it.row();
				int global_j = block_i*numElementsPerBlock + it.col();
				builder.coeff(global_i,global_j) += it.value();
			}
		}
	}

	return builder.build(numBlocks*numElementsPerBlock, numBlocks*numElementsPerBlock);
}

// this function creates the S matrix from the starmap paper. It can be used to convert complex
// coefficients to real coefficients. The problem here is that it differs from sph::complex2RealConversionMatrix
// which requires us to first convert from real to complex using the inverse of S (from starmap) and then convert
// from complex to real using sph::complex2RealConversionMatrix.
// TODO: use sph::complex2RealConversionMatrix as S during stencil generation, in order to get rid of any conversion
// alltogether.
Eigen::MatrixXcd createComplexToRealConversionMatrix( int order )
{
	int numCoeffs = (order+1)*(order+1);

	Eigen::MatrixXcd S = Eigen::MatrixXcd::Zero(numCoeffs, numCoeffs);


	int count = 0;
	for( int l=0;l<=order;++l )
		for( int m=l;m>=0;--m )
		{
			// computes the real part coefficients for a row (defined by l,m) in the S matrix
			// (see bottom of p.5 in the starmap paper)
			if( m==0 )
				S.coeffRef(count, sph::index(l, m)) = 1.0;
			else
			{
				S.coeffRef(count, sph::index(l, m)) = std::pow(-1.0, double(m))/std::sqrt(2.0);
				S.coeffRef(count, sph::index(l, -m)) = std::pow(-1.0, 2.0*m)/std::sqrt(2.0);
			}

			++count;

			// computes the imaginary part coefficients for a row (defined by l,m) in the S matrix
			// (see bottom of p.5 in the starmap paper)
			if( m>0 )
			{
				S.coeffRef(count, sph::index(l, m)) = std::complex<double>(0.0, std::pow(-1.0, double(m))/std::sqrt(2.0));
				S.coeffRef(count, sph::index(l, -m)) = std::complex<double>(0.0, -std::pow(-1.0, 2.0*m)/std::sqrt(2.0));
				++count;
			}
		}

	return S;
}

Eigen::MatrixXcd createRealToComplexConversionMatrix( int order )
{
	return createComplexToRealConversionMatrix( order ).inverse();
}

// the solution vector x is the result from solving the system represented by sys.
// in order to use it for rendering we have to apply the following steps:
// -remove the staggering: change coefficients, such that they are all defined at the cell centers
// -get rid of coefficients which live in boundary voxels outside the domain (these were required with staggered grids)
// -convert coefficients from real to complex, to get the true SH coefficients of the radiance field (see p5 in starmap paper)
void save_solution( const std::string& filename, PNSystem& sys, const Eigen::VectorXd& x )
{
	if(sys.getResolution()[2] == 1)
	{
		std::cout << "save_solution: unable to create and save a PNSolution file for 2d problem due to different number of coefficients.\n";
		return;
	}

	// the solution is given in real valued coefficients. However, the complex->real conversion matrix used during stencil
	// generation is that from the starmap paper, which is different than our conversion matrix. Therefore we need to first
	// convert from real to complex using the inverse starmap transform, and then we can use our conversion matrix to convert
	// from complex to real again.This requires us to first convert from real to complex using the inverse of S (from starmap)
	// and then convert from complex to real using sph::complex2RealConversionMatrix.

	// TODO: use our own conversion matrix during stencil generation straight away...
	Eigen::VectorXcd x_complex =buildBlockDiagonalMatrix( createRealToComplexConversionMatrix(sys.getStencil().order), sys.getNumVoxels() )*
								sys.removeStaggering(x);

	Eigen::VectorXd x_real = (buildBlockDiagonalMatrix( sph::buildComplexToRealConversionMatrix(sys.getStencil().order), sys.getNumVoxels() )*
							 x_complex).real();
	PNSolution solution( sys.getOrder(), sys.getResolution(), Box3d(sys.getDomain().getBoundMin(), sys.getDomain().getBoundMax()), x_real.data() );
	solution.save(filename);
}


PNSolution::Ptr load_solution( const std::string& filename )
{
	return std::make_shared<PNSolution>( PNSolution(filename) );
}


VoxelGridField3d::Ptr load_voxelgridfield3d(const std::string& filename)
{
	return std::make_shared<VoxelGridField3d>(filename);
}


/*
Eigen::VectorXd getSolutionVector( PNSolution::Ptr pns )
{
	V3i resolution = pns->getResolution();
	int numVoxels = resolution[0]*resolution[1]*resolution[2];
	return (buildBlockDiagonalMatrix( createComplexToRealConversionMatrix(pns->getOrder()), numVoxels )*
									 Eigen::Map<Eigen::VectorXcd>( pns->data(), numVoxels*pns->getNumCoeffs() )).real();

}
*/


py::array getCoefficientArray( PNSystem& sys, const Eigen::VectorXd& x )
{
	Eigen::VectorXd x_unstaggered = sys.removeStaggering(x);
	int numCoeffs = sys.getNumCoefficients();
	V3i res = sys.getResolution();
	std::vector<double> data(res[0]*res[1]*res[2]*numCoeffs);

	for( int i=0;i<res[0];++i )
		for( int j=0;j<res[1];++j )
			for( int k=0;k<res[2];++k )
				for( int l=0;l<numCoeffs;++l )
				{
					int index_dst = (i*res[2]*res[1] + j*res[2] + k)*numCoeffs + l;
					int index_src = (i*res[2]*res[1] + j*res[2] + k)*numCoeffs + l;
					data[index_dst] = x_unstaggered.data()[index_src];
				}

	py::array b( //py::dtype(std::string("float64")),
				py::dtype::of<double>(),
			   {res[0], res[1], res[2], numCoeffs},
			   {int(sizeof(double))*res[2]*res[1]*numCoeffs, int(sizeof(double))*res[2]*numCoeffs, int(sizeof(double))*numCoeffs, int(sizeof(double))},
			   data.data());
	return b;
}

PYBIND11_MODULE(pnsolver, m)
{
	/*
	m.def( "solve_multigrid", &solve_multigrid);
	m.def( "solve_multigrid2", &solve_multigrid2);
	*/
	m.def( "solve_cg", &solve_cg);
	m.def( "solve_cg_eigen", &solve_cg_eigen);
	m.def( "solve_ls_cg", &solve_ls_cg);
	m.def( "solve_ls_cg_eigen", &solve_ls_cg_eigen);
	m.def( "solve_gs", &solve_gs);
	m.def( "solve_mg", &solve_mg);
	m.def( "solve_ls_lscg", &solve_ls_lscg);


	/*
	m.def( "solve_cg_eigen", &solve_cg_eigen);
	m.def( "solve_gs", &solve_gs);
	m.def( "solve_blockgs", &solve_blockgs);
	m.def( "solve_sparseLU", &solve_sparseLU);
	*/
	//m.def( "solve_lscg", &solve_lscg);

	/*
	m.def( "createComplexToRealConversionMatrix", &createComplexToRealConversionMatrix);
	m.def( "createRealToComplexConversionMatrix", &createRealToComplexConversionMatrix);
	*/
	m.def( "save_solution", &save_solution);
	m.def( "getCoefficientArray", &getCoefficientArray);

	m.def( "load_solution", &load_solution);
	m.def( "load_voxelgridfield3d", &load_voxelgridfield3d );
	//m.def( "getSolutionVector", &getSolutionVector);



	/*
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
	*/

	// MultigridSolver ================================
	py::class_<MultigridSolver> class_multigridsolver(m, "MultigridSolver");
	class_multigridsolver
	.def("__init__",
	[](MultigridSolver &m, int numLevels)
	{
		new (&m) MultigridSolver(numLevels);
	})
	.def("setMultigridLevel",
	[]( MultigridSolver &m, int lvl_index,
							const MultigridSolver::RealMatrix& A,
							const MultigridSolver::RealMatrix& downsample,
							const MultigridSolver::RealMatrix& upsample )
	{
		m.getLevel(lvl_index)->A = A;
		m.getLevel(lvl_index)->b = MultigridSolver::RealVector(A.cols());
		m.getLevel(lvl_index)->b.fill(0.0);
		m.getLevel(lvl_index)->x = MultigridSolver::RealVector(A.cols());
		m.getLevel(lvl_index)->x.fill(0.0);
		m.getLevel(lvl_index)->downsample = downsample;
		m.getLevel(lvl_index)->upsample = upsample;
	})
	.def("setb",
	[]( MultigridSolver &m, const MultigridSolver::RealVector& b )
	{
		m.getLevel(0)->b = b;
	})
	.def("solve", &MultigridSolver::solve)
	;



	// PNSolution ==============================
	py::class_<PNSolution, PNSolution::Ptr> class_pnsolution(m, "PNSolution");
	class_pnsolution
	.def("__init__",
	///*
	[](PNSolution &m,
		int order,
		const Eigen::Matrix<int, 3, 1>& resolution,
		const Eigen::Matrix<double, 3, 1>& bound_min,
		const Eigen::Matrix<double, 3, 1>& bound_max,
		const Eigen::VectorXd& data
		)
	{
		new (&m) PNSolution( order, resolution, Box3d(bound_min, bound_max), data.data() );
	})
	.def("eval",
	 [](PNSolution &m,
		const Eigen::Matrix<double, 3, 1>& pWS,
		const Eigen::Matrix<double, 3, 1>& direction )
	 {
		return m.eval(pWS, direction);
	 })
	.def("evalCoefficient",
	 [](PNSolution &m,
		const Eigen::Matrix<double, 3, 1>& pWS,
		int coeff_index)
	 {
		return m.evalCoefficient(pWS, coeff_index);
	 })
	.def("save", &PNSolution::save)
	.def("getResolution", &PNSolution::getResolution )
	.def("getNumCoeffs", &PNSolution::getNumCoeffs )
	.def("getBoundMin", &PNSolution::getBoundMin )
	.def("getBoundMax", &PNSolution::getBoundMax )
	.def("localToWorld",
	[]( PNSolution &m,
		const Eigen::Matrix<double, 3, 1>& pLS)
	{
		return m.localToWorld(pLS);
	})
	.def("getCoefficientField",
		[]( PNSolution &m, int coeff_index )
		{
			int numCoeffs = m.getNumCoeffs();
			V3i res = m.getResolution();
			std::vector<double> data(res[0]*res[1]*res[2]);

			for( int i=0;i<res[0];++i )
				for( int j=0;j<res[1];++j )
					for( int k=0;k<res[2];++k )
					{
						int index_dst = i*res[2]*res[1] + j*res[2] + k;
						int index_src = m.getIndex(V3i(i, j, k))+coeff_index;
						data[index_dst] = m.data()[index_src];
					}
			py::array b( //py::dtype(std::string("float64")),
						py::dtype::of<double>(),
					   {res[0], res[1], res[2]},
					   {int(sizeof(double))*res[2]*res[1], int(sizeof(double))*res[2], int(sizeof(double))},
					   data.data());
			return b;
		})
	;

	// PNSystem ==============================
	py::class_<PNSystem> class_pnsystem(m, "PNSystem");
	class_pnsystem
	.def("__init__",
	[](PNSystem &m, const std::string& stencil_name, PNVolume::Ptr problem, bool neumannBC)
	{
		PNSystem::Stencil stencil;

		if( stencil_name == "noop" )
			stencil = PNSystem::Stencil::noop();
		else
			stencil = PNSystem::findStencil(stencil_name);
		new (&m) PNSystem(stencil, problem, neumannBC);
	})
	.def("getNumCoefficients", &PNSystem::getNumCoefficients )
	.def("getResolution", &PNSystem::getResolution )
	.def("getNumVoxels", &PNSystem::getNumVoxels )
	.def("getOrder", &PNSystem::getOrder )
	.def("build", &PNSystem::build )
	//.def("setField", &PNSystem::setField )
	.def("get_A", &PNSystem::get_A )
	.def("get_b", &PNSystem::get_b )
	.def("get_A_real_test", &PNSystem::get_A_real_test )
	.def("get_b_real_test", &PNSystem::get_b_real_test )
	.def("setDebugVoxel", &PNSystem::setDebugVoxel )
	//.def("computeGroundtruth", &PNSystem::computeGroundtruth )
	.def("getVoxelInfo", &PNSystem::getVoxelInfo )
	.def("removeStaggering", &PNSystem::removeStaggering )
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
	.def("getSize",
		[](Domain &m)
		{
			return m.getSize();
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


	/*
	// Field ============================================================
	py::class_<Field, Field::Ptr> class_field(m, "Field");
	class_field
	.def("__call__",
	[](Field &m, const Eigen::Matrix<double, 3, 1>& pWS)
	{
		return m.eval(pWS);
	})
	;
	*/

	/*
	//VoxelGrid ============================================================
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
	*/

	// Field ============================================================

	// V3d ---
	py::class_<Field3d, Field3d::Ptr> class_field3d(m, "Field3d");
	class_field3d
	;

	// complex valued ---
	py::class_<Fieldcd, Fieldcd::Ptr> class_fieldcd(m, "Fieldcd");
	class_fieldcd
	.def("eval",
	 [](Fieldcd &m,
		const Eigen::Matrix<double, 3, 1>& pWS )
	 {
		return m.eval(pWS);
	 })
	;

	// ConstantField ============================================================

	// V3d ---
	py::class_<ConstantField3d, ConstantField3d::Ptr> class_ConstantField3d(m, "ConstantField3d", class_field3d);
	class_ConstantField3d
		.def("__init__",
		[](ConstantField3d &m, const Eigen::Vector3d& v)
		{
			new (&m) ConstantField3d(v);
		})
	;

	// complex ---

	// VoxelGridField ========================================================

	// V3d valued ---
	py::class_<VoxelGridField3d, VoxelGridField3d::Ptr> class_VoxelGridField3d(m, "VoxelGridField3d", class_field3d);
	class_VoxelGridField3d
	.def("__init__",
	[](VoxelGridField3d &m, py::array b)
	{
		typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> Strides;

		// Request a buffer descriptor from Python
		py::buffer_info info = b.request();

		// Some sanity checks ...
		if (info.format != py::format_descriptor<double>::format())
			throw std::runtime_error("Incompatible format: expected a double array!");

		if (info.ndim != 4)
			throw std::runtime_error("Incompatible buffer dimension!");

		if (info.shape[3] != 3)
			throw std::runtime_error("Incompatible buffer dimension!");

		V3i resolution( int(info.shape[0]),
						int(info.shape[1]),
						int(info.shape[2]) );

		auto data = static_cast<double*>(info.ptr);
		new (&m) VoxelGridField3d(resolution, (V3d*)data);
	})
	.def("save", &VoxelGridField3d::save)
	.def("asArray",
		[]( VoxelGridField3d &m )
		{
			V3i res = m.getResolution();
			py::array b( py::dtype::of<double>(),
						{res[0], res[1], res[2], 3},
						{int(sizeof(double))*res[2]*res[1]*3, int(sizeof(double))*res[2]*3, int(sizeof(double))*3, int(sizeof(double))},
						m.getData());
			return b;
		})
	;





	// PNVolume =====================================================
	py::class_<PNVolume, PNVolume::Ptr> class_PNVolume(m, "PNVolume");
	class_PNVolume
		.def("__init__",
		[](PNVolume &m, const Domain& domain)
		{
			new (&m) PNVolume(domain);
		})
		.def("setExtinctionAlbedo", &PNVolume::setExtinctionAlbedo)
		.def("setExtinctionMinimumThreshold", &PNVolume::setExtinctionMinimumThreshold)
		.def("setEmission",
			[]( PNVolume &pnvolume, int l, int m, Field3d::Ptr field )
			{
				pnvolume.setEmission(l, m, field);
			})
		.def("setPhase",
			[]( PNVolume &pnvolume, int l, int m, Field3d::Ptr field )
			{
				pnvolume.setPhase(l, m, field);
			})
	;

}
