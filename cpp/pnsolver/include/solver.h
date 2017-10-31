#pragma once


#include <iostream>
#include <map>

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

#include <common/Domain.h>
#include <field/Field.h>
#include <field/Constant.h>
#include <PNSystem.h>
#include <util/timer.h>

Eigen::VectorXd to_vector( const std::vector<double>& values );
void buildUpAndDownsamplingMatrices( PNSystem& sys_fine, PNSystem::RealMatrix& downsampleMatrix, PNSystem::RealMatrix& upsampleMatrix );


struct Solver
{

	struct System
	{
		System( int index = 0 ):
			coarse(0),
			index(index),
			numPreSmoothingIterations(1),
			numPostSmoothingIterations(1)
		{
		}

		int numPreSmoothingIterations;
		int numPostSmoothingIterations;

		Eigen::SparseMatrix<double, Eigen::RowMajor> A; // we have rowmajor matrix so that we can iterate over row elements
		Eigen::VectorXd      b;
		Eigen::VectorXd      x;
		Eigen::VectorXd      r;

		double computeRMSE()
		{
			return (b-A*x).norm();
		}


		System* coarse; // used for multigrid
		PNSystem::RealMatrix upsample; // converts x to the next coarser grid
		PNSystem::RealMatrix downsample; // converts from the next coarser grid to this grid

		// used for blockgs test
		PNSystem::VoxelManager m_voxelManager;

		int index;
	};

	std::vector<Eigen::LLT<Eigen::MatrixXd>> direct_solvers;
	Solver( int numLevels = 1 )
	{
		m_levels = std::vector<System>(numLevels);
		for( int i=0;i<m_levels.size()-1;++i )
		{
			m_levels[i].coarse = &m_levels[i+1];
			m_levels[i].index = 0;
		}
		m_levels.back().index = numLevels-1;


		for( int i=0;i<3;++i )
		{
			direct_solvers.push_back(Eigen::LLT<Eigen::MatrixXd>());
		}

	}

	void setMultigridLevel( int level, const Eigen::SparseMatrix<double, Eigen::RowMajor>& A, const Eigen::MatrixXd& u, const Eigen::MatrixXd& b, const PNSystem::RealMatrix& downsample, const PNSystem::RealMatrix& upsample )
	{
		m_levels[level].A = A;
		m_levels[level].x = u;
		m_levels[level].b = b;
		m_levels[level].downsample = downsample;
		m_levels[level].upsample = upsample;
	}

	// runs a number of Gauss-Seidel iterations on the given problem
	void run_gs_iterations( Eigen::SparseMatrix<double, Eigen::RowMajor> A, Eigen::VectorXd& b, Eigen::VectorXd& x, int numIterations )
	{
		for( int k=0;k<numIterations;++k )
		{
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
	}

	// runs a number of block-Gauss-Seidel iterations on the given problem
	void run_blockgs_iterations( PNSystem::VoxelManager& vm, Eigen::SparseMatrix<double, Eigen::RowMajor> A, Eigen::VectorXd& b, Eigen::VectorXd& x, int numIterations )
	{
		Eigen::LLT<Eigen::MatrixXd> llt;

		for( int k=0;k<numIterations;++k )
		{
			// iterate all voxels
			for( auto&voxel:vm.getVoxels() )
			{
				// assemble local matrix system
				int numCoeffs = vm.getNumCoeffs(voxel);

				if(!numCoeffs)
					continue;

				int globalOffset = voxel.globalOffset;
				//Eigen::MatrixXd A_voxel = A.block(globalOffset, globalOffset, globalOffset+numCoeffs, globalOffset+numCoeffs);
				Eigen::MatrixXd A_voxel(numCoeffs, numCoeffs);
				Eigen::VectorXd b_voxel(numCoeffs);// = b.segment(globalOffset, globalOffset+numCoeffs);

				//std::cout << b_voxel.rows() << " " << b_voxel.cols() << " globalOffset=" << globalOffset << " numCoeffs=" << numCoeffs << std::endl;

				// iterate all rows for the current voxel-block
				for( int i=0;i<numCoeffs;++i )
				{
					int globalIndex_i = i+globalOffset;
					double& b_voxel_coeff = b_voxel.coeffRef(i);
					b_voxel_coeff = b.coeffRef(globalIndex_i);
					// iterate all non-zero column elements (this is why we need RowMajor storage order for A)
					for( Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A, globalIndex_i);it;++it )
					{
						int globalIndex_j = it.col();

						if( (globalIndex_j<globalOffset)||(globalIndex_j>=globalOffset+numCoeffs) )
							b_voxel_coeff -= it.value()*x.coeffRef(globalIndex_j);
						else
							A_voxel.coeffRef( i, globalIndex_j - globalOffset ) = it.value();
					}
				}

				// solve
				llt.compute(A_voxel);
				Eigen::VectorXd x_voxel = llt.solve(b_voxel);
				for( int i=0;i<numCoeffs;++i )
					x.coeffRef(globalOffset+i) = x_voxel.coeffRef(i);
			}// for each voxel
		}
	}

	// runs a number of CG iterations on the given problem
	// returns the square root of the residual
	double run_cg_iterations( Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& b, Eigen::VectorXd& x, Eigen::VectorXd& r, int numIterations, double tol )
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

	std::string indent_str( int level )
	{
		std::string result = "";
		for( int i=0;i<level;++i )
			result+="\t";
		return result;
	}

	void multigrid_cycle( System* lvl_fine )
	{
		System* lvl_coarse = lvl_fine->coarse;

		// pre smoothing
		run_gs_iterations( lvl_fine->A, lvl_fine->b, lvl_fine->x, 1);

		// compute residual on fine level
		lvl_fine->r = lvl_fine->b - lvl_fine->A*lvl_fine->x;

		// restriction
		lvl_coarse->b = lvl_fine->downsample*lvl_fine->r;

		// compute approximate solution to the correction equation on the coarser grid
		if( lvl_coarse->coarse == 0 )
		{

			// the coarse level is the last level...fully solve the thing
			run_cg_iterations( lvl_coarse->A, lvl_coarse->b, lvl_coarse->x, lvl_coarse->r, 1000, 1.0e-10 );
			//run_gs_iterations( lvl_coarse->A, lvl_coarse->b, lvl_coarse->x, 1000, 1.0e-10 );
		}else
		{
			lvl_coarse->x.fill(0.0);
			for( int i=0;i<1;++i )
				multigrid_cycle( lvl_coarse );
		}

		// upsample correction and apply
		lvl_fine->x = lvl_fine->x + lvl_fine->upsample*lvl_coarse->x;

		// post smoothing
		//run_cg_iterations( lvl_fine->A, lvl_fine->b, lvl_fine->x, lvl_fine->r, 1, 0.0 );
		run_gs_iterations( lvl_fine->A, lvl_fine->b, lvl_fine->x, 1);
	}


	Eigen::VectorXd to_vector( const std::vector<double>& values )
	{
		Eigen::RowVectorXd x = Eigen::VectorXd( values.size() );
		for( int i=0;i<values.size();++i )
			x(i) = values[i];
		return x;
	}

	std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd>  solve(int maxIter)
	{
		std::cout << "solving mg\n";std::flush(std::cout);
		Timer timer;

		std::vector<double> solve_convergence;
		std::vector<double> solve_convergence_timestamps;

		timer.start();

		m_levels[0].x.fill(0.0);

		double tol = 1.0e-10; // convergence error tolerance
		for( int i=0;i<maxIter;++i )
		{
			double rmse = (m_levels[0].b-m_levels[0].A*m_levels[0].x).norm();
			//double rmse = (m_levels[0].x - u_ref).norm();
			solve_convergence.push_back(rmse);
			solve_convergence_timestamps.push_back(timer.elapsedSeconds());
			if(rmse<tol)
				break;
			multigrid_cycle(&m_levels[0]);
		}
		timer.stop();
		std::cout << "solve_multigrid: " << timer.elapsedSeconds() << "s #iterations=" << solve_convergence.size() << "\n";
		return std::make_tuple( m_levels[0].x,
		                        to_vector(solve_convergence),
		                        to_vector(solve_convergence_timestamps));

	}


	// runs a number of Gauss-Seidel iterations on the given problem
	void run_gs_iterations2( Eigen::SparseMatrix<double, Eigen::RowMajor> A, Eigen::VectorXd& b, Eigen::VectorXd& x, int numIterations = 100 )
	{
		for( int k=0;k<numIterations;++k )
		{
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

			//double rel_error = std::sqrt((b-A*x).norm()/b.norm());
			//if( rel_error < 1.0e-1 )
			//	return rel_error;
		}
	}

	void multigrid_cycle2( System* lvl_fine )
	{
		System* lvl_coarse = lvl_fine->coarse;

		// pre smoothing
		//run_gs_iterations2( lvl_fine->A, lvl_fine->b, lvl_fine->x, lvl_fine->numPreSmoothingIterations);
		run_blockgs_iterations( lvl_fine->m_voxelManager, lvl_fine->A, lvl_fine->b, lvl_fine->x, lvl_fine->numPreSmoothingIterations);

		if( lvl_fine->index == 0 )
		{
			//std::cout << "rel_error=" << rel_error << std::endl;
		}

		// compute residual on fine level
		lvl_fine->r = lvl_fine->b - lvl_fine->A*lvl_fine->x;

		// restriction
		lvl_coarse->b = lvl_fine->downsample*lvl_fine->r;

		// compute approximate solution to the correction equation on the coarser grid
		if( lvl_coarse->coarse == 0 )
		{

			// the coarse level is the last level...fully solve the thing
			run_cg_iterations( lvl_coarse->A, lvl_coarse->b, lvl_coarse->x, lvl_coarse->r, 1000, 1.0e-10 );
			//run_gs_iterations( lvl_coarse->A, lvl_coarse->b, lvl_coarse->x, 1000, 1.0e-10 );
		}else
		{
			lvl_coarse->x.fill(0.0);
			for( int i=0;i<1;++i )
				multigrid_cycle2( lvl_coarse );
		}

		// upsample correction and apply
		lvl_fine->x = lvl_fine->x + lvl_fine->upsample*lvl_coarse->x;

		// post smoothing
		//run_cg_iterations( lvl_fine->A, lvl_fine->b, lvl_fine->x, lvl_fine->r, 1, 0.0 );
		//run_gs_iterations2( lvl_fine->A, lvl_fine->b, lvl_fine->x, lvl_fine->numPostSmoothingIterations);
		run_blockgs_iterations( lvl_fine->m_voxelManager, lvl_fine->A, lvl_fine->b, lvl_fine->x, lvl_fine->numPostSmoothingIterations);
	}


	std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd>  solve2(int maxIter)
	{
		std::cout << "solving mg\n";std::flush(std::cout);
		Timer timer;

		std::vector<double> solve_convergence;
		std::vector<double> solve_convergence_timestamps;

		timer.start();

		m_levels[0].x.fill(0.0);

		double tol = 1.0e-10; // convergence error tolerance
		for( int i=0;i<maxIter;++i )
		{
			double rmse = (m_levels[0].b-m_levels[0].A*m_levels[0].x).norm();
			//double rmse = (m_levels[0].x - u_ref).norm();
			solve_convergence.push_back(rmse);
			solve_convergence_timestamps.push_back(timer.elapsedSeconds());
			if(rmse<tol)
				break;
			multigrid_cycle2(&m_levels[0]);
		}
		timer.stop();
		std::cout << "solve_multigrid: " << timer.elapsedSeconds() << "s #iterations=" << solve_convergence.size() << "\n";
		return std::make_tuple( m_levels[0].x,
								to_vector(solve_convergence),
								to_vector(solve_convergence_timestamps));

	}


	std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_gs()
	{
		std::cout << "solving gs\n";std::flush(std::cout);
		Timer timer;
		Eigen::SparseMatrix<double, Eigen::RowMajor> A = m_levels[0].A;
		Eigen::VectorXd b = m_levels[0].b;
		//Eigen::VectorXd x = b; // initial guess
		Eigen::VectorXd x = Eigen::VectorXd(b.rows()); // initial guess
		x.fill(0.0);

		std::vector<double> solve_convergence;
		std::vector<double> solve_convergence_timestamps;
		timer.start();


		int numIterations = 1000;
		double tol = 1.0e-10; // convergence error tolerance
		for( int k=0;k<numIterations;++k )
		{
			double rmse = (b - A*x).norm();
			solve_convergence.push_back(rmse);
			solve_convergence_timestamps.push_back(timer.elapsedSeconds());
			if( rmse < tol )
				break;

			run_gs_iterations(A, b, x, 1);
		}

		timer.stop();

		return std::make_tuple( m_levels[0].x,
		                        to_vector(solve_convergence),
		                        to_vector(solve_convergence_timestamps));
	}


	std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_cg()
	{
		std::cout << "solving cg\n";std::flush(std::cout);

		Timer timer;
		std::vector<double> solve_convergence;
		std::vector<double> solve_convergence_timestamps;


		PNSystem::RealMatrix A = m_levels[0].A;
		Eigen::VectorXd b = m_levels[0].b;
		//Eigen::VectorXd x = b; // initial guess
		Eigen::VectorXd x = Eigen::VectorXd(b.rows()); // initial guess
		x.fill(0.0);
		Eigen::VectorXd r( x.rows() );

		double tol = 1.0e-10; // convergence error tolerance
		int numIterations = b.rows();

		//solve_convergence.reserve(numIterations);
		//solve_convergence_timestamps.reserve(numIterations);


		timer.start();




		// cg solver
		{
			r = b-A*x;
			Eigen::VectorXd p = r;
			Eigen::VectorXd Ap;
			double rsold = r.squaredNorm();
			for( int i=0;i<numIterations;++i )
			{
				//solve_convergence.push_back(std::sqrt(rsold));
				//solve_convergence_timestamps.push_back(timer.elapsedSeconds());

				Ap = A*p;
				double alpha = rsold/p.dot(Ap);
				x = x + alpha*p;
				r = r - alpha*Ap;
				double rsnew = r.squaredNorm();
				if( std::sqrt(rsnew) < tol )
					break;
				p = r + (rsnew/rsold)*p;
				rsold = rsnew;
			}
		} // end of cg solver

		timer.stop();
		std::cout << "solve_cg: " << timer.elapsedSeconds() << "s #iterations=" << solve_convergence.size() << "\n";
		//return std::make_tuple( x,
		//                        to_vector(solve_convergence),
		//                        to_vector(solve_convergence_timestamps));
		return std::make_tuple( x,
								Eigen::VectorXd(),
								Eigen::VectorXd());

	}


	std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_cg_eigen()
	{
		std::cout << "solving cg eigen\n";std::flush(std::cout);

		std::vector<double> solve_convergence;
		std::vector<double> solve_convergence_timestamps;

		m_levels[0].x.fill(0.0);
		double rmse = (m_levels[0].b - m_levels[0].A*m_levels[0].x).norm();
		solve_convergence.push_back(rmse);
		solve_convergence_timestamps.push_back(0.0);


		Timer timer;
		timer.start();

		//Eigen::SparseLU<PNSystem::RealMatrix> solver;
		//typedef Eigen::ConjugateGradient<PNSystem::RealMatrix,Eigen::Lower, Eigen::IncompleteCholesky<double>> ICCG;
		//ICCG solver;

		Eigen::ConjugateGradient<PNSystem::RealMatrix> solver;
		solver.compute(m_levels[0].A);
		if(solver.info()!=Eigen::Success)
		{
			throw std::runtime_error("solve_cg_eigen decomposition failed");
		}
		m_levels[0].x = solver.solve(m_levels[0].b);
		if(solver.info()!=Eigen::Success)
		{
			throw std::runtime_error("solve_cg_eigen solve failed");
		}

		timer.stop();

		rmse = (m_levels[0].b - m_levels[0].A*m_levels[0].x).norm();
		solve_convergence.push_back(rmse);
		solve_convergence_timestamps.push_back(timer.elapsedSeconds());


		std::cout << "solve_cg_eigen: " << timer.elapsedSeconds() << "s #iterations=" <<  solver.iterations() << "\n";
		return std::make_tuple( m_levels[0].x,
								to_vector(solve_convergence),
								to_vector(solve_convergence_timestamps));
	}

	std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_sparseLU()
	{
		std::cout << "solving sparseLU(eigen)\n";std::flush(std::cout);

		std::vector<double> solve_convergence;
		std::vector<double> solve_convergence_timestamps;

		m_levels[0].x.fill(0.0);
		double rmse = (m_levels[0].b - m_levels[0].A*m_levels[0].x).norm();
		solve_convergence.push_back(rmse);
		solve_convergence_timestamps.push_back(0.0);


		Timer timer;
		timer.start();

		Eigen::SparseLU<PNSystem::RealMatrix> solver;
		solver.compute(m_levels[0].A);
		if(solver.info()!=Eigen::Success)
		{
			throw std::runtime_error("solve_cg_eigen decomposition failed");
		}
		m_levels[0].x = solver.solve(m_levels[0].b);
		if(solver.info()!=Eigen::Success)
		{
			throw std::runtime_error("solve_cg_eigen solve failed");
		}

		timer.stop();

		rmse = (m_levels[0].b - m_levels[0].A*m_levels[0].x).norm();
		solve_convergence.push_back(rmse);
		solve_convergence_timestamps.push_back(timer.elapsedSeconds());


		std::cout << "solve_sparseLU: " << timer.elapsedSeconds() << "\n";
		return std::make_tuple( m_levels[0].x,
								to_vector(solve_convergence),
								to_vector(solve_convergence_timestamps));
	}


	std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_blockgs()
	{
		std::cout << "solving blockgs\n";std::flush(std::cout);
		Timer timer;
		Eigen::SparseMatrix<double, Eigen::RowMajor> A = m_levels[0].A;
		Eigen::VectorXd b = m_levels[0].b;
		//Eigen::VectorXd x = b; // initial guess
		Eigen::VectorXd x = Eigen::VectorXd(b.rows()); // initial guess
		x.fill(0.0);

		std::vector<double> solve_convergence;
		std::vector<double> solve_convergence_timestamps;
		timer.start();


		int numIterations = 100;
		double tol = 1.0e-10; // convergence error tolerance
		for( int k=0;k<numIterations;++k )
		{
			double rmse = (b - A*x).norm();
			solve_convergence.push_back(rmse);
			solve_convergence_timestamps.push_back(timer.elapsedSeconds());
			if( rmse < tol )
				break;


			run_blockgs_iterations(m_levels[0].m_voxelManager, A, b, x, 1);
		}

		timer.stop();

		return std::make_tuple( m_levels[0].x,
								to_vector(solve_convergence),
								to_vector(solve_convergence_timestamps));
	}


	std::vector<System> m_levels;


};
