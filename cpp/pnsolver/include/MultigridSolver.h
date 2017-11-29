#pragma once
#include <iostream>
#include <util/timer.h>
#include <Eigen/Sparse>



#include <PNSystem.h>






struct MultigridSolver
{
	// we have rowmajor matrix so that we can iterate over row elements
	typedef Eigen::SparseMatrix<double, Eigen::RowMajor> RealMatrix;
	typedef Eigen::VectorXd RealVector;

	struct Level
	{
		Level( int index = 0 ):
		    coarse(0),
		    index(index),
		    numPreSmoothingIterations(1),
		    numPostSmoothingIterations(1)
		{
		}

		double computeRMSE()
		{
			//std::cout << A.rows() << " " << A.cols() << std::endl;
			//std::cout << x.rows() << " " << x.cols() << std::endl;
			//std::cout << b.rows() << " " << b.cols() << std::endl;
			return (b-A*x).norm();
		}

		int index;
		int numPreSmoothingIterations;
		int numPostSmoothingIterations;
		RealMatrix A;
		RealVector b;
		RealVector x;
		RealVector r;
		RealMatrix upsample; // converts x to the next coarser grid
		RealMatrix downsample; // converts from the next coarser grid to this grid
		Level* coarse; // link to the next coarser level


		// debugging
		PNSystem* sys_level;
	};

	MultigridSolver()
	{
	}

	MultigridSolver( int numLevels )
	{
		for( int i=0;i<numLevels;++i )
			addLevel();
	}

	~MultigridSolver();

	Level* addLevel();
	Level* getLevel(int index);
	int getNumLevels()const;

	// runs a number of Gauss-Seidel iterations on the given problem
	void run_gs_iterations( RealMatrix& A, RealVector& b, RealVector& x, int numIterations )
	{
		for( int k=0;k<numIterations;++k )
		{
			// iterate all rows
			for( int i=0;i<A.rows();++i )
			{
				double aii = 0.0;
				double sum = b(i);
				// iterate all non-zero column elements (this is why we need RowMajor storage order for A)
				for( RealMatrix::InnerIterator it(A, i);it;++it )
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

	// runs a number of CG iterations on the given problem
	// returns the square root of the residual
	double run_cg_iterations( RealMatrix& A, RealVector& b, RealVector& x, RealVector& r, int numIterations, double tol )
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

	void multigrid_cycle( Level* lvl_fine )
	{
		Level* lvl_coarse = lvl_fine->coarse;
		//std::cout << "multigrid_cycle lvl=" << lvl_fine->index << std::endl;

		// pre smoothing
		run_gs_iterations( lvl_fine->A, lvl_fine->b, lvl_fine->x, 1);

		// compute residual on fine level
		lvl_fine->r = lvl_fine->b - lvl_fine->A*lvl_fine->x;

		// restriction
		lvl_coarse->b = lvl_fine->downsample*lvl_fine->r;

		// compute approximate solution to the correction equation on the coarser grid
		if( lvl_coarse->coarse == 0 )
		{
			//std::cout << "multigrid_cycle solving coarse level" << std::endl;
			// the coarse level is the last level...fully solve the thing
			run_cg_iterations( lvl_coarse->A, lvl_coarse->b, lvl_coarse->x, lvl_coarse->r, 1000, 1.0e-14 );
			//run_gs_iterations( lvl_coarse->A, lvl_coarse->b, lvl_coarse->x, 1000, 1.0e-10 );
			/*
			Eigen::SparseLU<RealMatrix> solver;
			solver.compute(lvl_coarse->A);
			if(solver.info()!=Eigen::Success)
			{
				throw std::runtime_error("multigrid_cycle decomposition failed");
			}
			lvl_coarse->x = solver.solve(lvl_coarse->b);
			if(solver.info()!=Eigen::Success)
			{
				throw std::runtime_error("multigrid_cycle solve failed");
			}
			*/
			//std::cout << "multigrid_cycle solving coarse level...done" << std::endl;
		}else
		{
			// recurse to the next coarser level
			lvl_coarse->x.fill(0.0);
			for( int i=0;i<1;++i )
				multigrid_cycle( lvl_coarse );
		}

		// upsample correction and apply it
		lvl_fine->x = lvl_fine->x + lvl_fine->upsample*lvl_coarse->x;

		// post smoothing
		run_gs_iterations( lvl_fine->A, lvl_fine->b, lvl_fine->x, 1);
	}



	std::tuple<RealVector, RealVector, RealVector> solve(int maxIter)
	{
		std::cout << "Multigridsolver::solve\n";std::flush(std::cout);
		Timer timer;

		std::vector<double> list_rmse;
		std::vector<double> list_time;

		timer.start();

		m_levels[0]->x.fill(0.0);

		///*
		double tol = 1.0e-10; // convergence error tolerance
		int iteration = 0;

		for( iteration=0;iteration<maxIter;++iteration )
		{
			double rmse = m_levels[0]->computeRMSE();
			//std::cout << "rmse=" << rmse << std::endl;
			list_rmse.push_back(rmse);
			list_time.push_back(timer.elapsedSeconds());
			if(rmse<tol)
				break;
			multigrid_cycle(m_levels[0]);
		}
		//*/

		timer.stop();

		std::cout << "Multigridsolver::solve: " << timer.elapsedSeconds() << "s #iterations=" << iteration << " rmse=" << m_levels[0]->computeRMSE() << "\n";

		return std::make_tuple( m_levels[0]->x,
		                        to_vector(list_rmse),
		                        to_vector(list_time));
	}


	Eigen::VectorXd to_vector( const std::vector<double>& values )
	{
		Eigen::RowVectorXd x = Eigen::VectorXd( values.size() );
		for( int i=0;i<values.size();++i )
			x(i) = values[i];
		return x;
	}
private:
	std::vector<Level*> m_levels;
};


void buildUpAndDownsamplingMatrices2( PNSystem& sys_fine, MultigridSolver::RealMatrix& downsampleMatrix, MultigridSolver::RealMatrix& upsampleMatrix );
