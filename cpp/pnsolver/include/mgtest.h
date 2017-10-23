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




struct MGTEST
{

	struct MultigridLevel
	{
		MultigridLevel( int index = 0 ):
			next(0),
			index(index)
		{

		}

		size_t getMemory()
		{
			size_t result = 0;
			result+=b.cols()*b.rows()*sizeof(double);
			result+=x.cols()*x.rows()*sizeof(double);
			result+=r.cols()*r.rows()*sizeof(double);
			result+=A.data().size();
			return result;
		}

		Eigen::SparseMatrix<double, Eigen::RowMajor> A;
		Eigen::VectorXd      b;
		Eigen::VectorXd      x;
		Eigen::VectorXd      r;

		PNSystem::RealMatrix upsample; // converts x to the next coarser grid
		PNSystem::RealMatrix downsample; // converts from the next coarser grid to this grid

		MultigridLevel* next;
		int index;
	};

	Eigen::VectorXd u_ref;

	MGTEST( int numLevels = 1 )
	{
		m_levels = std::vector<MultigridLevel>(numLevels);
		for( int i=0;i<m_levels.size()-1;++i )
		{
			//domain = domain.downsample();
			//fields = fields.createRestricted();
			//PNSystem sys_coarse(stencil, domain, boundaryConditions);
			//sys_coarse.setFields( fields );
			//sys_coarse.build();

			//levels[i+1].A = sys_coarse.get_A_real().transpose()*sys_coarse.get_A_real();
			//levels[i+1].b = sys_coarse.get_A_real().transpose()*sys_coarse.get_b_real();
			//levels[i+1].x = Eigen::VectorXd(sys_coarse.getVoxelManager().getNumUnknowns());
			//levels[i+1].x.fill(0.0);
			//levels[i+1].r = Eigen::VectorXd(sys_coarse.getVoxelManager().getNumUnknowns());
			//buildUpAndDownsamplingMatrices(sys_coarse, levels[i+1].downsample, levels[i+1].upsample);

			m_levels[i].next = &m_levels[i+1];
			m_levels[i].index = i;
		}

	}

	void setMultigridLevel( int level, const Eigen::SparseMatrix<double, Eigen::RowMajor>& A, const Eigen::MatrixXd& u, const PNSystem::RealMatrix& downsample, const PNSystem::RealMatrix& upsample )
	{
		m_levels[level].A = A;
		m_levels[level].x = u;
		m_levels[level].downsample = downsample;
		m_levels[level].upsample = upsample;
	}



	void memoryTest()
	{
		int level = 0;
		size_t sum = 0;
		for( auto&lvl:m_levels )
		{
			size_t mem = lvl.getMemory();
			sum += mem;
			std::cout << "level " << level << " :" << mem << std::endl;
			++level;
		}

		std::cout << "sum=" << sum << std::endl;
	}

	void setRef( const Eigen::MatrixXd& u_ref2 )
	{
		u_ref = u_ref2;
	}

	void setb( const Eigen::VectorXd& b )
	{
		m_levels[0].b = b;
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

	void multigrid_cycle( MultigridLevel* lvl_fine )
	{
		MultigridLevel* lvl_coarse = lvl_fine->next;

		// pre smoothing
		run_gs_iterations( lvl_fine->A, lvl_fine->b, lvl_fine->x, 1);

		// compute residual on fine level
		lvl_fine->r = lvl_fine->b - lvl_fine->A*lvl_fine->x;

		// restriction
		lvl_coarse->b = lvl_fine->downsample*lvl_fine->r;

		// compute approximate solution to the correction equation on the coarser grid
		if( lvl_coarse->next == 0 )
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
		timer.start();

		std::vector<double> solve_convergence;
		std::vector<double> solve_convergence_timestamps;



		PNSystem::RealMatrix A = m_levels[0].A;
		Eigen::VectorXd b = m_levels[0].b;
		//Eigen::VectorXd x = b; // initial guess
		Eigen::VectorXd x = Eigen::VectorXd(b.rows()); // initial guess
		x.fill(0.0);
		Eigen::VectorXd r( x.rows() );




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
				solve_convergence.push_back(std::sqrt(rsold));
				solve_convergence_timestamps.push_back(timer.elapsedSeconds());

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
		return std::make_tuple( x,
		                        to_vector(solve_convergence),
		                        to_vector(solve_convergence_timestamps));
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

	//std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_multigrid(PNSystem& sys, int numLevels)
	//{
	    /*
		ProfileTimer timers;
		double tol = 1.0e-10; // convergence error tolerance
		//int numLevels = 9; // for res=512
		//int numLevels = 6; // for res=64
		Timer timer;
		std::vector<double> solve_convergence;
		std::vector<double> solve_convergence_timestamps;
		std::vector<MultigridLevel> levels(numLevels);

		levels[0].A = sys.get_A_real().transpose()*sys.get_A_real();
		levels[0].b = sys.get_A_real().transpose()*sys.get_b_real();
		levels[0].x = levels[0].b;
		levels[0].r = levels[0].b-levels[0].A*levels[0].x;
		buildUpAndDownsamplingMatrices(sys, levels[0].downsample, levels[0].upsample);


		PNSystem::Stencil& stencil = sys.getStencil();
		int boundaryConditions = sys.getBoundaryConditions();
		Domain domain = sys.getDomain();
		PNSystem::Fields fields = sys.getFields();

		// link up
		for( int i=0;i<levels.size()-1;++i )
		{
			domain = domain.downsample();
			fields = fields.createRestricted();
			PNSystem sys_coarse(stencil, domain, boundaryConditions);
			sys_coarse.setFields( fields );
			sys_coarse.build();

			levels[i+1].A = sys_coarse.get_A_real().transpose()*sys_coarse.get_A_real();
			levels[i+1].b = sys_coarse.get_A_real().transpose()*sys_coarse.get_b_real();
			levels[i+1].x = Eigen::VectorXd(sys_coarse.getVoxelManager().getNumUnknowns());
			levels[i+1].x.fill(0.0);
			levels[i+1].r = Eigen::VectorXd(sys_coarse.getVoxelManager().getNumUnknowns());
			buildUpAndDownsamplingMatrices(sys_coarse, levels[i+1].downsample, levels[i+1].upsample);

			levels[i].next = &levels[i+1];
		}

		timer.start();
		int maxIter = 1000;
		for( int i=0;i<maxIter;++i )
		{
			multigrid_cycle( &levels[0], &timers );
			double rmse = (levels[0].b-levels[0].A*levels[0].x).norm();
			solve_convergence.push_back( rmse );
			solve_convergence_timestamps.push_back(timer.elapsedSeconds());
			if(rmse<tol)
				break;
		}

		timer.stop();
		std::cout << "solve_multigrid_test: " << timer.elapsedSeconds() << "s #iterations=" << solve_convergence.size() << "\n";
		timers.print();
		return std::make_tuple( sys.stripBoundary(levels[0].x),
								to_vector(solve_convergence),
								to_vector(solve_convergence_timestamps));
		*/
	//}


	std::vector<MultigridLevel> m_levels;


};
