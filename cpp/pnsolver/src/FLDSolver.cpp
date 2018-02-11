#include <FLDSolver.h>
#include <PNSystem.h>

#include <util/voxelgrid.h>
#include <util/threadpool.h>

int getVoxelIndex( const V3i& res, int i, int j, int k )
{
	return i*res[2]*res[1] + j*res[2] + k;
}

void getVoxelCoord( const V3i& res, int index, int& i, int& j, int& k )
{
	div_t divresult;

	divresult = div( index, res.z()*res.y() );

	i = divresult.quot;
	divresult = div( divresult.rem, res.z() );
	j = divresult.quot;
	k = divresult.rem;
}

struct GlobalInfo
{
	GlobalInfo():
		problem(0)
	{

	}

	int redblack;
	//int numVoxels;

	//VoxelGrid<double> phi;
	mutable VoxelGrid<double> D;
	mutable Eigen::VectorXd x;
	mutable Eigen::VectorXd r;
	PNVolume* problem;
};

struct TaskInfo
{
	TaskInfo( int taskid, int numTasks ):
		taskid(taskid),
		numTasks(numTasks)
	{

	}

	int taskid;
	int numTasks;
};

void gs_iteration(TaskInfo& ti, const GlobalInfo* gi)
{
	int color_channel = 0;
	Domain& domain = gi->problem->getDomain();
	V3i res = domain.getResolution();
	V3d vs = domain.getVoxelSize();
	V3d h_inv( 1.0/vs[0], 1.0/vs[1], 1.0/vs[2] );
	V3d h_inv2( 1.0/(vs[0]*vs[0]), 1.0/(vs[1]*vs[1]), 1.0/(vs[2]*vs[2]) );

	// iterate all rows
	for( int row=ti.taskid;row<domain.numVoxels();row+=ti.numTasks )
	{
		int i, j, k;
		getVoxelCoord(res, row, i, j, k);

		int redblack =(k%2+j%2+i%2)%2;
		if( redblack != gi->redblack )
			continue;

		V3d voxel_center_VS = V3d(i+0.5, j+0.5, k+0.5);
		P3d voxel_center_WS = domain.voxelToWorld(voxel_center_VS);

		// sample environment values -------
		double D_xp, D_xm, D_yp, D_ym, D_zp, D_zm;
		double phi_xp, phi_xm, phi_yp, phi_ym, phi_zp, phi_zm;
		D_xp = D_xm = D_yp = D_ym = D_zp = D_zm = 0.0;
		phi_xp = phi_xm = phi_yp = phi_ym = phi_zp = phi_zm = 0.0;

		if( i < res[0]-1 )
		{
			D_xp = gi->D.sample(i+1, j, k);
			phi_xp = gi->x.coeffRef(getVoxelIndex(res, i+1, j, k));
		}
		if( i > 0 )
		{
			D_xm = gi->D.sample(i-1, j, k);
			phi_xm = gi->x.coeffRef(getVoxelIndex(res, i-1, j, k));
		}
		if( j < res[1]-1 )
		{
			D_yp = gi->D.sample(i, j+1, k);
			phi_yp = gi->x.coeffRef(getVoxelIndex(res, i, j+1, k));
		}
		if( j > 0 )
		{
			D_ym = gi->D.sample(i, j-1, k);
			phi_ym = gi->x.coeffRef(getVoxelIndex(res, i, j-1, k));
		}
		if( k < res[2]-1 )
		{
			D_zp = gi->D.sample(i, j, k+1);
			phi_zp = gi->x.coeffRef(getVoxelIndex(res, i, j, k+1));
		}
		if( k > 0 )
		{
			D_zm = gi->D.sample(i, j, k-1);
			phi_zm = gi->x.coeffRef(getVoxelIndex(res, i, j, k-1));
		}

		// compute(and update) center diffusion coefficient -------
		double extinction_c = gi->problem->evalExtinction(voxel_center_WS)[color_channel];
		double F;
		//F = 1.0/3.0; // flux limiter for CDA
		///*
		// LP flux-limiter
		double phi_c = gi->x.coeffRef(row);
		V3d gradPhi;
		gradPhi[0] = (phi_xp - phi_xm)/(2.0*vs[0]);
		gradPhi[1] = (phi_yp - phi_ym)/(2.0*vs[1]);
		gradPhi[2] = (phi_zp - phi_zm)/(2.0*vs[2]);
		const double R = std::max(gradPhi.norm(), 1.0e-4) / std::max(phi_c*extinction_c, 1.0e-4);
		if(R<1.0e-4)
		{
			double Rsqr = R*R;
			F = 1.0f/3.0f - Rsqr/45.0 + 2.0*Rsqr*Rsqr/945.0;
		}else
			F = (1.0/tanh(R) - 1.0f/R)/R;
		//*/

		//double D_c = gi->D.sample(i, j, k);
		double D_c = F/extinction_c; // equ. 33 in FLD paper
		gi->D.lvalue(i, j, k) = D_c;


		double D_xph = 0.5*(D_c+D_xp);
		double D_xmh = 0.5*(D_c+D_xm);
		double D_yph = 0.5*(D_c+D_yp);
		double D_ymh = 0.5*(D_c+D_ym);
		double D_zph = 0.5*(D_c+D_zp);
		double D_zmh = 0.5*(D_c+D_zm);

		double denominator = 0.0;
		denominator += -h_inv2[0]*(D_xph+D_xmh);
		denominator += -h_inv2[1]*(D_yph+D_ymh);
		denominator += -h_inv2[2]*(D_zph+D_zmh);
		denominator += -gi->problem->evalAbsorption(voxel_center_WS)[color_channel];

		double numerator = -gi->problem->evalEmission(0, 0, voxel_center_WS)[color_channel];
		if( i < res[0]-1 )
			numerator -= h_inv2[0]*D_xph*gi->x.coeffRef(getVoxelIndex(res,i+1,j,k));
		if( i > 0 )
			numerator -= h_inv2[0]*D_xmh*gi->x.coeffRef(getVoxelIndex(res,i-1,j,k));
		if( j < res[1]-1 )
			numerator -= h_inv2[1]*D_yph*gi->x.coeffRef(getVoxelIndex(res,i,j+1,k));
		if( j > 0 )
			numerator -= h_inv2[1]*D_ymh*gi->x.coeffRef(getVoxelIndex(res,i,j-1,k));
		if( k < res[2]-1 )
			numerator -= h_inv2[2]*D_zph*gi->x.coeffRef(getVoxelIndex(res,i,j,k+1));
		if( k > 0 )
			numerator -= h_inv2[2]*D_zmh*gi->x.coeffRef(getVoxelIndex(res,i,j,k-1));

		double residual = numerator - gi->x.coeffRef(row)*denominator;
		gi->r.coeffRef(row) = residual;

		gi->x.coeffRef(row) = numerator/denominator;
	}
}


std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> FLDSolver::solve(PNVolume& problem, double tol, int maxIterations)
{
	Domain& domain = problem.getDomain();
	V3i res = domain.getResolution();
	V3d vs = domain.getVoxelSize();
	V3d h_inv( 1.0/vs[0], 1.0/vs[1], 1.0/vs[2] );
	V3d h_inv2( 1.0/(vs[0]*vs[0]), 1.0/(vs[1]*vs[1]), 1.0/(vs[2]*vs[2]) );

	int color_channel = 0;

	//Eigen::VectorXd x = Eigen::VectorXd(domain.numVoxels()); // initial guess
	//x.fill(0.0);


	GlobalInfo gi;
	gi.problem = &problem;
	gi.x = Eigen::VectorXd(domain.numVoxels()); // initial guess
	gi.x.fill(0.0);
	gi.r = Eigen::VectorXd(domain.numVoxels()); // initial guess
	gi.r.fill(0.0);
	gi.D.resize( res );

	// we precompute D for CDA
	/*
	for( int row=0;row<domain.numVoxels();++row )
	{
		int i, j, k;
		getVoxelCoord(res, row, i, j, k);
		P3d voxel_center_WS = domain.voxelToWorld(V3d(i+0.5, j+0.5, k+0.5));
		gi.D.lvalue(i, j, k) = 1.0/(3.0*problem.evalExtinction(voxel_center_WS)[color_channel]);
	}
	*/

	Timer timer;


	int numThreads = ThreadPool::getNumSystemCores();
	std::vector<GenericTask<TaskInfo, GlobalInfo>*> tasks;
	for (int i = 0; i < numThreads; i++)
		tasks.push_back(new GenericTask<TaskInfo, GlobalInfo>(TaskInfo(i, numThreads), &gi, gs_iteration));


	std::cout << "solving...\n";std::flush(std::cout);

	int iteration = 0;

	timer.start();

	double rmse = -1.0;
	for( iteration=0;iteration<maxIterations;++iteration )
	{
		//TaskInfo ti(0,1);
		gi.redblack = 0;
		//gs_iteration(ti, &gi);
		ThreadPool::instance()->enqueueAndWait(tasks);
		gi.redblack = 1;
		//gs_iteration(ti, &gi);
		ThreadPool::instance()->enqueueAndWait(tasks);
		rmse = gi.r.norm()/std::sqrt(problem.getDomain().numVoxels());

		std::cout << "iteration=" << iteration << " rmse=" << rmse << std::endl;
		if( rmse < tol )
			break;
	}

	timer.stop();
	std::cout << "solve_gs: " << timer.elapsedSeconds() << "s #iterations=" << iteration << " rmse=" << rmse <<  "\n";

	return std::make_tuple( gi.x,
							Eigen::VectorXd(),
							Eigen::VectorXd());
}

/*




void gs_iteration( TaskInfo& ti, const GlobalInfo* gi )
{
	// for each interior voxel
	for( int v = 0;v<gi->numVoxels;++v )
	{
		// work out index


		bool debug = false;
		//if ( i==29 && j==53 && k==16 )
		//	debug = true;
		// CDA CODE ------------------------------
		int color_channel = channel;
		V3i vi(i-1, j-1, k-1);
		V3d vd = vi.cast<double>();
		const Domain& domain = volume.getDomain();
		const PNVolume& problem = volume;
		V3d vs = domain.getVoxelSize();
		V3d h_inv( 1.0/vs[0], 1.0/vs[1], 1.0/vs[2] );
		V3d h_inv2( 1.0/(vs[0]*vs[0]), 1.0/(vs[1]*vs[1]), 1.0/(vs[2]*vs[2]) );

		V3d voxel_center_VS = vd + V3d(0.5, 0.5, 0.5);
		P3d voxel_center_WS = domain.voxelToWorld(voxel_center_VS);


		double phi_xp = phi.sample(i+1, j, k);
		double phi_xm = phi.sample(i-1, j, k);
		double phi_yp = phi.sample(j, j+1, k);
		double phi_ym = phi.sample(j, j-1, k);
		double phi_zp = phi.sample(i, j, k+1);
		double phi_zm = phi.sample(i, j, k-1);

		double D_xph = 1.0/(3.0*problem.evalExtinction(voxel_center_WS+V3d(vs[0]*0.5, 0.0, 0.0))[color_channel]);
		double D_xmh = 1.0/(3.0*problem.evalExtinction(voxel_center_WS-V3d(vs[0]*0.5, 0.0, 0.0))[color_channel]);
		double D_yph = 1.0/(3.0*problem.evalExtinction(voxel_center_WS+V3d(0.0, vs[1]*0.5, 0.0))[color_channel]);
		double D_ymh = 1.0/(3.0*problem.evalExtinction(voxel_center_WS-V3d(0.0, vs[1]*0.5, 0.0))[color_channel]);
		double D_zph = 1.0/(3.0*problem.evalExtinction(voxel_center_WS+V3d(0.0, 0.0, vs[2]*0.5))[color_channel]);
		double D_zmh = 1.0/(3.0*problem.evalExtinction(voxel_center_WS-V3d(0.0, 0.0, vs[2]*0.5))[color_channel]);

		// --

		double numerator = 0.0;
		numerator += h_inv2[0]*D_xph*phi_xp;
		numerator += h_inv2[0]*D_xmh*phi_xm;
		numerator += h_inv2[1]*D_yph*phi_yp;
		numerator += h_inv2[1]*D_ymh*phi_ym;
		numerator += h_inv2[2]*D_zph*phi_zp;
		numerator += h_inv2[2]*D_zmh*phi_zm;
		numerator += problem.evalEmission(0,0,voxel_center_WS)[color_channel];

		double denominator = 0.0;
		denominator += h_inv2[0]*D_xph;
		denominator += h_inv2[0]*D_xmh;
		denominator += h_inv2[1]*D_yph;
		denominator += h_inv2[1]*D_ymh;
		denominator += h_inv2[2]*D_zph;
		denominator += h_inv2[2]*D_zmh;
		numerator += problem.evalAbsorption(voxel_center_WS)[color_channel];

		if( debug )
		{
			std::cout << "numerator=" << numerator << std::endl;
			std::cout << "denominator=" << denominator << std::endl;
			std::cout << "D " << D_xph << std::endl;
			std::cout << "D " << D_xmh << std::endl;
		}
		//phi.lvalue(i, j, k) = problem.evalEmission(0,0,voxel_center_WS)[color_channel];

		double residual = numerator - phi.lvalue(i, j, k)*denominator;
		if(residual>max_residual)
		{
			max_residual = residual;
			max_residual_voxel = V3i(i, j, k);
		}
		rmse += (residual*residual-rmse)/double(voxel_index+1);

		// update phi value
		phi.lvalue(i, j, k) = numerator/denominator;
	} // for each voxel
}
*/

std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> FLDSolver::solve2(PNVolume& volume, double tol, int maxIterations)
{
	V3i res = volume.getDomain().getResolution();


	/*
	GlobalInfo gi;
	gi.volume = &volume;
	gi.phi.resize(res[0]+2, res[1]+2, res[2]+2);
	gi.D.resize(res[0]+2, res[1]+2, res[2]+2);

	int numThreads = ThreadPool::getNumSystemCores();

	std::vector<GenericTask<TaskInfo, GlobalInfo>*> tasks;
	for (int i = 0; i < numThreads; i++)
		tasks.push_back(new GenericTask<TaskInfo, GlobalInfo>(TaskInfo(i, numThreads), &gi, gs_iteration));

	//for( term.reset();term.keepRunning();term.advance() )
	{
		//std::cout << term.m_sample <<" " << term.m_numSamples <<  std::endl;


//		gi.redblack = 0;
//		ThreadPool::instance()->enqueueAndWait(tasks);
//		gi.redblack = 1;
//		ThreadPool::instance()->enqueueAndWait(tasks);

		// single threaded run
		TaskInfo ss_ti(0, 1);
		gs_iteration(&ss_ti, &gi);

		// compute residuum and check for tolerance...
	}

	for( auto task:tasks )
		delete task;
	*/

	VoxelGrid<double> phi;
	VoxelGrid<double> D;



	// holds the diffusion coefficient for each voxel
	// this is updated in-place for each iteration
	//Eigen::VectorXd D( volume.getDomain().numVoxels() );
	//D.fill(0.0);




	// including layer of boundary voxel
	phi.resize(res[0]+2, res[1]+2, res[2]+2);
	D.resize(res[0]+2, res[1]+2, res[2]+2);

	// here we assume all channels are equal
	// extending FLD to rgb channels is straightforward...
	int channel = 0;

	int numIterations = maxIterations;

	// for each iteration
	for( int it=0;it<numIterations;++it )
	{
		int voxel_index = 0;
		double rmse = 0.0;
		double max_residual = 0.0;
		V3i max_residual_voxel(0,0,0);
		// for each interior voxel
		for( int k=1;k<res.z()+1;++k )
			for( int j=1;j<res.y()+1;++j )
				for( int i=1;i<res.x()+1;++i, ++voxel_index )
				{
					bool debug = false;
					//if ( i==29 && j==53 && k==16 )
					//	debug = true;
					// CDA CODE ------------------------------
					int color_channel = channel;
					V3i vi(i-1, j-1, k-1);
					V3d vd = vi.cast<double>();
					const Domain& domain = volume.getDomain();
					const PNVolume& problem = volume;
					V3d vs = domain.getVoxelSize();
					V3d h_inv( 1.0/vs[0], 1.0/vs[1], 1.0/vs[2] );
					V3d h_inv2( 1.0/(vs[0]*vs[0]), 1.0/(vs[1]*vs[1]), 1.0/(vs[2]*vs[2]) );

					V3d voxel_center_VS = vd + V3d(0.5, 0.5, 0.5);
					P3d voxel_center_WS = domain.voxelToWorld(voxel_center_VS);


					double phi_xp = phi.sample(i+1, j, k);
					double phi_xm = phi.sample(i-1, j, k);
					double phi_yp = phi.sample(j, j+1, k);
					double phi_ym = phi.sample(j, j-1, k);
					double phi_zp = phi.sample(i, j, k+1);
					double phi_zm = phi.sample(i, j, k-1);

					double D_xph = 1.0/(3.0*problem.evalExtinction(voxel_center_WS+V3d(vs[0]*0.5, 0.0, 0.0))[color_channel]);
					double D_xmh = 1.0/(3.0*problem.evalExtinction(voxel_center_WS-V3d(vs[0]*0.5, 0.0, 0.0))[color_channel]);
					double D_yph = 1.0/(3.0*problem.evalExtinction(voxel_center_WS+V3d(0.0, vs[1]*0.5, 0.0))[color_channel]);
					double D_ymh = 1.0/(3.0*problem.evalExtinction(voxel_center_WS-V3d(0.0, vs[1]*0.5, 0.0))[color_channel]);
					double D_zph = 1.0/(3.0*problem.evalExtinction(voxel_center_WS+V3d(0.0, 0.0, vs[2]*0.5))[color_channel]);
					double D_zmh = 1.0/(3.0*problem.evalExtinction(voxel_center_WS-V3d(0.0, 0.0, vs[2]*0.5))[color_channel]);

					// --

					double numerator = 0.0;
					numerator += h_inv2[0]*D_xph*phi_xp;
					numerator += h_inv2[0]*D_xmh*phi_xm;
					numerator += h_inv2[1]*D_yph*phi_yp;
					numerator += h_inv2[1]*D_ymh*phi_ym;
					numerator += h_inv2[2]*D_zph*phi_zp;
					numerator += h_inv2[2]*D_zmh*phi_zm;
					numerator += problem.evalEmission(0,0,voxel_center_WS)[color_channel];

					double denominator = 0.0;
					denominator += h_inv2[0]*D_xph;
					denominator += h_inv2[0]*D_xmh;
					denominator += h_inv2[1]*D_yph;
					denominator += h_inv2[1]*D_ymh;
					denominator += h_inv2[2]*D_zph;
					denominator += h_inv2[2]*D_zmh;
					numerator += problem.evalAbsorption(voxel_center_WS)[color_channel];

					if( debug )
					{
						std::cout << "numerator=" << numerator << std::endl;
						std::cout << "denominator=" << denominator << std::endl;
						std::cout << "D " << D_xph << std::endl;
						std::cout << "D " << D_xmh << std::endl;
					}
					//phi.lvalue(i, j, k) = problem.evalEmission(0,0,voxel_center_WS)[color_channel];

					double residual = numerator - phi.lvalue(i, j, k)*denominator;
					if(residual>max_residual)
					{
						max_residual = residual;
						max_residual_voxel = V3i(i, j, k);
					}
					rmse += (residual*residual-rmse)/double(voxel_index+1);

					// update phi value
					phi.lvalue(i, j, k) = numerator/denominator;



					//double phi_new =

					/*
					// Assembling global system =============
					double aii = 0.0;
					aii += -h_inv2[0]*(D_xph+D_xmh);
					aii += -h_inv2[1]*(D_yph+D_ymh);
					aii += -h_inv2[2]*(D_zph+D_zmh);
					aii += -problem.evalAbsorption(voxel_center_WS)[color_channel];
					ctx.coeff_A( 0, vi + V3i(0,0,0), 0 ) += aii;

					ctx.coeff_A( 0, vi + V3i(1,0,0), 0 ) += h_inv2[0]*D_xph;
					ctx.coeff_A( 0, vi + V3i(-1,0,0), 0 ) += h_inv2[0]*D_xmh;
					ctx.coeff_A( 0, vi + V3i(0,1,0), 0 ) += h_inv2[1]*D_yph;
					ctx.coeff_A( 0, vi + V3i(0,-1,0), 0 ) += h_inv2[1]*D_ymh;
					ctx.coeff_A( 0, vi + V3i(0,0,1), 0 ) += h_inv2[2]*D_zph;
					ctx.coeff_A( 0, vi + V3i(0,0,-1), 0 ) += h_inv2[2]*D_zmh;

					ctx.coeff_b( 0 ) += -problem.evalEmission(0, 0, voxel_center_WS)[color_channel];
					*/





























					// FLD CODE ----------------------------------

					/*
					// single component index
					int idx = k*res.x*res.y + j*res.x + i;
					// read kt
					double kt = volume->evalExtinction(pWS)[channel];

					// vector3 component index
					idx = idx*3+c;
					int idx_xp = idx + 3;
					int idx_xm = idx - 3;
					int idx_yp = idx + res.x*3;
					int idx_ym = idx - res.x*3;
					int idx_zp = idx + res.x*res.y*3;
					int idx_zm = idx - res.x*res.y*3;

					float phi = _phi[idx];

					// diffusion coefficient at voxelcenter
					double dc(0.0f);

					// phi stencil
					double phi_xp = _phi[idx_xp];
					double phi_xm = _phi[idx_xm];
					double phi_yp = _phi[idx_yp];
					double phi_ym = _phi[idx_ym];
					double phi_zp = _phi[idx_zp];
					double phi_zm = _phi[idx_zm];

					// numerator
					float numerator = h*_j[idx];
					// denominator
					float denom = h * (1.0 - albedo) * kt;

					// knudsen number ----------
					cumath::V3f gradPhi;
					gradPhi.x = (phi_xp - phi_xm)/(2.0f*h);
					gradPhi.y = (phi_yp - phi_ym)/(2.0f*h);
					gradPhi.z = (phi_zp - phi_zm)/(2.0f*h);
					const float R = max(gradPhi.getLength(), epsilon_j) / max(phi*kt, epsilon_j);

					// flux limiter -------
					double F;

					// CDA
					//F = 1.0f/3.0f;
					// Larsen
					//F = std::pow(std::pow(3,dN) + std::pow(R,dN), -1.0/dN);
					// Kershaw
					//F = 2.0f/(3.0f + sqrt(9.0f + 4.0f*R*R));
					// Wilson Max
					//F = 1.0/max(3.0, R);
					// Wilson
					//F = 1.0/(3.0 + R);
					// Grosjean
					//F = (2.0f-albedo)/3.0f;




					// LP ---
					if(R<1.0e-4)
					{
						double Rsqr = R*R;
						F = 1.0f/3.0f - Rsqr/45.0f + 2.0f*Rsqr*Rsqr/945.0f;
					}else
						F = (1.0f/tanh(R) - 1.0f/R)/R;


					dc = hInv*0.5f*(1.0f/kt)*F;


					// we update diffusion coefficient in place
					_D[idx] = dc;

					float D_xph = dc + _D[idx_xp];
					float D_xmh = dc + _D[idx_xm];
					float D_yph = dc + _D[idx_yp];
					float D_ymh = dc + _D[idx_ym];
					float D_zph = dc + _D[idx_zp];
					float D_zmh = dc + _D[idx_zm];

					// numerator ---
					numerator += D_xph * phi_xp;
					numerator += D_xmh * phi_xm;
					numerator += D_yph * phi_yp;
					numerator += D_ymh * phi_ym;
					numerator += D_zph * phi_zp;
					numerator += D_zmh * phi_zm;

					// denominator ---
					denom += D_xph; // (dimensionless) terms coming from diffusion operator
					denom += D_xmh;
					denom += D_yph;
					denom += D_ymh;
					denom += D_zph;
					denom += D_zmh;

					//if(_residual)
					//	_residual[idx] = (numerator - denom*phi)/h;

					// standard
					_phi[idx] = numerator / denom;
					// sor
					//_phi[idx] = sor*(numerator / denom) + (1.0-sor)*phi;
					*/
				}
		std::cout << "iteration=" << it << " rmse=" << rmse << std::endl;
		if(rmse < tol)
			break;
		//std::cout << "max_residual=" << max_residual << " @" << max_residual_voxel[0] << ", " << max_residual_voxel[1] << ", " << max_residual_voxel[2] << std::endl;
	} // iteration

	// extract solution vector from phi data
	// this is the solution vector (single fluence value for each voxel)
	Eigen::VectorXd u( volume.getDomain().numVoxels() );
	for( int k=0;k<res.z();++k )
		for( int j=0;j<res.y();++j )
			for( int i=0;i<res.x();++i )
			{
				int index_dst = i*res[2]*res[1] + j*res[2] + k;
				u.coeffRef(index_dst) = phi.sample(i+1, j+1, k+1);
			}

	return std::make_tuple(u, u, u);
}
