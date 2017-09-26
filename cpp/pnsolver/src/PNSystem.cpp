#include <PNSystem.h>
#include <field/Function.h>
#include <field/VoxelGrid.h> // used only for getting maximum value of sigma_t

#include<Eigen/IterativeLinearSolvers>
#include<Eigen/SparseLU>
//#include<Eigen/SparseCholesky>

#include <math/rng.h>
#include <util/threadpool.h>

std::map<std::string, PNSystem::Stencil> PNSystem::g_stencils;

PNSystem::PNSystem( Stencil stencil, const Domain& domain)
	:
	m_domain(domain),
	m_fields(stencil.order),
	m_stencil(stencil)
{
	m_numCoeffs = 0;
	for( int l=0;l<=m_stencil.order;++l )
		for( int m=-l;m<=l;++m )
			// in 2d, we only need to solve for moments where l+m is even
			if( (l+m) % 2 == 0 )
			{
				m_lm_to_index[std::make_pair(l,m)] = m_numCoeffs;
				m_index_to_lm[m_numCoeffs] = std::make_pair(l,m);
				++m_numCoeffs;
			}

	/*
	// temp: boundary test ---
	m_numBoundaryVoxels = (domain.resolution()[0]+4)*(domain.resolution()[1]+4) - domain.numVoxels();
	m_boundary_A_complex = ComplexMatrix((domain.numVoxels()+m_numBoundaryVoxels)*m_numCoeffs, (domain.numVoxels()+m_numBoundaryVoxels)*m_numCoeffs);
	m_boundary_b_complex = ComplexVector((domain.numVoxels()+m_numBoundaryVoxels)*m_numCoeffs, 1);
	m_boundary_A_real = RealMatrix((domain.numVoxels()+m_numBoundaryVoxels)*m_numCoeffs, (domain.numVoxels()+m_numBoundaryVoxels)*m_numCoeffs);
	m_boundary_b_real = RealVector((domain.numVoxels()+m_numBoundaryVoxels)*m_numCoeffs, 1);

	//std::map<std::pair<int,int>, int> voxel_to_global_boundary_index;

	V2i current(-1,-1);
	V2i step(0, 1);

	//V2i res(4,4);
	//int numBoundaryVoxels = (res[0]+4)*(res[1]+4) - res[0]*res[1];
	V2i res = domain.resolution();
	int numBoundaryVoxels = m_numBoundaryVoxels;
	int layer = 0;
	int count = 0;
	int sidecount = 0;
	int direction = 0;


	for( int i=0;i<numBoundaryVoxels;++i,++count )
	{
		m_voxel_to_global_boundary_index[std::make_pair(current[0], current[1])] = i;

		//std::cout << "index=" << i << " current=" << current[0] << " " << current[1] <<  " count=" << count << std::endl;

		if( count == res[direction] + layer*2 + 1 )
		{
			// rotate clockwise
			step = V2i(step[1], -step[0]);
			direction = !direction;
			count = 0;
			++sidecount;
		}

		if( (count == res[direction] + layer*2)&&
			(sidecount==3))
		{
			count = -1;
			sidecount = 0;
			++layer;
			current = V2i( -layer - 1, -layer - 1 );
			step = V2i(step [1], -step [0]);
		}else
			current += step;
	}
	// end of temp
	*/


	debugVoxel = V2i(-1, -1);

	// TEMP: testing with non-rasterized rte parameters
	Function::FunctionType sigma_a = [](const P2d& pWS)
	{
		double x = pWS[0];
		double y = pWS[1];
		double cx = std::ceil(x);
		double cy = std::ceil(y);
		double g = 0;
		if( (std::ceil((x+y)/2.0)*2.0 == (cx+cy)) &&
			(cx > 1.0) && (cx < 7.0) &&
			(cy > 1.0) &&
			(cy-2.0*std::abs(cx-4.0) < 4) )
			g = 1;
		return (1.0-g)*0 + g*10;
	};
	Function::FunctionType sigma_s = [](const P2d& pWS)
	{
		double x = pWS[0];
		double y = pWS[1];
		double cx = std::ceil(x);
		double cy = std::ceil(y);
		double g = 0;
		if( (std::ceil((x+y)/2.0)*2.0 == (cx+cy)) &&
			(cx > 1.0) && (cx < 7.0) &&
			(cy > 1.0) &&
			(cy-2.0*std::abs(cx-4.0) < 4) )
			g = 1;
		return (1.0-g)*1 + g*0;
	};
	Function::FunctionType sigma_t = [=](const P2d& pWS)
	{
		return sigma_a(pWS)+sigma_s(pWS);
	};
	Function::FunctionType phase_shcoeff = [](const P2d& pWS)
	{
		return 1.0;
	};
	Function::FunctionType source_shcoeffs = [](const P2d& pWS)
	{
		double x = pWS[0];
		double y = pWS[1];
		if( (x > 3.0) && (x < 4.0) && (y > 3.0) && (y < 4.0) )
			return 1.0;
		return 0.0;

	};

	/*
	m_fields.sigma_a = std::make_shared<Function>(sigma_a);
	m_fields.sigma_s = std::make_shared<Function>(sigma_s);
	m_fields.sigma_t = std::make_shared<Function>(sigma_t);
	m_fields.f_p->setField(0,0,std::make_shared<Function>(phase_shcoeff));
	m_fields.q->setField(0,0,std::make_shared<Function>(source_shcoeffs));
	*/
}



struct ParticleTracer2D : public Task
{
	Eigen::MatrixXd accumulation_buffer;
	PNSystem::Fields* fields;
	const Domain* domain;
	RNGd rng;
	int numSamples;

	ParticleTracer2D()
	{
	}

	virtual void run() override
	{
		for( int i=0;i<numSamples;++i )
			sample();
	}

	void sample()
	{
		// particle state
		P2d pWS;
		V2d d;
		//double pdf_over_throughput = 1.0;


		// create initial position

		// randomly sample emission field to find the next spawning position
		// assuming checkerboard emission field
		pWS = P2d( rng.next1D() + 3.0, rng.next1D() + 3.0 );
		// currently we sample a unit square, for which the pdf is 1.0
		//pdf_over_throughput /= 1.0;

		//P2d pV;
		//pWS = m_domain.voxelToWorld(pVS);
		// TODO: envmap sampling

		// create initial direction (assuming isotropic emission)
		{
			V3d d_sphere = sampleSphere<double>(rng);
			d = V2d(d_sphere[0], d_sphere[1]);
			//pdf_over_throughput /= sampleSpherePDF();
		}


		// light tracing
		//int depth = 0;
		//while( depth++ == 0 )
		while(true)
		{
			// PROPAGATE ==============================
			// woodcock/delta tracking
			// TODO: assuming VoxelGrid field
			double sigma_t_max = std::dynamic_pointer_cast<VoxelGrid>(fields->sigma_t)->getMax();

			double t = 0.0;
			int event = -1; // 0 = absortion, 1=scattering

			//int numSteps = 0;
			//while( numSteps++ == 0 )
			while(true)
			{
				double step = -log( 1.0-rng.next1D() )/sigma_t_max;
				t += step;

				pWS += t*d;

				if( !domain->getBound().contains(pWS) )
					break;

				double sigma_t = fields->sigma_t->eval(pWS).real();

				// russian roulette
				double dice = rng.next1D();

				if(dice<sigma_t/sigma_t_max)
				{
					// scattering or absorbtion
					double sigma_a = fields->sigma_a->eval(pWS).real();
					if(dice < sigma_a/sigma_t_max)
						event = 0;
					else
						event = 1;
					break;
				}
				else
				{
					// virtual particle: accumulate
					accumulate(pWS);
				}
			}



			if(event==-1)
				// propagation step has moved particle out of domain-we stop tracing this particle
				break;

			if(event==0)
				// particle has been absorbed
				break;

			// ACCUMULATE ==============================
			accumulate(pWS);


			///*
			// SCATTERING ==============================
			// create new direction (assuming isotropic scattering)
			{
				V3d d_sphere = sampleSphere<double>(rng);
				d = V2d(d_sphere[0], d_sphere[1]);
				//pdf_over_throughput /= sampleSpherePDF();
			}
			//*/
		} // path tracing while loop
	}

	void accumulate( const P2d& pWS )
	{
		P2d pVS = domain->worldToVoxel(pWS);
		P2i voxel;
		voxel = P2i( int(pVS[0]), int(pVS[1]) );

		//accumlation_buffer
		if( domain->contains_voxel(voxel) )
		{
			accumulation_buffer.coeffRef( voxel[0], voxel[1] ) += 1.0;
		}
	}
};


Eigen::MatrixXd PNSystem::computeGroundtruth(int numSamples)
{
	int numTasks = ThreadPool::instance()->getNumWorkers();
	std::cout << "numTasks=" << numTasks << std::endl;

	std::vector<ParticleTracer2D*> tasks;
	for( int i=0;i<numTasks;++i )
	{
		ParticleTracer2D* pt = new ParticleTracer2D();
		pt->accumulation_buffer = Eigen::MatrixXd::Zero(m_domain.resolution()[0], m_domain.resolution()[1]);
		pt->fields = &m_fields;
		pt->domain = &m_domain;
		pt->rng = RNGd(123+i);
		pt->numSamples = numSamples/numTasks;

		tasks.push_back(pt);
	}
	//tasks[0]->run();
	ThreadPool::instance()->enqueueAndWait(tasks);

	Eigen::MatrixXd result = Eigen::MatrixXd::Zero(m_domain.resolution()[0], m_domain.resolution()[1]);
	int count = 0;

	for( auto t:tasks )
	{
		result += t->accumulation_buffer;
		count += t->numSamples;
		delete t;
	}

	return result/double(count);
}

void PNSystem::setField( const std::string& id, Field::Ptr field )
{
	//std::cout << "WARNING: PNSystem::setField currently ignored. Using explicit RTE functions.\n";
	///*
	if( id == "sigma_t" )
		m_fields.sigma_t = field;
	else
	if( id == "sigma_a" )
		m_fields.sigma_a = field;
	else
	if( id == "sigma_s" )
		m_fields.sigma_s = field;
	else
	if( id == "f_p" )
		// we assume the zero coefficient is to be set
		m_fields.f_p->setField(0,0,field);
	else
	if( id == "q" )
		// we assume the zero coefficient is to be set
		m_fields.q->setField(0,0,field);
	else
		throw std::runtime_error("PNSystem::setField unable to set field " + id);
	//*/
}


PNSystem::VoxelSystem PNSystem::getVoxelSystem( const V2i& voxel )
{
	return PNSystem::VoxelSystem(this, voxel);
}

const Domain& PNSystem::getDomain()const
{
	return m_domain;
}

/*
PNSystem::ComplexMatrix& PNSystem::get_A_complex()
{
	return m_A_complex;
}

PNSystem::ComplexVector& PNSystem::get_b_complex()
{
	return m_b_complex;
}
*/

PNSystem::RealMatrix& PNSystem::get_A_real()
{
	return m_builder_A.matrix;
}

PNSystem::RealVector& PNSystem::get_b_real()
{
	return m_builder_b.matrix;
}

void PNSystem::setDebugVoxel(const V2i &dv)
{
	debugVoxel = dv;

}

PNSystem::RealMatrix &PNSystem::get_boundary_A_real()
{
	return m_boundary_A_real;
}

PNSystem::RealMatrix &PNSystem::get_boundary_b_real()
{
	return m_boundary_b_real;
}


int PNSystem::getGlobalIndex( V2i voxel, int coeff )const
{
	int voxel_index = voxel[1]*m_domain.resolution()[0] + voxel[0];
	return voxel_index*m_numCoeffs + coeff;
}

void PNSystem::getVoxelAndCoefficient( int global_index, V2i& voxel, int& coeff )const
{
	std::div_t divresult;
	divresult = std::div( global_index, m_numCoeffs );
	coeff = divresult.rem;


	divresult = div( divresult.quot, m_domain.resolution()[0] );
	voxel.x() = divresult.rem;
	voxel.y() = divresult.quot;
}


int PNSystem::getIndex(int l, int m) const
{
	return m_lm_to_index.at(std::make_pair(l, m));
}

void PNSystem::getLM(int sh_index, int &l, int &m) const
{
	auto it = m_index_to_lm.find(sh_index);
	if(it != m_index_to_lm.end())
	{
		l = it->second.first;
		m = it->second.second;
		return;
	}
	l=-1;
	m=-1;
}

int PNSystem::getNumCoefficients()const
{
	return m_numCoeffs;
}

int PNSystem::getNumVoxels()const
{
	return m_domain.numVoxels();
}

int PNSystem::getOrder() const
{
	return m_stencil.order;
}

int PNSystem::getGlobalBoundaryIndex( const V2i boundaryVoxel, int coeff )
{
	auto it = m_voxel_to_global_boundary_index.find(std::make_pair(boundaryVoxel[0], boundaryVoxel[1]));
	if( it != m_voxel_to_global_boundary_index.end() )
		return (m_domain.numVoxels()+it->second)*m_numCoeffs + coeff;
	std::cout << "boundaryVoxel=" << boundaryVoxel[0] << " " << boundaryVoxel[1] << std::endl;
	throw std::runtime_error("out of bound access of boundary voxels");
	return -1;
}

/*
PNSystem::MatrixAccessHelper PNSystem::A( V2i voxel_i,
										  int coefficient_i,
										  V2i voxel_j,
										  int coefficient_j )
{
	// boundary test ----
	bool boundary_i = !m_domain.contains_voxel(voxel_i);
	bool boundary_j = !m_domain.contains_voxel(voxel_j);

	if(boundary_i)
		throw std::runtime_error("voxel_i is boundary which is unexpected");

	if(boundary_j)
	{
		MatrixAccessHelper mah(&m_boundary_triplets_A);
		mah.m_global_i = getGlobalIndex(voxel_i, coefficient_i);
		mah.m_global_j = getGlobalBoundaryIndex(voxel_j, coefficient_j);
		return mah;
	}



	MatrixAccessHelper mah(&m_triplets_A);
	mah.m_global_i = getGlobalIndex(voxel_i, coefficient_i);
	mah.m_global_j = getGlobalIndex(voxel_j, coefficient_j);
	return mah;
}
*/


void PNSystem::build()
{
	std::cout << "PNSystem::build building system matrices A and b...\n";
	int min_x = 0;
	int min_y = 0;
	int max_x = m_domain.resolution()[0];
	int max_y = m_domain.resolution()[1];

	if( (debugVoxel.x() != -1)&&(debugVoxel.y() != -1) )
	{
		min_x = debugVoxel.x();
		max_x = min_x+1;
		min_y = debugVoxel.y();
		max_y = min_y+1;
	}

	// iterate all voxels and apply stencil which has been generated from python script ---
	for( int i=min_x;i<max_x;++i )
		for( int j=min_y;j<max_y;++j )
			m_stencil.apply( getVoxelSystem(V2i(i,j)), m_fields );

	m_builder_A.build(m_domain.numVoxels()*m_numCoeffs, m_domain.numVoxels()*m_numCoeffs);
	m_builder_b.build(m_domain.numVoxels()*m_numCoeffs, 1);
}

PNSystem::RealVector PNSystem::solve()
{
	std::cout << "PNSystem::solve solving for x...\n";

	//Eigen::ConjugateGradient<RealMatrix> solver;
	//Eigen::BiCGSTAB<RealMatrix> solver;
	Eigen::SparseLU<RealMatrix> solver;
	//solver.setMaxIterations(100000);
	std::cout << "PNSystem::solve compute...\n";std::flush(std::cout);
	solver.compute(get_A_real());
	if(solver.info()!=Eigen::Success)
	{
		throw std::runtime_error("PNSystem::solve decomposition failed");
	}
	std::cout << "PNSystem::solve solve...\n";std::flush(std::cout);
	RealVector x = solver.solve(get_b_real());
	//std::cout << "iterations=" << solver.iterations() << std::endl;
	if(solver.info()!=Eigen::Success)
	{
		throw std::runtime_error("PNSystem::solve solve failed");
	}

	return x;

}

PNSystem::RealVector PNSystem::solveWithGuess(PNSystem::RealVector &x0)
{
	std::cout << "PNSystem::solve solving with guess for x...\n";

	//Eigen::ConjugateGradient<RealMatrix> solver;
	Eigen::BiCGSTAB<RealMatrix> solver;
	//Eigen::SparseLU<RealMatrix> solver;
	//solver.setMaxIterations(100000);
	std::cout << "PNSystem::solve compute...\n";std::flush(std::cout);
	solver.compute(get_A_real());
	if(solver.info()!=Eigen::Success)
	{
		throw std::runtime_error("PNSystem::solve decomposition failed");
	}
	std::cout << "PNSystem::solve solve...\n";std::flush(std::cout);
	Eigen::VectorXd b = get_b_real();
	Eigen::VectorXd res = x0;
	//RealVector x = solver.solveWithGuess(get_b_real(), res);
	Eigen::VectorXd x = solver.solveWithGuess(b, res);
	//RealVector x = solver.solve(get_b_real());

	//std::cout << "iterations=" << solver.iterations() << std::endl;
	if(solver.info()!=Eigen::Success)
	{
		throw std::runtime_error("PNSystem::solve solve failed");
	}

	return x.sparseView();
}

/*
PNSystem::ComplexVector PNSystem::solve2()
{
	std::cout << "PNSystem::solve solving for x...\n";

	//Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver;
	//Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
	Eigen::SparseLU<Eigen::SparseMatrix<Complex> > solver;
	//solver.setMaxIterations(100000);
	solver.compute(get_A_complex());
	if(solver.info()!=Eigen::Success)
	{
		throw std::runtime_error("PNSystem::solve decomposition failed");
	}
	ComplexVector x = solver.solve(get_b_complex());
	//std::cout << "iterations=" << solver.iterations() << std::endl;
	if(solver.info()!=Eigen::Success)
	{
		throw std::runtime_error("PNSystem::solve solve failed");
	}

	return x;
}
*/

PNSystem::RealVector PNSystem::solve_boundary()
{
	std::cout << "PNSystem::solve solving for x...\n";

	// temp: boundary test:
	//Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver;
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
	solver.compute(get_boundary_A_real());
	if(solver.info()!=Eigen::Success)
	{
		throw std::runtime_error("PNSystem::solve decomposition failed");
	}
	RealVector x = solver.solve(get_boundary_b_real());
	if(solver.info()!=Eigen::Success)
	{
		throw std::runtime_error("PNSystem::solve solve failed");
	}

	return x;
}

//PNSystem::VoxelSystem ==============================================
PNSystem::VoxelSystem::VoxelSystem(PNSystem* pns,
				const V2i& voxel_i ):
	m_pns(pns),
	m_voxel_i(voxel_i)
{
}

const V2i& PNSystem::VoxelSystem::getVoxel()const
{
	return m_voxel_i;
}

P2d PNSystem::VoxelSystem::getVoxelSize() const
{
	return m_pns->getDomain().voxelSize();
}


PNSystem::MatrixBuilderd::MatrixAccessHelper PNSystem::VoxelSystem::coeff_A(int coefficient_i, V2i voxel_j, int coefficient_j)
{
	int global_i = m_pns->getGlobalIndex(m_voxel_i, coefficient_i);
	int global_j = m_pns->getGlobalIndex(voxel_j, coefficient_j);

	if( m_pns->getDomain().contains_voxel(m_voxel_i) &&
		m_pns->getDomain().contains_voxel(voxel_j))
		return m_pns->m_builder_A.coeff(global_i, global_j);

	// dirichlet BC
	return PNSystem::MatrixBuilderd::MatrixAccessHelper(0);
}

PNSystem::MatrixBuilderd::MatrixAccessHelper PNSystem::VoxelSystem::coeff_b(int coefficient_i)
{
	int global_i = m_pns->getGlobalIndex(m_voxel_i, coefficient_i);
	return m_pns->m_builder_b.coeff(global_i, 0);
}


V2d PNSystem::VoxelSystem::voxelToWorld(const V2d& pVS )const
{
	return m_pns->getDomain().voxelToWorld(pVS);
}

