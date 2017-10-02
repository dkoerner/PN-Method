#include <PNSystem.h>
#include <field/Function.h>
#include <field/VoxelGrid.h> // used only for getting maximum value of sigma_t

#include<Eigen/IterativeLinearSolvers>
#include<Eigen/SparseLU>
//#include<Eigen/SparseCholesky>

#include <math/rng.h>
#include <util/threadpool.h>

std::map<std::string, PNSystem::Stencil> PNSystem::g_stencils;

PNSystem::PNSystem( Stencil stencil, const Domain& domain, bool neumannBC)
	:
	m_domain(domain),
	m_fields(stencil.order),
	m_stencil(stencil),
	m_neumannBC(neumannBC),
	m_voxelManager()
{
	// find out number of coefficients and mapping between sh bands l,m and linear indices ---
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

	// voxelmanager needs to be initialized after the number of coefficients has been determined
	m_voxelManager.init(this);


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

PNSystem::MatrixBuilderd::MatrixAccessHelper PNSystem::coeff_A(V2i voxel_i, int coefficient_i, V2i voxel_j, int coefficient_j)
{
	if( m_voxelManager.voxelIsValid(voxel_i) && m_voxelManager.voxelIsValid(voxel_j) )
	{
		int globalIndex_i = m_voxelManager.getGlobalIndex(m_voxelManager.getVoxel(voxel_i), coefficient_i);
		int globalIndex_j = m_voxelManager.getGlobalIndex(m_voxelManager.getVoxel(voxel_j), coefficient_j);

		if( (globalIndex_i >= 0)&&(globalIndex_j >= 0) )
			return m_builder_A.coeff(globalIndex_i, globalIndex_j);
	}

	return PNSystem::MatrixBuilderd::MatrixAccessHelper();
}

PNSystem::MatrixBuilderd::MatrixAccessHelper PNSystem::coeff_b(V2i voxel_i, int coefficient_i)
{
	int globalIndex_i = m_voxelManager.getGlobalIndex(m_voxelManager.getVoxel(voxel_i), coefficient_i);

	if( globalIndex_i >= 0 )
		return m_builder_b.coeff(globalIndex_i, 0);
	return PNSystem::MatrixBuilderd::MatrixAccessHelper();
}


PNSystem::VoxelManager& PNSystem::getVoxelManager()
{
	return m_voxelManager;
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

PNSystem::Fields &PNSystem::getFields()
{
	return m_fields;
}

const Domain& PNSystem::getDomain()const
{
	return m_domain;
}

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
	//for( int i=min_x;i<max_x;++i )
	//	for( int j=min_y;j<max_y;++j )
	//		m_stencil.apply( getVoxelSystem(V2i(i,j)), m_fields );
	for( auto& v: m_voxelManager.getVoxels() )
		m_stencil.apply(Stencil::Context(*this, v));

	m_builder_A.build(m_voxelManager.getNumUnknowns(), m_voxelManager.getNumUnknowns());
	m_builder_b.build(m_voxelManager.getNumUnknowns(), 1);
}



PNSystem::RealMatrix PNSystem::getVoxelInfo(const std::string &info)
{
	MatrixBuilderd indexMatrix;
/*
	V2i res = m_domain.resolution();
	int numLayers = m_numBoundaryLayers;
	for( int i=-numLayers;i<res[0]+numLayers;++i )
		for( int j=-numLayers;j<res[1]+numLayers;++j )
		{
			V2i coord(i, j);
			Voxel& v = getVoxel2(coord);
			double value = 0.123;

			if( info == "globalOffset" )
				value = v.globalOffset;
			else
			if( info == "type" )
				value = double(int(v.type));
			else
			if( info == "coord.x" )
				value = double(int(v.coord[0]));
			else
			if( info == "coord.y" )
				value = double(int(v.coord[1]));
			else
			if( info == "tmp0" )
				value = double(int(v.tmp0));
			else
			if( info == "tmp1" )
				value = double(int(v.tmp1));

			indexMatrix.coeff(i+numLayers,j+numLayers) += value;
		}

	indexMatrix.build(res[0]+numLayers*2, res[1]+numLayers*2);
*/
	//std::cout << indexMatrix.matrix.cols() << " " << indexMatrix.matrix.rows() << std::endl;
	return indexMatrix.matrix;
}

PNSystem::RealVector PNSystem::solve()
{
	std::cout << "PNSystem::solve solving for x...\n";

	//Eigen::ConjugateGradient<RealMatrix> solver;
	//Eigen::BiCGSTAB<RealMatrix> solver;
	Eigen::SparseLU<RealMatrix> solver;
	//solver.setMaxIterations(1000);
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



	// ..................
	RealVector x2( m_domain.numVoxels()*m_numCoeffs,1);
	std::vector<RealTriplet> triplets;
	std::vector<Voxel>& voxels = m_voxelManager.getVoxels();
	for( auto&v:voxels )
	{
		if( !m_domain.contains_voxel(v.coord) )
			// boundary or mixed boundary voxel
			continue;
		for( int i=0;i<m_numCoeffs;++i )
		{
			int new_index = getGlobalIndex(v.coord, i);
			int index = m_voxelManager.getGlobalIndex(v, i);
			if(index==-1)
				throw std::runtime_error("PNSystem::solve didnt expect invalid coefficient index");
			triplets.push_back(RealTriplet(new_index, 0, x.coeffRef(index, 0)));
		}
	}
	x2.setFromTriplets(triplets.begin(), triplets.end());

	return x2;
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

int PNSystem::getBoundaryConditions()
{
	return int(m_neumannBC);
}

PNSystem::Stencil& PNSystem::getStencil()
{
	return m_stencil;
}


// Stencil::Context --------------------------------------------

PNSystem::MatrixBuilderd::MatrixAccessHelper PNSystem::Stencil::Context::coeff_A(int coeff_i, const V2i &voxel_j, int coeff_j)
{
	// check if row for coefficient i is actually supposed to be set by the stencil
	// This is not the case with corner boundary voxels, where the equation for some coefficients are
	// descibing the average with neighbouring coefficients. In this case we prevent the stencil from
	// writing into those rows
	if(sys.m_voxelManager.isBoundaryCoefficient(voxel, coeff_i))
		return PNSystem::MatrixBuilderd::MatrixAccessHelper();

	return sys.coeff_A(voxel.coord, coeff_i, voxel_j, coeff_j);
}

PNSystem::MatrixBuilderd::MatrixAccessHelper PNSystem::Stencil::Context::coeff_b(int coeff_i)
{
	return sys.coeff_b(voxel.coord, coeff_i);
}


// VoxelManager -------------------------------------------------

PNSystem::VoxelManager::VoxelManager()
{
}

void PNSystem::VoxelManager::init(PNSystem *sys)
{
	this->sys = sys;
	m_resolution = sys->getDomain().resolution();
	m_numBoundaryLayers = sys->getStencil().width;


	mixedTypes[ std::make_pair(1, 0) ] = Voxel::EVT_N;
	mixedTypes[ std::make_pair(0, 1) ] = Voxel::EVT_E;
	mixedTypes[ std::make_pair(-1, 0) ] = Voxel::EVT_S;
	mixedTypes[ std::make_pair(0, -1) ] = Voxel::EVT_W;
	mixedTypes[ std::make_pair(1, 1) ] = Voxel::EVT_NE;
	mixedTypes[ std::make_pair(-1, 1) ] = Voxel::EVT_SE;
	mixedTypes[ std::make_pair(-1, -1) ] = Voxel::EVT_SW;
	mixedTypes[ std::make_pair(1, -1) ] = Voxel::EVT_NW;

	// initialize all coefficients to boundary coefficients
	int numCoeffs = sys->getNumCoefficients();
	m_localCoefficientIndices = std::vector<std::vector<int>>( 10, std::vector<int>(numCoeffs, -1) );
	m_isBoundaryCoefficient = std::vector<std::vector<bool>>( 10, std::vector<bool>(numCoeffs, true) );
	m_numNonBoundaryCoefficients = std::vector<int>(10, 0);

	for( int i=0;i<numCoeffs;++i )
	{
		V2i offset = sys->m_stencil.getOffset(i);

		// interior cells have all voxels defines, therefore their indices are all set
		m_localCoefficientIndices[Voxel::EVT_Interior][i] = m_numNonBoundaryCoefficients[Voxel::EVT_Interior]++;
		m_isBoundaryCoefficient[Voxel::EVT_Interior][i] = false;

		// here we make those coefficients non-boundary coefficients, which sit directly on the boundary
		// I think this is required for both, neumann and dirichlet BC
		if( (offset[0] == 0)&&(offset[1] == 0) )
		{
			m_localCoefficientIndices[Voxel::EVT_N][i] = m_numNonBoundaryCoefficients[Voxel::EVT_N]++;
			m_isBoundaryCoefficient[Voxel::EVT_N][i] = false;
			m_localCoefficientIndices[Voxel::EVT_NE][i] = m_numNonBoundaryCoefficients[Voxel::EVT_NE]++;
			m_isBoundaryCoefficient[Voxel::EVT_NE][i] = false;
			m_localCoefficientIndices[Voxel::EVT_E][i] = m_numNonBoundaryCoefficients[Voxel::EVT_E]++;
			m_isBoundaryCoefficient[Voxel::EVT_E][i] = false;
		}else
		if( (offset[0] == 0)&&(offset[1] == 1) )
		{
			m_localCoefficientIndices[Voxel::EVT_N][i] = m_numNonBoundaryCoefficients[Voxel::EVT_N]++;
			m_isBoundaryCoefficient[Voxel::EVT_N][i] = false;
		}else
		if( (offset[0] == 1)&&(offset[1] == 0) )
		{
			m_localCoefficientIndices[Voxel::EVT_E][i] = m_numNonBoundaryCoefficients[Voxel::EVT_E]++;
			m_isBoundaryCoefficient[Voxel::EVT_E][i] = false;
		}
	}


	/*
	for(int i=0;i<10;++i)
	{
		int numValidCoeffs = 0;
		for( int j=0;j<numCoeffs;++j )
			if( m_localCoefficientIndices[i][j] >= 0 )
				++numValidCoeffs;
		std::cout << "type=" << typeStr(i) << " numValidCoeffs=" <<  numValidCoeffs << " numCoeffs=" << m_numNonBoundaryCoefficients[i] << std::endl;
		//m_numNonBoundaryCoefficients[i] = numValidCoeffs;
	}
	*/

	// create all voxels ----------------------
	int globalOffset = 0;
	for( int i=-m_numBoundaryLayers;i<m_resolution[0]+m_numBoundaryLayers;++i )
		for( int j=-m_numBoundaryLayers;j<m_resolution[1]+m_numBoundaryLayers;++j )
		{
			V2i coord(i, j);
			Voxel v;
			v.coord = coord;
			v.globalOffset = globalOffset;


			// determine voxeltype ==================
			if( sys->getDomain().contains_voxel(coord) )
				v.type = Voxel::EVT_Interior;
			else
			{
				V2i layer = getBoundaryLayer(v);
				if( (std::abs(layer[0]) > 1)||(std::abs(layer[1]) > 1) )
					v.type = Voxel::EVT_Boundary;
				else
				{
					// voxel is part of the first layer around the interior part
					// therefore it is of a mixed type and may contain boundary and
					// interior coefficients
					auto it = mixedTypes.find(std::make_pair(layer[0], layer[1]));
					if(it==mixedTypes.end())
					{
						std::cout << layer[0] << "  " << layer[1] << std::endl;
						throw std::runtime_error("unexpected");
					}
					v.type = it->second;
				}
			}


			m_voxels.push_back(v);
			globalOffset += getNumCoeffs(m_voxels.back());
		}
	// the value of globalOffset equals the number of cols and rows of our global system A
	m_numUnknowns = globalOffset;
}

PNSystem::Voxel& PNSystem::VoxelManager::getVoxel( const V2i& voxel )
{
	if( (voxel[0] < -m_numBoundaryLayers)||
		(voxel[1] < -m_numBoundaryLayers)||
		(voxel[0] >= m_resolution[0]+m_numBoundaryLayers)||
		(voxel[1] >= m_resolution[1]+m_numBoundaryLayers) )
	{
		std::cout << "test: " << voxel[0] << " " << voxel[1] << std::endl;
		throw std::runtime_error("PNSystem::getVoxel2 invalid voxel access");
	}
	int index = (voxel[0]+m_numBoundaryLayers)*(m_resolution[1]+2*m_numBoundaryLayers) + voxel[1] + m_numBoundaryLayers;
	return m_voxels[index];
}

int PNSystem::VoxelManager::getNumCoeffs(const Voxel& voxel)
{
	return m_numNonBoundaryCoefficients[voxel.type];
}


int PNSystem::VoxelManager::getGlobalIndex(const Voxel &voxel, int coeff)
{
	int localOffset = m_localCoefficientIndices[voxel.type][coeff];
	if( localOffset >= 0 )
		return voxel.globalOffset + localOffset;


	if( sys->getBoundaryConditions() == 0 )
		// dirichlet BC ---------------------
		// we simply lookup the valid coefficients, the invalid ones are ignored when setting
		// coefficients in A (which equates them to zero)
		return -1;
	else
	if( sys->getBoundaryConditions() == 1 )
	{
		// neumann bc -----------------
		// a boundary coefficient was accessed
		// we return the coefficient of the closest voxel within the domain
		V2i layer = getBoundaryLayer(voxel);

		if( (layer[0] == 0)&&(layer[1] == 0) )
			// means voxel is inside the domain where all coefficients should be defined
			throw std::runtime_error("unexpected layer coordinate");
		else
		{
			V2i step = V2i(-sgn(layer[0]), -sgn(layer[1]));
			return getGlobalIndex(sys->m_voxelManager.getVoxel(voxel.coord+step), coeff);
		}
	}

	throw std::runtime_error("unexpected");
}



bool PNSystem::VoxelManager::isBoundaryCoefficient(const Voxel& voxel, int coeff)
{
	return m_isBoundaryCoefficient[voxel.type][coeff];
}

V2i PNSystem::VoxelManager::getBoundaryLayer(const PNSystem::Voxel &v)
{
	V2i dist; // distance to the boundary
	for( int i=0;i<2;++i )
	{
		if( v.coord[i] < 0 )
			dist[i] = v.coord[i];
		else
		if( v.coord[i] >= sys->getDomain().resolution()[i] )
			dist[i] = v.coord[i] - sys->getDomain().resolution()[i] + 1;
		else
			dist[i] = 0;
	}
	return dist;
}
