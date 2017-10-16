#include <PNSystem.h>
#include <field/Function.h>
#include <field/VoxelGridField.h> // used only for getting maximum value of sigma_t

#include<Eigen/IterativeLinearSolvers>
#include<Eigen/SparseLU>
//#include<Eigen/SparseCholesky>

#include <math/rng.h>
#include <util/threadpool.h>



std::map<std::string, PNSystem::Stencil> PNSystem::g_stencils;


V3i PNSystem::Stencil::Context::g_grid_offsets[8] ={V3i(0, 0, 1),
													V3i(1, 0, 1),
													V3i(1, 1, 1),
													V3i(0, 1, 1),
													V3i(0, 0, 0),
													V3i(1, 0, 0),
													V3i(1, 1, 0),
													V3i(0, 1, 0)};


V3d PNSystem::Stencil::Context::g_grid_offsets2[8] ={V3d(0.0, 0.0, 0.5),
													V3d(0.5, 0.0, 0.5),
													V3d(0.5, 0.5, 0.5),
													V3d(0.0, 0.5, 0.5),
													V3d(0.0, 0.0, 0.0),
													V3d(0.5, 0.0, 0.0),
													V3d(0.5, 0.5, 0.0),
													V3d(0.0, 0.5, 0.0)};


PNSystem::PNSystem( Stencil stencil, const Domain& domain, bool neumannBC)
	:
	m_domain(domain),
	m_fields(stencil.order),
	m_stencil(stencil),
	m_neumannBC(neumannBC),
	m_voxelManager()
{
	// find out number of coefficients and mapping between sh bands l,m and linear indices ---
	//m_numCoeffs = 0;


	// voxelmanager needs to be initialized after the number of coefficients has been determined
	m_voxelManager.init(this);


	debugVoxel = V3i(-1, -1);

	// TEMP: testing with non-rasterized rte parameters
	Function::FunctionType sigma_a = [](const P3d& pWS)
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
		//return (1.0-g)*1 + g*10;
		//return 6.0;
	};
	Function::FunctionType sigma_s = [](const P3d& pWS)
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
		//return 0.0;
	};
	Function::FunctionType sigma_t = [=](const P3d& pWS)
	{
		return sigma_a(pWS)+sigma_s(pWS);
	};
	Function::FunctionType phase_shcoeff = [](const P3d& pWS)
	{
		return 1.0;
	};
	Function::FunctionType source_shcoeffs = [](const P3d& pWS)
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


PNSystem::MatrixBuilderd::MatrixAccessHelper PNSystem::coeff_A(V3i voxel_i, int coefficient_i, V3i voxel_j, int coefficient_j, bool debug)
{
	if( m_voxelManager.voxelIsValid(voxel_i) && m_voxelManager.voxelIsValid(voxel_j) )
	{
		int globalIndex_i = m_voxelManager.getGlobalIndex(m_voxelManager.getVoxel(voxel_i), coefficient_i);
		int globalIndex_j = m_voxelManager.getGlobalIndex(m_voxelManager.getVoxel(voxel_j), coefficient_j);

		if( (globalIndex_i >= 0)&&(globalIndex_j >= 0) )
			return m_builder_A.coeff(globalIndex_i, globalIndex_j, debug);
	}

	return PNSystem::MatrixBuilderd::MatrixAccessHelper();
}

PNSystem::MatrixBuilderd::MatrixAccessHelper PNSystem::coeff_b(V3i voxel_i, int coefficient_i, bool debug)
{
	int globalIndex_i = m_voxelManager.getGlobalIndex(m_voxelManager.getVoxel(voxel_i), coefficient_i);

	if( globalIndex_i >= 0 )
		return m_builder_b.coeff(globalIndex_i, 0, debug);
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
/*
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
			double sigma_t_max = std::dynamic_pointer_cast<VoxelGridField>(fields->sigma_t)->getMax();

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


			// SCATTERING ==============================
			// create new direction (assuming isotropic scattering)
			{
				V3d d_sphere = sampleSphere<double>(rng);
				d = V2d(d_sphere[0], d_sphere[1]);
				//pdf_over_throughput /= sampleSpherePDF();
			}
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
*/
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

void PNSystem::setDebugVoxel(const V3i &dv)
{
	debugVoxel = dv;

}


int PNSystem::getNumCoefficients()const
{
	return m_stencil.numCoeffs;
}

int PNSystem::getNumVoxels()const
{
	return m_domain.numVoxels();
}

int PNSystem::getOrder() const
{
	return m_stencil.order;
}

void PNSystem::build()
{
	std::cout << "PNSystem::build building system matrices A and b...\n";
	int min_x = 0;
	int max_x = m_domain.getResolution()[0];
	int min_y = 0;
	int max_y = m_domain.getResolution()[1];
	int min_z = 0;
	int max_z = m_domain.getResolution()[2];

	if( (debugVoxel.x() != -1)&&(debugVoxel.y() != -1)&&(debugVoxel.z() != -1) )
	{
		min_x = debugVoxel.x();
		max_x = min_x+1;
		min_y = debugVoxel.y();
		max_y = min_y+1;
		min_z = debugVoxel.z();
		max_z = min_z+1;
	}

	// set debug voxels
	//m_voxelManager.getVoxel(V3i(21, 38, 0)).debug = true;
	//m_voxelManager.getVoxel(V3i(71-21-1 + 1, 38, 0)).debug = true;


	// iterate all voxels and apply stencil which has been generated from python script ---
	//for( int i=min_x;i<max_x;++i )
	//	for( int j=min_y;j<max_y;++j )
	//		m_stencil.apply( getVoxelSystem(V2i(i,j)), m_fields );
	for( auto& v: m_voxelManager.getVoxels() )
	{
		if(v.debug)
			std::cout << "Voxel=" << v.coord[0] << " " << v.coord[1] << " " << v.coord[2] << "==============================\n";
		m_stencil.apply(Stencil::Context(*this, v));
	}

	m_builder_A.build(m_voxelManager.getNumUnknowns(), m_voxelManager.getNumUnknowns());
	m_builder_b.build(m_voxelManager.getNumUnknowns(), 1);
}



PNSystem::RealMatrix PNSystem::getVoxelInfo(const std::string &info)
{
	MatrixBuilderd indexMatrix;

	V3i res = m_domain.getResolution();
	V3i numBoundaryLayers = m_voxelManager.getNumBoundaryLayers();
	for( int i=-numBoundaryLayers[0];i<res[0]+numBoundaryLayers[0];++i )
		for( int j=-numBoundaryLayers[1];j<res[1]+numBoundaryLayers[1];++j )
			for( int k=-numBoundaryLayers[2];k<res[2]+numBoundaryLayers[2];++k )
			{
				//std::cout << "i=" << i << " j=" << j << " k=" << k << std::endl; std::flush(std::cout);
				V3i coord(i, j, k);
				Voxel& v = m_voxelManager.getVoxel(coord);

				//std::cout << "voxel=" << v.coord[0] << " " << v.coord[1] << " " << v.coord[2] << std::endl;

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
				if( info == "coord.z" )
					value = double(int(v.coord[2]));

				indexMatrix.coeff(i+numBoundaryLayers[0],j+numBoundaryLayers[1]) += value;
			}

	//std::cout << "!!!!!?????\n"; std::flush(std::cout);
	indexMatrix.build(res[0]+numBoundaryLayers[0]*2, res[1]+numBoundaryLayers[1]*2);

	//std::cout << indexMatrix.matrix.cols() << " " << indexMatrix.matrix.rows() << std::endl;
	return indexMatrix.matrix;
}

PNSystem::RealVector PNSystem::solve()
{
	std::cout << "PNSystem::solve solving for x...\n";

	///*
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
	//*/

	/*
	Eigen::ConjugateGradient<RealMatrix> solver;
	//Eigen::BiCGSTAB<RealMatrix> solver;
	//Eigen::SparseLU<RealMatrix> solver;
	//solver.setMaxIterations(1000);
	std::cout << "PNSystem::solve compute...\n";std::flush(std::cout);
	solver.compute(get_A_real().transpose()*get_A_real());
	if(solver.info()!=Eigen::Success)
	{
		throw std::runtime_error("PNSystem::solve decomposition failed");
	}
	std::cout << "PNSystem::solve solve...\n";std::flush(std::cout);
	RealVector x = solver.solve(get_A_real().transpose()*get_b_real());
	//std::cout << "iterations=" << solver.iterations() << std::endl;
	if(solver.info()!=Eigen::Success)
	{
		throw std::runtime_error("PNSystem::solve solve failed");
	}
	*/



	// ..................
	RealVector x2( m_domain.numVoxels()*m_stencil.numCoeffs,1);
	std::vector<RealTriplet> triplets;
	std::vector<Voxel>& voxels = m_voxelManager.getVoxels();
	for( auto&v:voxels )
	{
		if( !m_domain.contains_voxel(v.coord) )
			// boundary or mixed boundary voxel
			continue;
		for( int i=0;i<m_stencil.numCoeffs;++i )
		{
			// this is the new index, which wont include boundary voxels
			V3i coord = v.coord;
			V3i res = m_domain.getResolution();
			int new_voxel_index = coord[0]*res[2]*res[1] + coord[1]*res[2] + coord[2];
			// since we store only internal voxels, we know that they all have all coefficients defined
			int new_global_index = new_voxel_index*m_stencil.numCoeffs + i;


			// this is the current index, which includes boundary voxels
			int index = m_voxelManager.getGlobalIndex(v, i);

			if(index==-1)
				throw std::runtime_error("PNSystem::solve didnt expect invalid coefficient index");
			triplets.push_back(RealTriplet(new_global_index, 0, x.coeffRef(index, 0)));

			if( (coord[0] == 21)&&(coord[1] == 38) )
			{
				//triplets.push_back(RealTriplet(new_global_index, 0, 1.0));
				//std::cout <<
			}

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

PNSystem::MatrixBuilderd::MatrixAccessHelper PNSystem::Stencil::Context::coeff_A(int coeff_i, const V3i &voxel_j, int coeff_j)
{
	// check if row for coefficient i is actually supposed to be set by the stencil
	// This is not the case with corner boundary voxels, where the equation for some coefficients are
	// descibing the average with neighbouring coefficients. In this case we prevent the stencil from
	// writing into those rows
	if(sys.m_voxelManager.isBoundaryCoefficient(voxel, coeff_i))
		return PNSystem::MatrixBuilderd::MatrixAccessHelper();

	return sys.coeff_A(voxel.coord, coeff_i, voxel_j, coeff_j, voxel.debug);
}

PNSystem::MatrixBuilderd::MatrixAccessHelper PNSystem::Stencil::Context::coeff_b(int coeff_i)
{
	return sys.coeff_b(voxel.coord, coeff_i, voxel.debug);
}


// VoxelManager -------------------------------------------------

PNSystem::VoxelManager::VoxelManager()
{
}

void PNSystem::VoxelManager::init(PNSystem *sys)
{
	this->sys = sys;
	int numCoeffs = sys->getNumCoefficients();
	m_resolution = sys->getDomain().getResolution();


	if( m_resolution[2] == 1 )
		// in 2d we have zero boundary layers
		m_numBoundaryLayers = V3i(sys->getStencil().width, sys->getStencil().width, 0);
	else
		m_numBoundaryLayers = V3i(sys->getStencil().width, sys->getStencil().width, sys->getStencil().width);

	for( int i=-m_numBoundaryLayers[0];i<=m_numBoundaryLayers[0];++i )
		for( int j=-m_numBoundaryLayers[1];j<=m_numBoundaryLayers[1];++j )
			for( int k=-m_numBoundaryLayers[2];k<=m_numBoundaryLayers[2];++k )
			{
				auto boundaryLayer = std::make_tuple(i, j, k);

				// check if type has already been registered..
				for( auto&vt : m_voxelTypes )
					if( vt.getBoundaryLayer() == boundaryLayer )
						continue;

				// type hasnt been registered, so we register it
				m_voxelTypes.push_back(VoxelType(m_voxelTypes.size(), sys->getNumCoefficients(), boundaryLayer));
			}

	int numTypes = m_voxelTypes.size();

	// establish reverse mapping, which we use during voxel creation
	for( int i=0;i<numTypes;++i )
	{
		auto& vt = m_voxelTypes[i];
		m_layerToVoxelTypeIndex[vt.getBoundaryLayer()] = i;
	}


	for( int i=0;i<numCoeffs;++i )
	{
		V3i offset = sys->m_stencil.getOffset(i);

		// interior cells have all coefficients defined, therefore their indices are all set
		getVoxelType(0,0,0).registerActiveCoefficient(i);


		if(sys->getStencil().width > 0)
		{
			// here we make those coefficients non-boundary(active) coefficients, which sit directly on the boundary
			// I think this is required for both, neumann and dirichlet BC
			if( (offset[0] == 0)&&(offset[1] == 0)&&(offset[2] == 1) )
			{
				getVoxelType(1,1,0).registerActiveCoefficient(i);
				getVoxelType(1,0,0).registerActiveCoefficient(i);
				getVoxelType(0,1,0).registerActiveCoefficient(i);
			}else
			if( (offset[0] == 0)&&(offset[1] == 1)&&(offset[2] == 1) )
			{
				getVoxelType(1,0,0).registerActiveCoefficient(i);
			}else
			if( (offset[0] == 1)&&(offset[1] == 0)&&(offset[2] == 1) )
			{
				getVoxelType(0,1,0).registerActiveCoefficient(i);
			}else
			if( (offset[0] == 1)&&(offset[1] == 1)&&(offset[2] == 1) )
			{
				// coefficients in the center will always be boundary/passive coefficients
			}else
			if( (offset[0] == 0)&&(offset[1] == 0)&&(offset[2] == 0) )
			{
				// coefficients at the intersection are active for north, east and forward boundary cells
				getVoxelType(1,0,0).registerActiveCoefficient(i);
				getVoxelType(0,1,0).registerActiveCoefficient(i);
				getVoxelType(0,0,1).registerActiveCoefficient(i);
				// and edgecorner cells
				getVoxelType(0,1,1).registerActiveCoefficient(i);
				getVoxelType(1,0,1).registerActiveCoefficient(i);
				getVoxelType(1,1,0).registerActiveCoefficient(i);
				// and the corner cell
				getVoxelType(1,1,1).registerActiveCoefficient(i);
			}else
			if( (offset[0] == 0)&&(offset[1] == 1)&&(offset[2] == 0) )
			{
				getVoxelType(1,0,1).registerActiveCoefficient(i);
				getVoxelType(1,0,0).registerActiveCoefficient(i);
				getVoxelType(0,0,1).registerActiveCoefficient(i);
			}else
			if( (offset[0] == 1)&&(offset[1] == 0)&&(offset[2] == 0) )
			{
				getVoxelType(0,1,1).registerActiveCoefficient(i);
				getVoxelType(0,1,0).registerActiveCoefficient(i);
				getVoxelType(0,0,1).registerActiveCoefficient(i);
			}else
			if( (offset[0] == 1)&&(offset[1] == 1)&&(offset[2] == 0) )
			{
				getVoxelType(0,0,1).registerActiveCoefficient(i);
			}else
				throw std::runtime_error("unexpected");
		} // if stencil.widh > 0 (if we need to register active coefficients in boundary voxels)

	}

	// create all voxels ----------------------
	int voxelIndex = 0;
	int globalOffset = 0;
	for( int i=-m_numBoundaryLayers[0];i<m_resolution[0]+m_numBoundaryLayers[0];++i )
		for( int j=-m_numBoundaryLayers[1];j<m_resolution[1]+m_numBoundaryLayers[1];++j )
			for( int k=-m_numBoundaryLayers[2];k<m_resolution[2]+m_numBoundaryLayers[2];++k, ++voxelIndex )
			{
				V3i coord(i, j, k);
				Voxel v;
				v.coord = coord;
				v.globalOffset = globalOffset;

				{
					V3i coord_shifted = coord + V3i(m_numBoundaryLayers[0], m_numBoundaryLayers[1], m_numBoundaryLayers[2]);
					V3i res_shifted = m_resolution + V3i(2*m_numBoundaryLayers[0], 2*m_numBoundaryLayers[1], 2*m_numBoundaryLayers[2]);
					int index = coord_shifted[0]*res_shifted[2]*res_shifted[1] + coord_shifted[1]*res_shifted[2] + coord_shifted[2];
					//std::cout << "check: index=" << voxelIndex << " " << index << std::endl;
				}


				// determine voxeltype ==================
				V3i boundaryLayer = getBoundaryLayer(v);
				v.type = getVoxelType(boundaryLayer[0], boundaryLayer[1], boundaryLayer[2]).getIndex();


				// add voxel and shift global offset index ========
				m_voxels.push_back(v);
				globalOffset += getNumCoeffs(m_voxels.back());
			}
	// the value of globalOffset equals the number of cols and rows of our global system A
	m_numUnknowns = globalOffset;
}

PNSystem::Voxel& PNSystem::VoxelManager::getVoxel( const V3i& coord )
{
	if( (coord[0] < -m_numBoundaryLayers[0])||
		(coord[1] < -m_numBoundaryLayers[1])||
		(coord[2] < -m_numBoundaryLayers[2])||
		(coord[0] >= m_resolution[0]+m_numBoundaryLayers[0])||
		(coord[1] >= m_resolution[1]+m_numBoundaryLayers[1])||
		(coord[2] >= m_resolution[2]+m_numBoundaryLayers[2]))
	{
		std::cout << "test: " << coord[0] << " " << coord[1] << " " << coord[2] << std::endl;
		throw std::runtime_error("PNSystem::getVoxel invalid voxel access");
	}
	V3i coord_shifted = coord + V3i(m_numBoundaryLayers[0], m_numBoundaryLayers[1], m_numBoundaryLayers[2]);
	V3i res_shifted = m_resolution + V3i(2*m_numBoundaryLayers[0], 2*m_numBoundaryLayers[1], 2*m_numBoundaryLayers[2]);
	int index = coord_shifted[0]*res_shifted[2]*res_shifted[1] + coord_shifted[1]*res_shifted[2] + coord_shifted[2];
	return m_voxels[index];
}

int PNSystem::VoxelManager::getNumCoeffs(const Voxel& voxel)
{
	//return m_numNonBoundaryCoefficients[voxel.type];
	return m_voxelTypes[voxel.type].getNumActiveCoefficients();
}


int PNSystem::VoxelManager::getGlobalIndex(const Voxel &voxel, int coeff)
{
	VoxelType& vt = m_voxelTypes[voxel.type];

	//int localOffset = m_localCoefficientIndices[voxel.type][coeff];
	int localOffset = vt.getLocalCoefficientIndex(coeff);
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
		V3i layer = getBoundaryLayer(voxel);

		if( (layer[0] == 0)&&(layer[1] == 0)&&(layer[2] == 0) )
			// means voxel is inside the domain where all coefficients should be defined
			throw std::runtime_error("unexpected layer coordinate");
		else
		{
			V3i step = V3i(-sgn(layer[0]), -sgn(layer[1]), -sgn(layer[2]));
			return getGlobalIndex(sys->m_voxelManager.getVoxel(voxel.coord+step), coeff);
		}
	}

	throw std::runtime_error("unexpected");
}



bool PNSystem::VoxelManager::isBoundaryCoefficient(const Voxel& voxel, int coeff)
{
	return !m_voxelTypes[voxel.type].isActiveCoefficient( coeff );
}

V3i PNSystem::VoxelManager::getNumBoundaryLayers() const
{
	return m_numBoundaryLayers;
}

V3i PNSystem::VoxelManager::getBoundaryLayer(const PNSystem::Voxel &v)
{
	V3i dist; // distance to the boundary
	for( int i=0;i<3;++i )
	{
		if( v.coord[i] < 0 )
			dist[i] = v.coord[i];
		else
		if( v.coord[i] >= sys->getDomain().getResolution()[i] )
			dist[i] = v.coord[i] - sys->getDomain().getResolution()[i] + 1;
		else
			dist[i] = 0;
	}
	return dist;
}


PNSystem::VoxelManager::VoxelType& PNSystem::VoxelManager::getVoxelType(int boundaryLayer_x, int boundaryLayer_y, int boundaryLayer_z)
{
	return getVoxelType(std::make_tuple(boundaryLayer_x, boundaryLayer_y, boundaryLayer_z));
}

PNSystem::VoxelManager::VoxelType &PNSystem::VoxelManager::getVoxelType(const std::tuple<int, int, int> &boundaryLayer)
{
	auto it = m_layerToVoxelTypeIndex.find(boundaryLayer);
	if( it != m_layerToVoxelTypeIndex.end() )
		return m_voxelTypes[it->second];
	std::cout << "boundaryLayer=" << std::get<0>(boundaryLayer) << " " << std::get<1>(boundaryLayer) << " " << std::get<2>(boundaryLayer) << std::endl;
	throw std::runtime_error("PNSystem::VoxelManager::getVoxelType: type does not exist for given boundary layer");
}


PNSystem::VoxelManager::VoxelType::VoxelType(int index, int maxNumCoefficients, const std::tuple<int, int, int>& boundaryLayer)
	:m_index(index),
	 m_boundaryLayer(boundaryLayer)
{
	m_localCoefficientIndices = std::vector<int>(maxNumCoefficients, -1);
	m_isBoundaryCoefficient = std::vector<bool>(maxNumCoefficients, true);
	m_numNonBoundaryCoefficients = 0;
}

int PNSystem::VoxelManager::VoxelType::getIndex() const
{
	return m_index;
}

void PNSystem::VoxelManager::VoxelType::registerActiveCoefficient(int coeff_index)
{
	m_localCoefficientIndices[coeff_index] = m_numNonBoundaryCoefficients++;
	m_isBoundaryCoefficient[coeff_index] = false;
}

std::tuple<int, int, int> &PNSystem::VoxelManager::VoxelType::getBoundaryLayer()
{
	return m_boundaryLayer;
}

int PNSystem::VoxelManager::VoxelType::getNumActiveCoefficients() const
{
	return m_numNonBoundaryCoefficients;
}

bool PNSystem::VoxelManager::VoxelType::isActiveCoefficient(int coeff_index) const
{
	return !m_isBoundaryCoefficient[coeff_index];
}

int PNSystem::VoxelManager::VoxelType::getLocalCoefficientIndex(int coeff_index) const
{
	return m_localCoefficientIndices[coeff_index];
}
