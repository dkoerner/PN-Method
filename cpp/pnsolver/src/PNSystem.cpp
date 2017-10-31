#include <PNSystem.h>
#include <field/Function.h>
#include <field/VoxelGridField.h> // used only for getting maximum value of sigma_t

#include<Eigen/IterativeLinearSolvers>
#include<Eigen/SparseLU>
//#include<Eigen/SparseCholesky>

#include <math/rng.h>
#include <util/threadpool.h>



std::map<std::string, PNSystem::Stencil> PNSystem::g_stencils;




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
	m_voxelManager.init(m_domain.getResolution(), m_stencil, neumannBC);


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


PNSystem::VoxelManager& PNSystem::getVoxelManager()
{
	return m_voxelManager;
}


V3i PNSystem::getGridSpaceOffsetFromGrid(int grid_index)
{
	static const V3i g_grid_offsets[8] ={ V3i(0, 0, 1),
										  V3i(1, 0, 1),
										  V3i(1, 1, 1),
										  V3i(0, 1, 1),
										  V3i(0, 0, 0),
										  V3i(1, 0, 0),
										  V3i(1, 1, 0),
										  V3i(0, 1, 0)};
	return g_grid_offsets[grid_index];
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

void PNSystem::setFields(PNSystem::Fields fields)
{
	m_fields = fields;
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
	return m_A;
}

PNSystem::RealVector& PNSystem::get_b_real()
{
	return m_b;
}

PNSystem::RealMatrix PNSystem::get_A_real_test()
{
	return m_A.transpose()*m_A;
}

Eigen::VectorXd PNSystem::get_b_real_test()
{
	return m_A.transpose()*m_b;
}


void PNSystem::setDebugVoxel(const V3i &dv)
{
	debugVoxel = dv;

}


int PNSystem::getNumCoefficients()const
{
	return m_stencil.numCoeffs;
}

V3i PNSystem::getResolution()const
{
	return m_domain.getResolution();
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
	MatrixBuilderd builder_A;
	MatrixBuilderd builder_b;


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
		m_stencil.apply(Stencil::Context(*this, v, &builder_A, &builder_b));
	}

	m_A = builder_A.build(m_voxelManager.getNumUnknowns(), m_voxelManager.getNumUnknowns());
	m_b = builder_b.build(m_voxelManager.getNumUnknowns(), 1);
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
	RealMatrix m = indexMatrix.build(res[0]+numBoundaryLayers[0]*2, res[1]+numBoundaryLayers[1]*2);

	//std::cout << indexMatrix.matrix.cols() << " " << indexMatrix.matrix.rows() << std::endl;
	return m;
}

int getIndex( const V3i& coord, int coeff_index, const V3i& res, int numCoeffs )
{

}


void setCoeff( PNSystem::MatrixBuilderd& builder, int i, int j, double value )
{
	if( j<0 )
		return;
	builder.coeff(i,j) += value;
}

Eigen::VectorXd PNSystem::removeStaggering(const Eigen::VectorXd &x)
{
	/*
	Eigen::VectorXd x2( m_domain.numVoxels()*m_stencil.numCoeffs,1);
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
			int index = m_voxelManager.getGlobalIndexBC(v, i);

			if(index==-1)
				throw std::runtime_error("PNSystem::stripBoundary didnt expect invalid coefficient index");
			x2.coeffRef(new_global_index, 0) = x.coeffRef(index, 0);
		} // for each coefficient
	} // for each voxel
	return x2;
	*/

	///*
	PNSystem::MatrixBuilderd builder;
	VoxelManager& vm = getVoxelManager();

	Eigen::VectorXd x2( m_domain.numVoxels()*m_stencil.numCoeffs,1);
	std::vector<Voxel>& voxels = m_voxelManager.getVoxels();
	V3i res = m_domain.getResolution();
	for( auto&v:voxels )
	{
		if( !m_domain.contains_voxel(v.coord) )
			// boundary or mixed boundary voxel
			continue;
		for( int coeff_index=0;coeff_index<m_stencil.numCoeffs;++coeff_index )
		{
			// this is the new index, which wont include boundary voxels
			int new_voxel_index = v.coord[0]*res[2]*res[1] + v.coord[1]*res[2] + v.coord[2];
			// since we store only internal voxels, we know that they all have all coefficients defined
			int new_global_i = new_voxel_index*m_stencil.numCoeffs + coeff_index;


			// in addition, we interpolate such that the coefficient is located at the voxel center
			// for that we need the coefficient location
			int grid_index = m_stencil.getGridIndexFromCoefficient(coeff_index);
			std::vector<int> indices;
			switch(grid_index)
			{
				case 0: // 0, 0, 1
				{
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(0,0,0), coeff_index));
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(0,1,0), coeff_index));
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(1,0,0), coeff_index));
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(1,1,0), coeff_index));
				}break;
				case 1: // 1, 0, 1
				{
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(0,0,0), coeff_index));
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(0,1,0), coeff_index));
				}break;
				case 2: // 1, 1, 1
				{
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(0,0,0), coeff_index));
				}break;
				case 3: // 0, 1, 1
				{
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(0,0,0), coeff_index));
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(1,0,0), coeff_index));
				}break;
				case 4: // 0, 0, 0
				{
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(0,0,0), coeff_index));
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(1,0,0), coeff_index));
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(0,1,0), coeff_index));
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(0,0,1), coeff_index));
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(0,1,1), coeff_index));
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(1,0,1), coeff_index));
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(1,1,0), coeff_index));
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(1,1,1), coeff_index));
				}break;
				case 5: // 1, 0, 0
				{
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(0,0,0), coeff_index));
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(0,1,0), coeff_index));
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(0,0,1), coeff_index));
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(0,1,1), coeff_index));
				}break;
				case 6: // 1, 1, 0
				{
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(0,0,0), coeff_index));
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(0,0,1), coeff_index));
				}break;
				case 7: // 0, 1, 0
				{
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(0,0,0), coeff_index));
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(1,0,0), coeff_index));
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(0,0,1), coeff_index));
					indices.push_back(vm.getGlobalIndex(v.coord+V3i(1,0,1), coeff_index));
				}break;
			};

			///*
			int numValidIndices = 0;
			for( auto& index:indices )
				if( index>=0 )
					++numValidIndices;
			double weight = 1.0/double(numValidIndices);
			for( auto& index:indices )
			{
				if( index>=0 )
					builder.coeff(new_global_i, index) += weight;
			}
			//*/


			// this is the current index, which includes boundary voxels
			int index = m_voxelManager.getGlobalIndexBC(v, coeff_index);
			if(index==-1)
				throw std::runtime_error("PNSystem::stripBoundary didnt expect invalid coefficient index");
			x2.coeffRef(new_global_i, 0) = x.coeffRef(index, 0);
		} // for each coefficient
	} // for each voxel


	MatrixBuilderd::Matrix conv = builder.build( m_domain.numVoxels()*m_stencil.numCoeffs, getVoxelManager().getNumUnknowns() );
	return conv*x;
	//*/
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
	Voxel& voxel_i = voxel;

	if( sys.getVoxelManager().voxelIsValid(voxel_j) )
	{
		// check if row for coefficient i is actually supposed to be set by the stencil
		// This is not the case with corner boundary voxels, where the equation for some coefficients are
		// descibing the average with neighbouring coefficients. In this case we prevent the stencil from
		// writing into those rows
		if(!sys.getVoxelManager().isBoundaryCoefficient(voxel_i, coeff_i))
		{
			// getGlobalIndex applies BC, which means it either returns -1 or bends indices
			int globalIndex_i = sys.getVoxelManager().getGlobalIndexBC(voxel_i, coeff_i);
			int globalIndex_j = sys.getVoxelManager().getGlobalIndexBC(sys.getVoxelManager().getVoxel(voxel_j), coeff_j);

			if( (globalIndex_i >= 0)&&(globalIndex_j >= 0) )
				return m_builder_A->coeff(globalIndex_i, globalIndex_j);
		}
	}

	return PNSystem::MatrixBuilderd::MatrixAccessHelper();
}

PNSystem::MatrixBuilderd::MatrixAccessHelper PNSystem::Stencil::Context::coeff_b(int coeff_i)
{
	Voxel& voxel_i = voxel;

	// check if row for coefficient i is actually supposed to be set by the stencil
	// This is not the case with corner boundary voxels, where the equation for some coefficients are
	// descibing the average with neighbouring coefficients. In this case we prevent the stencil from
	// writing into those rows
	if(!sys.getVoxelManager().isBoundaryCoefficient(voxel_i, coeff_i))
	{
		int globalIndex_i = sys.getVoxelManager().getGlobalIndexBC(voxel_i, coeff_i);
		if( globalIndex_i >= 0 )
			return m_builder_b->coeff(globalIndex_i, 0);
	}
	return PNSystem::MatrixBuilderd::MatrixAccessHelper();
}


// VoxelManager -------------------------------------------------

PNSystem::VoxelManager::VoxelManager()
{
}

void PNSystem::VoxelManager::init(const V3i& resolution, const Stencil& stencil, int boundaryConditions)
{
	int numCoeffs = stencil.numCoeffs;
	int boundaryWidth = stencil.width;
	std::vector<V3i> coefficientGridSpaceOffsets;
	for( int i=0;i<numCoeffs;++i )
		coefficientGridSpaceOffsets.push_back(stencil.getOffset(i));
	init( resolution,
		  numCoeffs,
		  boundaryWidth,
		  coefficientGridSpaceOffsets,
		  boundaryConditions);
}

void PNSystem::VoxelManager::init(const V3i& resolution,
								  int numCoeffs,
								  int boundaryWidth,
								  const std::vector<V3i>& coefficientGridSpaceOffsets,
								  int boundaryConditions)
{
	m_numCoeffs = numCoeffs;
	m_boundaryWidth = boundaryWidth;
	m_coefficientGridSpaceOffsets = std::vector<V3i>(coefficientGridSpaceOffsets.begin(), coefficientGridSpaceOffsets.end());
	m_resolution = resolution;
	m_boundaryConditions = boundaryConditions;
	//m_stencil = &stencil;

	if( m_resolution[2] == 1 )
		// in 2d we have zero boundary layers
		m_numBoundaryLayers = V3i(boundaryWidth, boundaryWidth, 0);
	else
		m_numBoundaryLayers = V3i(boundaryWidth, boundaryWidth, boundaryWidth);

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
				m_voxelTypes.push_back(VoxelType(m_voxelTypes.size(), numCoeffs, boundaryLayer));
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
		//V3i offset = stencil.getOffset(i);
		V3i offset = coefficientGridSpaceOffsets[i];

		// interior cells have all coefficients defined, therefore their indices are all set
		getVoxelType(0,0,0).registerActiveCoefficient(i);


		if(boundaryWidth > 0)
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
				throw std::runtime_error("PNSystem::VoxelManager::init unexpected");
		} // if boundaryWidth > 0 (if we need to register active coefficients in boundary voxels)

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


int PNSystem::VoxelManager::getGlobalIndexBC(const Voxel &voxel, int coeff)
{
	VoxelType& vt = m_voxelTypes[voxel.type];

	//int localOffset = m_localCoefficientIndices[voxel.type][coeff];
	int localOffset = vt.getLocalCoefficientIndex(coeff);
	if( localOffset >= 0 )
		return voxel.globalOffset + localOffset;


	if( m_boundaryConditions == 0 )
		// dirichlet BC ---------------------
		// we simply lookup the valid coefficients, the invalid ones are ignored when setting
		// coefficients in A (which equates them to zero)
		return -1;
	else
	if( m_boundaryConditions == 1 )
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
			return getGlobalIndexBC(getVoxel(voxel.coord+step), coeff);
		}
	}

	throw std::runtime_error("PNSystem::VoxelManager::getGlobalIndex unexpected");
}

int PNSystem::VoxelManager::getGlobalIndex(const Voxel &voxel, int coeff)
{
	VoxelType& vt = m_voxelTypes[voxel.type];
	int localOffset = vt.getLocalCoefficientIndex(coeff);
	if( localOffset >= 0 )
		return voxel.globalOffset + localOffset;
	return -1;
}

// returns global index as is
int PNSystem::VoxelManager::getGlobalIndex(const V3i &coord, int coeff )
{
	return getGlobalIndex(getVoxel(coord), coeff);
}

/*
int PNSystem::VoxelManager::getGlobalIndex( const V3i& coord, int coeff )
{
	Voxel& voxel = getVoxel(coord);
	return getGlobalIndex();
}
*/


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
		if( v.coord[i] >= m_resolution[i] )
			dist[i] = v.coord[i] - m_resolution[i] + 1;
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
