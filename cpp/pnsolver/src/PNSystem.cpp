#include <PNSystem.h>
#include <field/Function.h>

#include<Eigen/IterativeLinearSolvers>


PNSystem::PNSystem(const Domain& domain,
				   int order)
	:
	m_domain(domain),
	m_order(order),
	m_fields(order)
{
	m_numCoeffs = 0;
	for( int l=0;l<=order;++l )
		for( int m=-l;m<=l;++m )
			// in 2d, we only need to solve for moments where l+m is even
			if( (l+m) % 2 == 0 )
			{
				m_lm_to_index[std::make_pair(l,m)] = m_numCoeffs;
				++m_numCoeffs;
			}

	m_A_complex = ComplexMatrix(domain.numVoxels()*m_numCoeffs, domain.numVoxels()*m_numCoeffs);
	m_b_complex = ComplexVector(domain.numVoxels()*m_numCoeffs, 1);
	m_A_real = RealMatrix(domain.numVoxels()*m_numCoeffs, domain.numVoxels()*m_numCoeffs);
	m_b_real = RealVector(domain.numVoxels()*m_numCoeffs, 1);
	build_S();

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

	m_fields.sigma_a = std::make_shared<Function>(sigma_a);
	m_fields.sigma_s = std::make_shared<Function>(sigma_s);
	m_fields.sigma_t = std::make_shared<Function>(sigma_t);
	m_fields.f_p->setField(0,0,std::make_shared<Function>(phase_shcoeff));
	m_fields.q->setField(0,0,std::make_shared<Function>(source_shcoeffs));

}

void PNSystem::setField( const std::string& id, Field::Ptr field )
{
	std::cout << "WARNING: PNSystem::setField currently ignored. Using explicit RTE functions.\n";
	/*
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
	*/
}

PNSystem::VoxelSystem PNSystem::getVoxelSystem( const V2i& voxel )
{
	return PNSystem::VoxelSystem(this, voxel);
}

const Domain& PNSystem::getDomain()const
{
	return m_domain;
}

PNSystem::ComplexMatrix& PNSystem::get_A_complex()
{
	return m_A_complex;
}

PNSystem::RealMatrix& PNSystem::get_A_real()
{
	return m_A_real;
}

PNSystem::ComplexVector& PNSystem::get_b_complex()
{
	return m_b_complex;
}

PNSystem::RealVector& PNSystem::get_b_real()
{
	return m_b_real;
}

void PNSystem::setDebugVoxel(const V2i &dv)
{
	debugVoxel = dv;

}


void PNSystem::build_S()
{
	std::vector<ComplexTriplet> triplets_S;

	int count = 0;
	for( int l=0;l<=m_order;++l )
		for( int m=l;m>=0;--m )
		{
			// in 2D, we skip coefficients for which l+m is odd
			if( (l+m)%2 != 0 )
				continue;

			// build S matrix, which converts solution from complex to real values

			// computes the real part coefficients for a row (defined by l,m) in the S matrix
			// (see bottom of p.5 in the starmap paper)

			if(m == 0)
			{
				triplets_S.push_back( ComplexTriplet(count, getIndex(l,m), 1.0) );
			}else
			{
				triplets_S.push_back( ComplexTriplet(count, getIndex(l,m), std::pow(-1.0, m)/std::sqrt(2)) );
				if( m_lm_to_index.find(std::make_pair(l,-m)) != m_lm_to_index.end() )
				{
					triplets_S.push_back( ComplexTriplet(count, getIndex(l,-m), std::pow(-1.0, 2.0*m)/std::sqrt(2) ) );
				}
			}
			count+=1;

			// computes the imaginary part coefficients for a row (defined by l,m) in the S matrix
			// (see bottom of p.5 in the starmap paper)
			if( m > 0)
			{
				triplets_S.push_back( ComplexTriplet(count, getIndex(l,m), Complex(0.0, std::pow(-1.0, m)/std::sqrt(2)) ) );
				if( m_lm_to_index.find(std::make_pair(l,-m)) != m_lm_to_index.end() )
				{
					triplets_S.push_back( ComplexTriplet(count, getIndex(l,-m), Complex(0.0, -std::pow(-1.0, 2*m)/std::sqrt(2)) ) );
				}
				count+=1;
			}
		}


	// here we compute the inverse of S
	// and extract the triplets for sparse construction
	std::vector<ComplexTriplet> triplets_S_inv;
	{
		ComplexMatrix S(m_numCoeffs, m_numCoeffs);
		S.setFromTriplets(triplets_S.begin(), triplets_S.end());

		ComplexMatrix S_inv = Eigen::Matrix<Complex,Eigen::Dynamic,Eigen::Dynamic>(S).inverse().sparseView();
		// extract triplets
		for(int k=0; k < S_inv.outerSize(); ++k)
			for (ComplexMatrix::InnerIterator it(S_inv,k); it; ++it)
				triplets_S_inv.push_back( ComplexTriplet(it.row(), it.col(), it.value()) );
	}

	int numVoxels = m_domain.numVoxels();
	m_S = ComplexMatrix(m_numCoeffs*numVoxels, m_numCoeffs*numVoxels);
	m_S_inv = ComplexMatrix(m_numCoeffs*numVoxels, m_numCoeffs*numVoxels);

	std::vector<ComplexTriplet> triplets_S_big;
	std::vector<ComplexTriplet> triplets_S_inv_big;
	for( int i=0;i<numVoxels;++i )
	{
		int global_i = i*m_numCoeffs;
		int global_j = i*m_numCoeffs;
		for( auto t : triplets_S )
			triplets_S_big.push_back( ComplexTriplet(global_i+t.row(), global_j+t.col(), t.value()) );
		for( auto t : triplets_S_inv )
			triplets_S_inv_big.push_back( ComplexTriplet(global_i+t.row(), global_j+t.col(), t.value()) );
	}
	m_S.setFromTriplets(triplets_S_big.begin(), triplets_S_big.end());
	m_S_inv.setFromTriplets(triplets_S_inv_big.begin(), triplets_S_inv_big.end());
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

int PNSystem::getNumCoefficients()const
{
	return m_numCoeffs;
}

int PNSystem::getNumVoxels()const
{
	return m_domain.numVoxels();
}

PNSystem::MatrixAccessHelper PNSystem::A( V2i voxel_i,
										  int coefficient_i,
										  V2i voxel_j,
										  int coefficient_j )
{
	if( !m_domain.contains_voxel(voxel_i) ||
		!m_domain.contains_voxel(voxel_j) )
	{
		// the voxel to be accessed by the stencil is out of bounds
		// we return a MatrixAccessHelper with a zero pointer to a
		// triplet list. All modifications/operations will have no effect.
		return MatrixAccessHelper(0);
	}

	MatrixAccessHelper mah(&m_triplets_A);
	mah.m_global_i = getGlobalIndex(voxel_i, coefficient_i);
	mah.m_global_j = getGlobalIndex(voxel_j, coefficient_j);
	return mah;
}


PNSystem::MatrixAccessHelper PNSystem::b( V2i voxel_i,
										  int coefficient_i )
{
	if( !m_domain.contains_voxel(voxel_i) )
	{
		// the voxel to be accessed by the stencil is out of bounds
		// we return a MatrixAccessHelper with a zero pointer to a
		// triplet list. All modifications/operations will have no effect.
		return MatrixAccessHelper(0);
	}

	MatrixAccessHelper mah(&m_triplets_b);
	mah.m_global_i = getGlobalIndex(voxel_i, coefficient_i);
	mah.m_global_j = 0;
	return mah;
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
	for( int i=min_x;i<max_x;++i )
		for( int j=min_y;j<max_y;++j )
			set_system_row( getVoxelSystem(V2i(i,j)), m_fields );

	// construct sparse matrix from all stencil operations ---
	m_A_complex.setFromTriplets(m_triplets_A.begin(), m_triplets_A.end());
	m_b_complex.setFromTriplets(m_triplets_b.begin(), m_triplets_b.end());

	// convert from complex to real valued system
	m_A_real = (m_S*m_A_complex*m_S_inv).real();
	m_b_real = (m_S*m_b_complex).real();
}

PNSystem::RealVector PNSystem::solve()
{
	std::cout << "PNSystem::solve solving for x...\n";

	Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver;
	solver.compute(get_A_real());
	if(solver.info()!=Eigen::Success)
	{
		throw std::runtime_error("PNSystem::solve decomposition failed");
	}
	RealVector x = solver.solve(get_b_real());
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

PNSystem::MatrixAccessHelper PNSystem::VoxelSystem::A( int coefficient_i, V2i voxel_j, int coefficient_j )
{
	return m_pns->A( m_voxel_i, coefficient_i, voxel_j, coefficient_j );
}

PNSystem::MatrixAccessHelper PNSystem::VoxelSystem::b( int coefficient_i )
{
	return m_pns->b( m_voxel_i, coefficient_i );
}

V2d PNSystem::VoxelSystem::voxelToWorld(const V2d& pVS )const
{
	return m_pns->getDomain().voxelToWorld(pVS);
}

