#include <PNSystem.h>



PNSystem::PNSystem(const Domain& domain,
					int order)
	:
	m_domain(domain),
	m_order(order)
{
	m_numCoeffs = 0;
	for( int l=0;l<=order;++l )
		for( int m=-l;m<=l;++m )
			// in 2d, we only need to solve for moments where l+m is even
			if( (l+m) % 2 == 0 )
				++m_numCoeffs;

	m_matrix = Matrix(domain.numVoxels()*m_numCoeffs, domain.numVoxels()*m_numCoeffs);
}

void PNSystem::setField( const std::string& id, Field::Ptr field )
{
	if( id == "sigma_t" )
		m_fields.sigma_t = field;
	else
	if( id == "sigma_a" )
		m_fields.sigma_a = field;
	else
	if( id == "sigma_s" )
		m_fields.sigma_s = field;
	else
		throw std::runtime_error("PNSystem::setField unable to set field " + id);
}

PNSystem::Voxel PNSystem::getVoxel( const V2i& voxel )
{
	return PNSystem::Voxel(this, voxel);
}

const Domain& PNSystem::getDomain()const
{
	return m_domain;
}

PNSystem::Matrix& PNSystem::get_A()
{
	return m_matrix;
}

int PNSystem::getGlobalIndex( V2i voxel, int coeff )const
{
	//int voxel_index = voxel[0]*m_domain.resolution()[1] + voxel[1];
	//TODO: is this correct? shouldnt it be voxel[0]*res[1] + voxel[1] ?
	int voxel_index = voxel[1]*m_domain.resolution()[0] + voxel[0];
	return voxel_index*m_numCoeffs + coeff;
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
	MatrixAccessHelper mah(this);
	mah.m_global_i = getGlobalIndex(voxel_i, coefficient_i);
	mah.m_global_j = getGlobalIndex(voxel_j, coefficient_j);
	return mah;

	//int index = global_index_i*m_domain.resolution()[1]*m_numCoeffs + global_index_j;
	//return m_A_ptr[index];
}

/*
PNSystem::MatrixAccessHelper PNSystem::b( V2i voxel_i,
								   int coefficient_i )
{
	int global_index_i = global_index(voxel_i, coefficient_i);
	return m_b_ptr[global_index_i];
}
*/

void PNSystem::build()
{
	int res_x = m_domain.resolution()[0];
	int res_y = m_domain.resolution()[1];
	//for( int i=0;i<res_x;++i )
	//	for( int j=0;j<res_y;++j )
	V2i debug_voxel(35,35);
	for( int i=debug_voxel.x();i<debug_voxel.x()+1;++i )
		for( int j=debug_voxel.y();j<debug_voxel.y()+1;++j )
		{
			Voxel voxel = getVoxel(V2i(i,j));
			set_system_row( voxel, m_fields );
		}
	m_matrix.setFromTriplets(m_triplets.begin(), m_triplets.end());
}

// ==============================================
PNSystem::Voxel::Voxel(PNSystem* pns,
				const V2i& voxel_i ):
	m_pns(pns),
	m_voxel_i(voxel_i)
{
}


const V2i& PNSystem::Voxel::getVoxel()const
{
	return m_voxel_i;
}



PNSystem::MatrixAccessHelper PNSystem::Voxel::A( int coefficient_i, V2i voxel_j, int coefficient_j )
{
	return m_pns->A( m_voxel_i, coefficient_i, voxel_j, coefficient_j );
}

/*
std::complex<double>& PNSystem::Voxel::b( int coefficient_i )
{
	return m_pns->b( m_voxel_i, coefficient_i );
}
*/

V2d PNSystem::Voxel::voxelToWorld(const V2d& pVS )const
{
	return m_pns->getDomain().voxelToWorld(pVS);
}

