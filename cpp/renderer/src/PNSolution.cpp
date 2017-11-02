#include <PNSolution.h>
#include <math/sph.h>

#include <iostream>
#include <fstream>




int PNSolution::getIndex( const V3i& voxel_coord )const
{
	const int coeff_index = 0;
	// this is the new index, which wont include boundary voxels
	int new_voxel_index = voxel_coord[0]*m_resolution[2]*m_resolution[1] + voxel_coord[1]*m_resolution[2] + voxel_coord[2];
	// since we store only internal voxels, we know that they all have all coefficients defined
	return new_voxel_index*m_numCoeffs + coeff_index;
}


double PNSolution::eval( const P3d& pWS, const V3d& direction )const
{
	P3d pVS = worldToVoxel(pWS);

	// take sample location within voxel into account
	pVS -= V3d(0.5, 0.5, 0.5);

	double tx = pVS.x() - std::floor(pVS.x());
	double ty = pVS.y() - std::floor(pVS.y());
	double tz = pVS.z() - std::floor(pVS.z());

	// lower left corner
	V3i c1;
	c1[0] = (int)std::floor(pVS.x());
	c1[1] = (int)std::floor(pVS.y());
	c1[2] = (int)std::floor(pVS.z());

	// upper right corner
	V3i c2 = c1+V3i(1);

	// clamp the indexing coordinates
	c1[0] = std::max(0, std::min(c1.x(), m_resolution.x()-1));
	c2[0] = std::max(0, std::min(c2.x(), m_resolution.x()-1));
	c1[1] = std::max(0, std::min(c1.y(), m_resolution.y()-1));
	c2[1] = std::max(0, std::min(c2.y(), m_resolution.y()-1));
	c1[2] = std::max(0, std::min(c1.z(), m_resolution.z()-1));
	c2[2] = std::max(0, std::min(c2.z(), m_resolution.z()-1));

	int voxel_indices[8] = {getIndex(V3i(c1[0], c1[1], c1[2])),
							getIndex(V3i(c2[0], c1[1], c1[2])),
							getIndex(V3i(c1[0], c2[1], c1[2])),
							getIndex(V3i(c2[0], c2[1], c1[2])),
							getIndex(V3i(c1[0], c1[1], c2[2])),
							getIndex(V3i(c2[0], c1[1], c2[2])),
							getIndex(V3i(c1[0], c2[1], c2[2])),
							getIndex(V3i(c2[0], c2[1], c2[2])),
						   };


	// evaluate SH
	P2d theta_phi = sphericalCoordinates(direction);

	std::complex<double> result = 0.0;
	int coeff_index = 0;
	for( int l=0;l<=m_order;++l )
		for( int m=-l;m<=l;++m,++coeff_index )
		{
			std::complex<double> coeff = 0.0;

			//lerp coefficient
			coeff += m_data[voxel_indices[0] + coeff_index]*(1.0-tx)*(1.0-ty)*(1.0-tz);
			coeff += m_data[voxel_indices[1] + coeff_index]*tx*(1.0-ty)*(1.0-tz);
			coeff += m_data[voxel_indices[2] + coeff_index]*(1.0-tx)*ty*(1.0-tz);
			coeff += m_data[voxel_indices[3] + coeff_index]*tx*ty*(1.0-tz);
			coeff += m_data[voxel_indices[4] + coeff_index]*(1.0-tx)*(1.0-ty)*tz;
			coeff += m_data[voxel_indices[5] + coeff_index]*tx*(1.0-ty)*tz;
			coeff += m_data[voxel_indices[6] + coeff_index]*(1.0-tx)*ty*tz;
			coeff += m_data[voxel_indices[7] + coeff_index]*tx*ty*tz;

			// now do the SH accumulation
			result += coeff*sph::sph_basis(l, m, theta_phi[0], theta_phi[1]);
		}



	// NB: here we apply a normalization. This is required to have the intergal over
	// the sphere match the zero coefficient
	return result.real()/std::sqrt(4.0*M_PI);
}

std::complex<double> PNSolution::evalCoefficient(const P3d &pWS, int coeff_index)const
{
	P3d pVS = worldToVoxel(pWS);

	// take sample location within voxel into account
	pVS -= V3d(0.5, 0.5, 0.5);

	double tx = pVS.x() - std::floor(pVS.x());
	double ty = pVS.y() - std::floor(pVS.y());
	double tz = pVS.z() - std::floor(pVS.z());

	// lower left corner
	V3i c1;
	c1[0] = (int)std::floor(pVS.x());
	c1[1] = (int)std::floor(pVS.y());
	c1[2] = (int)std::floor(pVS.z());

	// upper right corner
	V3i c2 = c1+V3i(1);

	// clamp the indexing coordinates
	c1[0] = std::max(0, std::min(c1.x(), m_resolution.x()-1));
	c2[0] = std::max(0, std::min(c2.x(), m_resolution.x()-1));
	c1[1] = std::max(0, std::min(c1.y(), m_resolution.y()-1));
	c2[1] = std::max(0, std::min(c2.y(), m_resolution.y()-1));
	c1[2] = std::max(0, std::min(c1.z(), m_resolution.z()-1));
	c2[2] = std::max(0, std::min(c2.z(), m_resolution.z()-1));

	int voxel_indices[8] = {getIndex(V3i(c1[0], c1[1], c1[2])),
							getIndex(V3i(c2[0], c1[1], c1[2])),
							getIndex(V3i(c1[0], c2[1], c1[2])),
							getIndex(V3i(c2[0], c2[1], c1[2])),
							getIndex(V3i(c1[0], c1[1], c2[2])),
							getIndex(V3i(c2[0], c1[1], c2[2])),
							getIndex(V3i(c1[0], c2[1], c2[2])),
							getIndex(V3i(c2[0], c2[1], c2[2])),
						   };


	std::complex<double> coeff = 0.0;

	//lerp coefficient
	coeff += m_data[voxel_indices[0] + coeff_index]*(1.0-tx)*(1.0-ty)*(1.0-tz);
	coeff += m_data[voxel_indices[1] + coeff_index]*tx*(1.0-ty)*(1.0-tz);
	coeff += m_data[voxel_indices[2] + coeff_index]*(1.0-tx)*ty*(1.0-tz);
	coeff += m_data[voxel_indices[3] + coeff_index]*tx*ty*(1.0-tz);
	coeff += m_data[voxel_indices[4] + coeff_index]*(1.0-tx)*(1.0-ty)*tz;
	coeff += m_data[voxel_indices[5] + coeff_index]*tx*(1.0-ty)*tz;
	coeff += m_data[voxel_indices[6] + coeff_index]*(1.0-tx)*ty*tz;
	coeff += m_data[voxel_indices[7] + coeff_index]*tx*ty*tz;

	return coeff;
}


P3d PNSolution::voxelToLocal(const P3d& pVS)const
{
	return P3d(pVS[0]/m_resolution[0], pVS[1]/m_resolution[1], pVS[2]/m_resolution[2]);
}
P3d PNSolution::localToVoxel(const P3d& pLS)const
{
	return P3d(pLS[0]*m_resolution[0], pLS[1]*m_resolution[1], pLS[2]*m_resolution[2]);
}
P3d PNSolution::localToWorld( const P3d& pLS )const
{
	return P3d(pLS[0]*m_extend[0] + m_bound.min[0], pLS[1]*m_extend[1] + m_bound.min[1], pLS[2]*m_extend[2] + m_bound.min[2]);
}

P3d PNSolution::worldToLocal( const P3d& pWS )const
{
	return P3d((pWS[0]-m_bound.min[0])/m_extend[0], (pWS[1]-m_bound.min[1])/m_extend[1], (pWS[2]-m_bound.min[2])/m_extend[2] );
}

P3d PNSolution::voxelToWorld( const P3d& pVS)const
{
	return localToWorld(voxelToLocal(pVS));
}


P3d PNSolution::worldToVoxel( const P3d& pWS)const
{
	return localToVoxel(worldToLocal(pWS));
}

bool PNSolution::contains_voxel( const P3i& voxel )const
{
	if( (voxel[0] < 0)||(voxel[0] >= m_resolution[0])||
		(voxel[1] < 0)||(voxel[1] >= m_resolution[1])||
		(voxel[2] < 0)||(voxel[2] >= m_resolution[2]))
		return false;
	return true;
}


PNSolution::PNSolution(int order, const V3i& resolution, const Box3d& bound, const std::complex<double> *data):
	m_order(order),
	m_numCoeffs((order + 1) * (order + 1)),
	m_resolution(resolution),
	m_bound(bound)
{
	m_data = std::vector<std::complex<double>>(data, data+m_resolution[0]*m_resolution[1]*m_resolution[2]*m_numCoeffs);
	m_extend = m_bound.getExtents();
}

PNSolution::PNSolution(const std::string& filename)
{
	std::ifstream file( filename.c_str(), std::ofstream::in|std::ofstream::binary );

	file.read((char *)&m_order, sizeof(int));
	file.read((char *)&m_resolution, sizeof(V3i));
	file.read((char *)&m_bound, sizeof(Box3d));
	m_numCoeffs = (m_order + 1) * (m_order + 1);
	m_extend = m_bound.getExtents();
	m_data.resize(m_resolution[0]*m_resolution[1]*m_resolution[2]*m_numCoeffs);
	file.read((char *)m_data.data(), sizeof(std::complex<double>)*m_data.size());
}

/*
PNSolution PNSolution::read( const std::string& filename )
{
	std::ifstream file( filename.c_str(), std::ofstream::in|std::ofstream::binary );
	int order=0;
	V3i resolution;
	Box3d bound;

	file.read((char *)&order, sizeof(int));
	file.read((char *)&resolution, sizeof(V3i));
	file.read((char *)&bound, sizeof(Box3d));

	int numCoeffs = (order + 1) * (order + 1);
	std::vector<std::complex<double>> data(resolution[0]*resolution[1]*resolution[2]*numCoeffs);
	file.read((char *)data.data(), sizeof(std::complex<double>)*data.size());

	return PNSolution(order, resolution, bound, data.data());
}
*/

void PNSolution::save( const std::string& filename )
{
	std::ofstream file( filename.c_str(), std::ofstream::out|std::ofstream::binary|std::ios::trunc );
	file.write( (const char*)&m_order, sizeof(int) );
	file.write( (const char*)&m_resolution, sizeof(V3i) );
	file.write( (const char*)&m_bound, sizeof(Box3d) );
	file.write( (const char*)m_data.data(), sizeof(std::complex<double>)*m_data.size() );
}

int PNSolution::getOrder()const
{
	return m_order;
}

int PNSolution::getNumCoeffs()const
{
	return m_numCoeffs;
}

V3i PNSolution::getResolution()const
{
	return m_resolution;
}

P3d PNSolution::getBoundMin()const
{
	return m_bound.min;
}

P3d PNSolution::getBoundMax()const
{
	return m_bound.max;
}

std::complex<double>* PNSolution::data()
{
	return m_data.data();
}
