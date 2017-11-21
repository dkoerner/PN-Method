#include <PNSolution.h>
#include <math/sph.h>

#include <iostream>
#include <fstream>


#include <boost/math/special_functions/factorials.hpp>



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
	/*
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
	double theta, phi;
	sphericalCoordinates(direction, theta, phi);


	double result = 0.0;
	int coeff_index = 0;
	for( int l=0;l<=m_order;++l )
		for( int m=-l;m<=l;++m,++coeff_index )
		{
			//lerp coefficient
			double coeff = 0.0;
			coeff += m_data[voxel_indices[0] + coeff_index]*(1.0-tx)*(1.0-ty)*(1.0-tz);
			coeff += m_data[voxel_indices[1] + coeff_index]*tx*(1.0-ty)*(1.0-tz);
			coeff += m_data[voxel_indices[2] + coeff_index]*(1.0-tx)*ty*(1.0-tz);
			coeff += m_data[voxel_indices[3] + coeff_index]*tx*ty*(1.0-tz);
			coeff += m_data[voxel_indices[4] + coeff_index]*(1.0-tx)*(1.0-ty)*tz;
			coeff += m_data[voxel_indices[5] + coeff_index]*tx*(1.0-ty)*tz;
			coeff += m_data[voxel_indices[6] + coeff_index]*(1.0-tx)*ty*tz;
			coeff += m_data[voxel_indices[7] + coeff_index]*tx*ty*tz;

			// now do the SH accumulation
			result += coeff*sph::basis_real(l, m, theta, phi);
		}

	// NB: here we apply a normalization. This is required to have the intergal over
	// the sphere match the zero coefficient
	return result/std::sqrt(4.0*M_PI);
	*/
	Eigen::VectorXd coeffs = evalCoefficients(pWS);

	if( !m_coeffs_filter.empty() )
		sph::convolve(coeffs.data(), m_coeffs_filter.data(), m_order);

	double theta, phi;
	sphericalCoordinates(direction, theta, phi);

	double result = sph::eval(theta, phi, coeffs.data(), m_order);
	return result/std::sqrt(4.0*M_PI);
}

double PNSolution::evalCoefficient(const P3d &pWS, int coeff_index)const
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

	//lerp coefficient
	double coeff = 0.0;
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

void PNSolution::evalCoefficients(const P3d &pWS, double* result)const
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



	int coeff_index = 0;
	for( int l=0;l<=m_order;++l )
		for( int m=-l;m<=l;++m,++coeff_index )
		{
			double coeff = 0.0;

			//lerp coefficient
			coeff += m_data[voxel_indices[0] + coeff_index]*(1.0-tx)*(1.0-ty)*(1.0-tz);
			coeff += m_data[voxel_indices[1] + coeff_index]*tx*(1.0-ty)*(1.0-tz);
			coeff += m_data[voxel_indices[2] + coeff_index]*(1.0-tx)*ty*(1.0-tz);
			coeff += m_data[voxel_indices[3] + coeff_index]*tx*ty*(1.0-tz);
			coeff += m_data[voxel_indices[4] + coeff_index]*(1.0-tx)*(1.0-ty)*tz;
			coeff += m_data[voxel_indices[5] + coeff_index]*tx*(1.0-ty)*tz;
			coeff += m_data[voxel_indices[6] + coeff_index]*(1.0-tx)*ty*tz;
			coeff += m_data[voxel_indices[7] + coeff_index]*tx*ty*tz;

			result[coeff_index] = coeff;
		}
}

Eigen::VectorXd PNSolution::evalCoefficients(const P3d &pWS)const
{
	Eigen::VectorXd result(m_numCoeffs);
	evalCoefficients(pWS, result.data());
	return result;
}


V3d PNSolution::sample( const P3d& pWS, double& pdf, RNGd& rng )const
{
	// evaluate coefficients at given world-space position
	Eigen::VectorXd coeffs = evalCoefficients(pWS);
	// generate sample in unit 2d space
	P2d sample( rng.next1D(), rng.next1D() );
	// warp it into distribution given by SH function
	pdf = m_shsampler.warp( coeffs.data(), sample );
	// and transform from spherical angles to cartesian vector
	return sphericalDirection( sample[0], sample[1] );
}

V3d PNSolution::sample( const P3d& pWS, double& pdf, const P2d& sample )const
{
	P2d sample2 = sample;
	// evaluate coefficients at given world-space position
	Eigen::VectorXd coeffs = evalCoefficients(pWS);
	// warp it into distribution given by SH function
	pdf = m_shsampler.warp( coeffs.data(), sample2 );
	// and transform from spherical angles to cartesian vector
	return sphericalDirection( sample2[0], sample2[1] );
}



double PNSolution::pdf( const P3d& pWS, const V3d& direction )const
{
	// evaluate coefficients at given world-space position
	Eigen::VectorXd coeffs = evalCoefficients(pWS);
	return m_shsampler.pdf( coeffs.data(), sphericalCoordinates(direction) );
}

Eigen::MatrixXd PNSolution::getBlocks( const P3d& pWS, int depth )const
{
	Eigen::VectorXd coeffs = evalCoefficients(pWS);
	//return m_shsampler.getBlocks(m_shsampler.m_depth, coeffs.data());
	return m_shsampler.getBlocks(depth, coeffs.data());
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


PNSolution::PNSolution(int order, const V3i& resolution, const Box3d& bound, const double *data):
	m_order(order),
	m_numCoeffs((order + 1) * (order + 1)),
	m_resolution(resolution),
	m_bound(bound),
	m_shsampler(order+1)
{
	sph::staticInit();

	m_data = std::vector<double>(data, data+m_resolution[0]*m_resolution[1]*m_resolution[2]*m_numCoeffs);
	m_extend = m_bound.getExtents();
}

PNSolution::PNSolution(const std::string& filename):
	m_shsampler(1)
{
	sph::staticInit();

	std::ifstream file( filename.c_str(), std::ofstream::in|std::ofstream::binary );
	file.read((char *)&m_order, sizeof(int));
	file.read((char *)&m_resolution, sizeof(V3i));
	file.read((char *)&m_bound, sizeof(Box3d));
	m_numCoeffs = (m_order + 1) * (m_order + 1);
	m_extend = m_bound.getExtents();
	m_data.resize(m_resolution[0]*m_resolution[1]*m_resolution[2]*m_numCoeffs);
	file.read((char *)m_data.data(), sizeof(double)*m_data.size());

	m_shsampler = SHSampler(m_order+1);





}

void PNSolution::setFilter( const std::string& name, double width )
{
	if( name == "none" )
	{
		m_coeffs_filter.clear();
		return;
	}

	// filter kernel test ---------
	std::function<double(double, double)> hanning =
	[=]( double theta, double phi )
	{
		double w = width;
		double value = M_PI*theta/w;
		return (1.0+std::cos(value))*0.5;
	};
	std::function<double(double, double)> lanzcos =
	[=]( double theta, double phi )
	{
		double w = width;
		double value = M_PI*theta/w;
		if(std::abs(value) < 1.0e-6)
			return 1.0;
		return std::sin(value)/value;
	};

	m_coeffs_filter = std::vector<double>(m_numCoeffs, 0.0);

	std::function<double(double, double)> filter;
	if( name == "hanning" )
		filter = hanning;
	else
	if( name == "lanzcos" )
		filter = lanzcos;
	else
		throw std::runtime_error("PNSolution::setFilter: unknown filter");

	sph::project( filter, m_order, 512, m_coeffs_filter.data() );
}

void PNSolution::save( const std::string& filename )
{
	std::ofstream file( filename.c_str(), std::ofstream::out|std::ofstream::binary|std::ios::trunc );
	file.write( (const char*)&m_order, sizeof(int) );
	file.write( (const char*)&m_resolution, sizeof(V3i) );
	file.write( (const char*)&m_bound, sizeof(Box3d) );
	file.write( (const char*)m_data.data(), sizeof(double)*m_data.size() );
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

double* PNSolution::data()
{
	return m_data.data();
}









// PNSolution::SHSampler -------------------------------

PNSolution::SHSampler::SHSampler(int bands, int depth) : m_bands(bands), m_depth(depth)
{
	m_phiMap = new double**[depth+1];
	m_legendreMap = new double**[depth+1];
	m_normalization = new double[m_bands*(m_bands+1)/2];
	m_dataSize = m_bands*(m_bands+1)/2;
	//Assert(depth >= 1);

	for (int i=0; i<=depth; ++i)
	{
		int res = 1 << i;
		double zStep  = -2 / (double) res;
		double phiStep = 2 * (double) M_PI / (double) res;
		m_phiMap[i] = new double*[res];
		m_legendreMap[i] = new double*[res];

		for (int j=0; j<res; ++j)
		{
			m_phiMap[i][j] = phiIntegrals(phiStep*j, phiStep*(j+1));
			m_legendreMap[i][j] = legendreIntegrals(1+zStep*j, 1+zStep*(j+1));
		}
	}

	for (int l=0; l<m_bands; ++l)
	{
		for (int m=0; m<=l; ++m)
		{
			double normFactor = boost::math::tgamma_delta_ratio(
				(double) (l - m + 1), (double) (2 * m), boost::math::policies::policy<>());
			normFactor = std::sqrt(normFactor * (2 * l + 1) / (4 * (double) M_PI));
			if (m != 0)
				normFactor *= SQRT_TWO;
			m_normalization[I(l, m)] = normFactor;
		}
	}
}


double PNSolution::SHSampler::integrate(int depth, int zBlock, int phiBlock, const double* coeffs) const
{
	double result = 0;
	for (int l=0; l<m_bands; ++l)
	{
		for (int m=-l; m<=l; ++m)
		{
			double basisIntegral = m_normalization[I(l, std::abs(m))]*lookupIntegral(depth, zBlock, phiBlock, l, m);
			result += basisIntegral * sph::get(l, m, coeffs);
		}
	}
	return result;
}



int PNSolution::SHSampler::indexofSmallestElement(double array[], int size)
{
	int index = 0;

	for(int i = 1; i < size; i++)
	{
		if(array[i] < array[index])
			index = i;
	}

	return index;
}

double PNSolution::SHSampler::integrateChilds(int depth, int i, int j, const double *coeffs) const
{
	if(depth<m_depth)
	{
		double q00_1 = std::max(integrate(depth+1, i*2, j*2, coeffs), 0.0);
		double q10_1 = std::max(integrate(depth+1, i*2, j*2+1, coeffs), 0.0);
		double q01_1 = std::max(integrate(depth+1, i*2+1, j*2, coeffs), 0.0);
		double q11_1 = std::max(integrate(depth+1, i*2+1, j*2+1, coeffs), 0.0);
		return q00_1+q10_1+q01_1+q11_1;
	}else
		return std::max(integrate(depth, i, j, coeffs), 0.0);

}

double PNSolution::SHSampler::warp(const double *coeffs, P2d &sample) const
{
	int i = 0, j = 0;
	double integral = 0;
	double integralRoot = 0.0;
	integralRoot = integrate(0, 0, 0, coeffs);

	for (int depth = 1; depth <= m_depth; ++depth)
	{
		/*
		// Original implementation which causes problems with negative areas...
		// Do not sample negative areas
		double q00_1 = std::max(integrate(depth, i, j, f), (double) 0);
		double q10_1 = std::max(integrate(depth, i, j+1, f), (double) 0);
		double q01_1 = std::max(integrate(depth, i+1, j, f), (double) 0);
		double q11_1 = std::max(integrate(depth, i+1, j+1, f), (double) 0);
		*/
		double q[4] = {integrateChilds(depth, i, j, coeffs),
					   integrateChilds(depth, i, j+1, coeffs),
					   integrateChilds(depth, i+1, j, coeffs),
					   integrateChilds(depth, i+1, j+1, coeffs)};

		double q00 = q[0];
		double q10 = q[1];
		double q01 = q[2];
		double q11 = q[3];

		double z1 = q00 + q10, z2 = q01 + q11, phi1, phi2;
		double zNorm = (double) 1 / (z1+z2);
		z1 *= zNorm; z2 *= zNorm;

		if (sample.x() < z1)
		{
			sample.x() /= z1;
			phi1 = q00; phi2 = q10;
			i <<= 1;
		}else
		{
			sample.x() = (sample.x() - z1) / z2;
			phi1 = q01; phi2 = q11;
			i = (i+1) << 1;
		}

		double phiNorm = (double) 1 / (phi1+phi2);
		double phi1Norm = phi1*phiNorm, phi2Norm = phi2*phiNorm;

		if (sample.y() <= phi1Norm)
		{
			sample.y() /= phi1Norm;
			j <<= 1;
			integral = phi1;
		}else
		{
			sample.y() = (sample.y() - phi1Norm) / phi2Norm;
			j = (j+1) << 1;
			integral = phi2;
		}
	}

	double zStep = -2 / (double) (1 << m_depth);
	double phiStep = 2 * (double) M_PI / (double) (1 << m_depth);
	i >>= 1; j >>= 1;

	double z = 1 + zStep * i + zStep * sample.x();
	sample.x() = std::acos(z);
	sample.y() = phiStep * j + phiStep * sample.y();

	// PDF of sampling the mip-map bin
	double pdfBin = integral/integralRoot;

	// Density within the bin
	double density = -1/(zStep*phiStep);

	return density*pdfBin;
}


double PNSolution::SHSampler::pdf( const double *coeffs, const P2d& sample )const
{
	int res = 1 << m_depth;
	double integralRoot = integrate(0, 0, 0, coeffs);
	double zStep = -2 / (double) (res);
	double phiStep = 2 * (double) M_PI / (double) (res);

	// find bin which has been selected by sample
	double u = (std::cos(sample.x())-1.0)/(-2.0);
	double v = sample.y()/(2.0*M_PI);
	int i = int(std::floor(u*(res-1)));
	int j = int(std::floor(v*(res-1)));
	double integral = integrate(m_depth, i, j, coeffs);

	// PDF of sampling the mip-map bin
	double pdfBin = integral/integralRoot;

	// Density within the bin
	double density = -1/(zStep*phiStep);

	return density*pdfBin;
}



PNSolution::SHSampler::~SHSampler()
{
	// TODO: for some reason the commented code crashes when closing the program
	for (int i=0; i<=m_depth; ++i)
	{
		int res = 1 << i;
		for (int j=0; j<res; ++j)
		{
			//delete[] m_phiMap[i][j];
			//delete[] m_legendreMap[i][j];
		}
		//delete[] m_phiMap[i];
		//delete[] m_legendreMap[i];
	}
	//delete[] m_phiMap;
	//delete[] m_legendreMap;
	//delete[] m_normalization;
}



double *PNSolution::SHSampler::phiIntegrals(double a, double b)
{
	double *sinPhiA = new double[m_bands+1];
	double *sinPhiB = new double[m_bands+1];
	double *cosPhiA = new double[m_bands+1];
	double *cosPhiB = new double[m_bands+1];
	double *result = new double[2*m_bands+1];
	m_dataSize += 2*m_bands+1;

	cosPhiA[0] = 1; sinPhiA[0] = 0;
	cosPhiB[0] = 1; sinPhiB[0] = 0;
	cosPhiA[1] = std::cos(a);
	sinPhiA[1] = std::sin(a);
	cosPhiB[1] = std::cos(b);
	sinPhiB[1] = std::sin(b);

	for (int m=2; m<=m_bands; ++m)
	{
		sinPhiA[m] = 2*sinPhiA[m-1]*cosPhiA[1] - sinPhiA[m-2];
		sinPhiB[m] = 2*sinPhiB[m-1]*cosPhiB[1] - sinPhiB[m-2];

		cosPhiA[m] = 2*cosPhiA[m-1]*cosPhiA[1] - cosPhiA[m-2];
		cosPhiB[m] = 2*cosPhiB[m-1]*cosPhiB[1] - cosPhiB[m-2];
	}

	for (int m=-m_bands; m<=m_bands; ++m)
	{
		if (m == 0)
			result[P(m)] = b-a;
		else if (m > 0)
			result[P(m)] = (sinPhiB[m]-sinPhiA[m])/m;
		else
			result[P(m)] = (cosPhiB[-m]-cosPhiA[-m])/m;
	}

	delete[] sinPhiA;
	delete[] sinPhiB;
	delete[] cosPhiA;
	delete[] cosPhiB;
	return result;
}

double *PNSolution::SHSampler::legendreIntegrals(double a, double b)
{
	double *P = new double[m_bands*(m_bands+1)/2];
	m_dataSize += m_bands*(m_bands+1)/2;

	P[I(0, 0)] = b-a;

	if (m_bands == 1)
		return P;

	double *Pa = new double[m_bands*(m_bands+1)/2];
	double *Pb = new double[m_bands*(m_bands+1)/2];

	for (int l=0; l<m_bands; ++l)
	{
		for (int m=0; m<=l; ++m)
		{
			Pa[I(l,m)] = sph::legendreP(l, m, a);
			Pb[I(l,m)] = sph::legendreP(l, m, b);
		}
	}

	P[I(1,0)] = (b*b - a*a)/2;
	P[I(1,1)] = .5f * (-b*std::sqrt(1-b*b) - std::asin(b) + a*std::sqrt(1-a*a) + std::asin(a));

	for (int l=2; l<m_bands; ++l)
	{
		for (int m=0; m<=l-2; ++m)
		{
			double ga = (2*l-1)*(1-a*a) * Pa[I(l-1,m)];
			double gb = (2*l-1)*(1-b*b) * Pb[I(l-1,m)];
			P[I(l, m)] = ((l-2)*(l-1+m)*P[I(l-2, m)]-gb+ga)/((l+1)*(l-m));
		}

		P[I(l, l-1)] = (2*l-1)/(double)(l+1) * ((1-a*a)*Pa[I(l-1, l-1)] - (1-b*b)*Pb[I(l-1, l-1)]);
		P[I(l, l)] = 1/(double)(l+1) * (l*(2*l-3)*(2*l-1) * P[I(l-2, l-2)] + b*Pb[I(l,l)] - a*Pa[I(l, l)]);
	}

	delete[] Pa;
	delete[] Pb;

	return P;
}




/*
double PNSolution::test_real_conversion( int l, int m, double theta, double phi )
{
	if (m < 0)
		return (std::complex<double>(0.0, 1.0)/std::sqrt(2.0)*SHVector::csp(m)*sph::sph_basis(l, m, theta, phi)-
				std::complex<double>(0.0, 1.0)/std::sqrt(2.0)*sph::sph_basis(l, -m, theta, phi) ).real();
	else
	if(m == 0)
		return (sph::sph_basis(l, 0, theta, phi)).real();
	else
	if(m > 0)
		return (1.0/std::sqrt(2.0)*SHVector::csp(m)*sph::sph_basis(l, -m, theta, phi)+
				1.0/std::sqrt(2.0)*sph::sph_basis(l, m, theta, phi) ).real();
	throw std::runtime_error("asdsad");
}
*/

/*
void PNSolution::test()
{
	P3d pWS(1.1, 1.0, 1.0);

	double theta = M_PI*0.4;
	double phi = 2.0*M_PI*0.4;
	V3d d = sphericalDirection(theta, phi);


	Eigen::VectorXd Y_real_vector(m_numCoeffs);
	Eigen::VectorXcd Y_complex_vector(m_numCoeffs);
	for( int l=0;l<=m_order;++l )
		for( int m=-l;m<=l;++m )
		{
			double sph_real_basis = SHVector::sph_basis_real(l, m, theta, phi);
			std::complex<double> sph_complex_basis = sph::sph_basis(l, m, theta, phi);
			Y_real_vector.coeffRef(SHVector::index(l, m)) = sph_real_basis;
			Y_complex_vector.coeffRef(SHVector::index(l, m)) = sph_complex_basis;
		}

	Eigen::VectorXd Y_real_vector_check = (this->m_complex2RealConversionMatrix*Y_complex_vector).real();

	std::cout << "conversion: " << std::endl;
	for( int i=0;i<m_numCoeffs;++i )
		std::cout << Y_real_vector.coeffRef(i) << " " << Y_real_vector_check.coeffRef(i) << std::endl;

}
*/

Eigen::MatrixXd PNSolution::SHSampler::getBlocks(int depth, const double* coeffs) const
{
	int res = 1 << depth;
	Eigen::MatrixXd	blocks = Eigen::MatrixXd::Zero(res, res);

	for( int block_i=0;block_i<res;++block_i )
		for( int block_j=0;block_j<res;++block_j )
			blocks.coeffRef(block_i, block_j) = integrate(depth, block_i, block_j, coeffs);

	return blocks;
}


