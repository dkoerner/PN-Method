#include <math/shvector.h>
#include <math/sph.h>
//#include <mitsuba/core/transform.h>
#include <boost/math/special_functions/factorials.hpp>

double *SHVector::m_normalization = NULL;
/*
SHVector::SHVector(Stream *stream) {
    m_bands = stream->readInt();
    unsigned int size = m_bands*m_bands;
    m_coeffs.resize(size);
    for (size_t i=0; i<size; ++i)
		m_coeffs[i] = stream->readdouble();
}


void SHVector::serialize(Stream *stream) const {
    stream->writeInt(m_bands);
    for (size_t i=0; i<(size_t) m_coeffs.size(); ++i)
		stream->writedouble(m_coeffs[i]);
}
*/


double legendreP(int l, int m, double x)
{
	double p_mm = 1;

	if (m > 0) {
		double somx2 = std::sqrt((1 - x) * (1 + x));
		double fact = 1;
		for (int i=1; i<=m; i++) {
			p_mm *= (-fact) * somx2;
			fact += 2;
		}
	}

	if (l == m)
		return p_mm;

	double p_mmp1 = x * (2*m + 1) * p_mm;
	if (l == m+1)
		return p_mmp1;

	double p_ll = 0;
	for (int ll=m+2; ll <= l; ++ll) {
		p_ll = ((2*ll-1)*x*p_mmp1 - (ll+m-1) * p_mm) / (ll-m);
		p_mm = p_mmp1;
		p_mmp1 = p_ll;
	}

	return p_ll;
}

bool SHVector::isAzimuthallyInvariant() const
{
    for (int l=0; l<m_bands; ++l) {
        for (int m=1; m<=l; ++m) {
            if (std::abs(operator()(l, -m)) > Epsilon
             || std::abs(operator()(l, m)) > Epsilon)
                return false;
        }
    }
    return true;
}

double SHVector::eval(double theta, double phi) const
{
	double result = 0;
	double cosTheta = std::cos(theta);
	double *sinPhi = (double *) alloca(sizeof(double)*m_bands),
		  *cosPhi = (double *) alloca(sizeof(double)*m_bands);

    for (int m=0; m<m_bands; ++m) {
        sinPhi[m] = std::sin((m+1) * phi);
        cosPhi[m] = std::cos((m+1) * phi);
    }

    for (int l=0; l<m_bands; ++l) {
        for (int m=1; m<=l; ++m) {
			double L = legendreP(l, m, cosTheta) * normalization(l, m);
            result += operator()(l, -m) * SQRT_TWO * sinPhi[m-1] * L;
            result += operator()(l, m)  * SQRT_TWO * cosPhi[m-1] * L;
        }

        result += operator()(l, 0) * legendreP(l, 0, cosTheta) * normalization(l, 0);
    }
    return result;
}

double SHVector::findMinimum(int res = 32) const
{
	double hExt = (double) M_PI / res, hInt = (2 * (double) M_PI)/(res*2);
	double minimum = std::numeric_limits<double>::infinity();

    for (int i=0; i<=res; ++i) {
		double theta = hExt*i;
        for (int j=0; j<=res*2; ++j) {
			double phi = hInt*j;
            minimum = std::min(minimum, eval(theta, phi));
        }
    }

    return minimum;
}

void SHVector::addOffset(double value) {
	operator()(0, 0) += 2 * value * (double) std::sqrt(M_PI);
}

double SHVector::eval(const Vector &v) const
{
	double result = 0;

	double cosTheta = v.z(), phi = std::atan2(v.y(), v.x());
    if (phi < 0) phi += 2*M_PI;
	double *sinPhi = (double *) alloca(sizeof(double)*m_bands),
		  *cosPhi = (double *) alloca(sizeof(double)*m_bands);

    for (int m=0; m<m_bands; ++m) {
        sinPhi[m] = std::sin((m+1) * phi);
        cosPhi[m] = std::cos((m+1) * phi);
    }
    for (int l=0; l<m_bands; ++l) {
        for (int m=1; m<=l; ++m) {
			double L = legendreP(l, m, cosTheta) * normalization(l, m);
            result += operator()(l, -m) * SQRT_TWO * sinPhi[m-1] * L;
            result += operator()(l, m)  * SQRT_TWO * cosPhi[m-1] * L;
        }

        result += operator()(l, 0) * legendreP(l, 0, cosTheta) * normalization(l, 0);
    }
    return result;
}

double SHVector::evalAzimuthallyInvariant(double theta, double phi) const
{
	double result = 0, cosTheta = std::cos(theta);
    for (int l=0; l<m_bands; ++l)
        result += operator()(l, 0) * legendreP(l, 0, cosTheta) * normalization(l, 0);
    return result;
}

double SHVector::evalAzimuthallyInvariant(const Vector &v) const
{
	double result = 0, cosTheta = v.z();
    for (int l=0; l<m_bands; ++l)
        result += operator()(l, 0) * legendreP(l, 0, cosTheta) * normalization(l, 0);
    return result;
}

void SHVector::normalize()
{
	double correction = 1/(2 * (double) std::sqrt(M_PI)*operator()(0,0));

    for (size_t i=0; i<(size_t) m_coeffs.size(); ++i)
        m_coeffs[i] *= correction;
}


double SHVector::computeNormalization(int l, int m)
{
	//SAssert(m>=0);
    return std::sqrt(
			((2*l+1) * boost::math::factorial<double>(l-m))
		/    (4 * (double) M_PI * boost::math::factorial<double>(l+m)));
}


void SHVector::staticInitialization() {
	m_normalization = new double[SH_NORMTBL_SIZE*(SH_NORMTBL_SIZE+1)/2];
    for (int l=0; l<SH_NORMTBL_SIZE; ++l)
        for (int m=0; m<=l; ++m)
            m_normalization[l*(l+1)/2 + m] = computeNormalization(l, m);
}

void SHVector::staticShutdown() {
    delete[] m_normalization;
    m_normalization = NULL;
}


SHSampler::SHSampler(int bands, int depth) : m_bands(bands), m_depth(depth)
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

std::string SHSampler::toString() const {
    std::ostringstream oss;
    oss << "SHSampler[bands=" << m_bands << ", depth=" << m_depth
        << ", size=" << (m_dataSize*sizeof(double))/1024 << " KiB]";
    return oss.str();
}





double SHSampler::integrate(int depth, int zBlock, int phiBlock, const SHVector &f) const
{
	double result = 0;
	for (int l=0; l<m_bands; ++l)
	{
		for (int m=-l; m<=l; ++m)
		{
			double basisIntegral = m_normalization[I(l, std::abs(m))]*lookupIntegral(depth, zBlock, phiBlock, l, m);
			result += basisIntegral * f(l, m);
		}
	}
	return result;
}



int indexofSmallestElement(double array[], int size)
{
	int index = 0;

	for(int i = 1; i < size; i++)
	{
		if(array[i] < array[index])
			index = i;
	}

	return index;
}

double SHSampler::integrateChilds(int depth, int i, int j, const SHVector &f) const
{
	if(depth<m_depth)
	{
		double q00_1 = std::max(integrate(depth+1, i*2, j*2, f), 0.0);
		double q10_1 = std::max(integrate(depth+1, i*2, j*2+1, f), 0.0);
		double q01_1 = std::max(integrate(depth+1, i*2+1, j*2, f), 0.0);
		double q11_1 = std::max(integrate(depth+1, i*2+1, j*2+1, f), 0.0);
		return q00_1+q10_1+q01_1+q11_1;
	}else
		return std::max(integrate(depth, i, j, f), 0.0);

}

double SHSampler::warp(const SHVector &f, P2d &sample) const
{
    int i = 0, j = 0;
	double integral = 0;
	double integralRoot = 0.0;
	integralRoot = integrate(0, 0, 0, f);

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
		double q[4] = {integrateChilds(depth, i, j, f),
					   integrateChilds(depth, i, j+1, f),
					   integrateChilds(depth, i+1, j, f),
					   integrateChilds(depth, i+1, j+1, f)};

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

SHSampler::~SHSampler() {
    for (int i=0; i<=m_depth; ++i) {
        int res = 1 << i;
        for (int j=0; j<res; ++j) {
            delete[] m_phiMap[i][j];
            delete[] m_legendreMap[i][j];
        }
        delete[] m_phiMap[i];
        delete[] m_legendreMap[i];
    }
    delete[] m_phiMap;
    delete[] m_legendreMap;
    delete[] m_normalization;
}



double *SHSampler::phiIntegrals(double a, double b)
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

double *SHSampler::legendreIntegrals(double a, double b)
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
            Pa[I(l,m)] = legendreP(l, m, a);
            Pb[I(l,m)] = legendreP(l, m, b);
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


Eigen::MatrixXd SHSampler::getBlocks(int depth, const SHVector& f, bool use_org) const
{
	int res = 1 << depth;
	Eigen::MatrixXd	blocks = Eigen::MatrixXd::Zero(res, res);

	for( int block_i=0;block_i<res;++block_i )
		for( int block_j=0;block_j<res;++block_j )
			blocks.coeffRef(block_i, block_j) = integrate(depth, block_i, block_j, f);

	return blocks;
}

