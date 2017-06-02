#pragma once

#include <cmath>
#include <math.h>
#include <ctgmath>

#include <random>
#include <sstream>
#include <vector>
#include <algorithm>

#include <math/vector.h>



template<typename T>
struct RNG
{
	RNG(unsigned long seed=123)
	{
		m_rng.seed(seed);
		m_distr = std::uniform_real_distribution<T>();
	}

	T next1D()
	{
		return m_distr(m_rng);
	}

	std::mt19937& generator()
	{
		return m_rng;
	}

	const std::mt19937& cgenerator()
	{
		return m_rng;
	}

private:
	std::mt19937                      m_rng;
	std::uniform_real_distribution<T> m_distr;
};

typedef RNG<double> RNGd;







// sampling routines --------------------------

template<typename T>
TVector<T, 3> sampleSphere(RNG<T>& rng)
{
	T r1 = T(rng.next1D());
	T r2 = T(rng.next1D());
	T z = T(1.0) - T(2.0) * r2;
	T r = std::sqrt(std::max(T(0.0), T(1.0) - z*z));
	T theta = T(2.0) * M_PI * r1;
	T sinPhi = sin(theta);
	T cosPhi = cos(theta);
	return TVector<T, 3>(r * cosPhi, r * sinPhi, z);
}


inline double sampleSpherePDF()
{
	return INV_FOURPI;
}

template<typename T>
TVector<T, 3> sampleHemisphereUniform(RNG<T>& rng)
{
	T cosTheta = T(rng.next1D());
	T sinTheta = std::sqrt(std::max((T) 0, 1-cosTheta*cosTheta));

	T sinPhi, cosPhi;
	sincos<T>(T(2.0) * M_PI * rng.next1D(), &sinPhi, &cosPhi);

	return TVector<T, 3>(cosPhi * sinTheta, sinPhi * sinTheta, cosTheta);
}

inline double sampleHemisphereUniformPDF()
{
	return 1.0/(2.0*M_PI);
}


/// Low-distortion concentric square to disk mapping by Peter Shirley (PDF: 1/PI)
//extern Point2f squareToUniformDiskConcentric(const Point2f &sample);

/*
template<typename T>
TPoint<T, 2> squareToUniformDiskConcentric( RNG<T>& rng )
{
	float r1 = 2.0f*sample.x() - 1.0f;
	float r2 = 2.0f*sample.y() - 1.0f;

	Point2f coords;
	if (r1 == 0 && r2 == 0) {
		coords = Point2f(0, 0);
	} else if (r1 > -r2) { // Regions 1/2
		if (r1 > r2)
			coords = Point2f(r1, (M_PI/4.0f) * r2/r1);
		else
			coords = Point2f(r2, (M_PI/4.0f) * (2.0f - r1/r2));
	} else { // Regions 3/4
		if (r1<r2)
			coords = Point2f(-r1, (M_PI/4.0f) * (4.0f + r2/r1));
		else
			coords = Point2f(-r2, (M_PI/4.0f) * (6.0f - r1/r2));
	}

	Point2f result;
	sincosf(coords.y(), &result[1], &result[0]);
	return result*coords.x();
}
*/

/// Low-distortion concentric square to disk mapping by Peter Shirley (PDF: 1/PI)
/// http://psgraphics.blogspot.ch/2011/01/improved-code-for-concentric-map.html
template<typename T>
TPoint<T, 2> sampleUniformDiskConcentric(RNG<T>& rng)
{
	T r1 = 2.0f*rng.next1D() - 1.0f;
	T r2 = 2.0f*rng.next1D() - 1.0f;

	TPoint<T, 2> coords;
	if (r1 == 0 && r2 == 0)
	{
		coords = TPoint<T, 2>(0, 0);
	} else if (r1 > -r2) { // Regions 1/2
		if (r1 > r2)
			coords = TPoint<T, 2>(r1, (M_PI/4.0f) * r2/r1);
		else
			coords = TPoint<T, 2>(r2, (M_PI/4.0f) * (2.0f - r1/r2));
	} else { // Regions 3/4
		if (r1<r2)
			coords = TPoint<T, 2>(-r1, (M_PI/4.0f) * (4.0f + r2/r1));
		else
			coords = TPoint<T, 2>(-r2, (M_PI/4.0f) * (6.0f - r1/r2));
	}

	TPoint<T, 2> result(0.0);
	sincos(coords.y(), &result.y(), &result.x());
	return result*coords.x();
}

template<typename T>
TVector<T, 3> sampleCosineHemisphere(RNG<T>& rng)
{
	TPoint<T, 2> p = sampleUniformDiskConcentric(rng);
	T z = std::sqrt(std::max((T) 0,
		T(1.0f) - p.x()*p.x() - p.y()*p.y()));
	return TVector<T, 3>(p.x(), p.y(), z);
}

inline double cosineHemispherePDF( const V3d& d )
{
	return INV_PI * std::abs(d.z());
}


template<typename T>
TVector<T, 3> sampleUniformCone( T cosCutoff, RNG<T>& rng )
{
	T x = rng.next1D();
	T cosTheta = (T(1)-x) + x * cosCutoff;
	T sinTheta = std::sqrt(std::max( T(0.0), T(1) - cosTheta * cosTheta));

	T sinPhi, cosPhi;
	sincos(T(2.0) * M_PI * rng.next1D(), &sinPhi, &cosPhi);

	return TVector<T, 3>(cosPhi * sinTheta,
						 sinPhi * sinTheta,
						 cosTheta);
}

template<typename T>
T sampleUniformConePDF(T cosCutoff)
{
	return INV_TWOPI / (T(1)-cosCutoff);
}

template<typename T>
T normalPDF( T x, T mean, T stddev )
{
	double nom = x-mean;
	double denom = T(2)*stddev*stddev;
	double coeff = T(1)/(stddev*std::sqrt(T(2)*M_PI));
	return coeff*std::exp( -nom*nom/denom );
}


template<typename T>
T sampleTruncatedExponential( T lambda, T max, RNG<T>& rng, T& pdf )
{
	T sigma_t = 1.0/lambda;
	T u = rng.next1D();

	T R = u*(1.0-std::exp(-max*sigma_t));
	T t = -std::log(1-R)*lambda;

	pdf = std::exp( -t*sigma_t )/(lambda*(1.0-std::exp(-max*sigma_t)));

	return t;
}

template<typename T>
T sampleTruncatedExponentialPDF( T lambda, T max, T t )
{
	T sigma_t = 1.0/lambda;
	return std::exp( -t*sigma_t )/(lambda*(1.0-std::exp(-max*sigma_t)));
}





