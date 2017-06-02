#pragma once

#define EIGEN_NO_DEBUG 

#if defined(_MSC_VER)
/* Disable some warnings on MSVC++ */
#pragma warning(disable : 4307) /*  "integral constant overflow" in Eigen */
#define WIN32_LEAN_AND_MEAN     /* Don't ever include MFC on Windows */
#define NOMINMAX                /* Don't override min/max */
#endif

/* Include the basics needed by any Nori file */
#include <stdexcept> // runtime_error
#include <iostream>
#include <algorithm>
#include <vector>
#include <Eigen/Core>
#include <stdint.h>
#include <math.h>
#include <util/string.h>

#ifndef M_PI
	#define M_PI    3.14159265358979323846
#endif

#ifndef M_PI_2
	#define M_PI_2  1.57079632679489661923 // pi/2
#endif


#if defined(__APPLE__)
#define PLATFORM_MACOS
#elif defined(__linux__)
#define PLATFORM_LINUX
#elif defined(WIN32)
#define PLATFORM_WINDOWS
#endif

/* "Ray epsilon": relative error threshold for ray intersection computations */
#define Epsilon 1e-4f

/* A few useful constants */
#define INV_PI       0.31830988618379067154f
#define INV_TWOPI    0.15915494309189533577f
#define INV_FOURPI   0.07957747154594766788f
#define SQRT_TWO     1.41421356237309504880f
#define INV_SQRT_TWO 0.70710678118654752440f

/* Optimization-related macros */
#if defined(__GNUC__)
#define EXPECT_TAKEN(a)        __builtin_expect(a, true)
#define EXPECT_NOT_TAKEN(a)    __builtin_expect(a, false)
#if defined(__linux) 
#define __restrict             __restrict__
#endif
#else
#define EXPECT_TAKEN(a)        a
#define EXPECT_NOT_TAKEN(a)    a
#endif

/* MSVC is missing a few C99 functions */
#if defined(_MSC_VER)
	/// No nextafterf()! -- an implementation is provided in support_win32.cpp
	extern float nextafterf(float x, float y);
#endif

#if !defined(_GNU_SOURCE)
	/// Emulate sincosf using sinf() and cosf()
	inline void sincosf(float theta, float *_sin, float *_cos) {
		*_sin = sinf(theta);
		*_cos = cosf(theta);
	}
#endif


/* Forward declarations */
template <typename Scalar, int Dimension>  struct TVector;
template <typename Scalar, int Dimension>  struct TPoint;
template <typename Point, typename Vector> struct TRay;
template <typename Point>                  struct TBoundingBox;

/* Basic Nori data structures (vectors, points, rays, bounding boxes,
   kd-trees) are oblivious to the underlying data type and dimension.
   The following list of typedefs establishes some convenient aliases
   for specific types. */
typedef TVector<float, 1>       Vector1f;
typedef TVector<float, 2>       Vector2f;
typedef TVector<float, 3>       Vector3f;
typedef TVector<float, 4>       Vector4f;
typedef TVector<double, 1>      Vector1d;
typedef TVector<double, 2>      Vector2d;
typedef TVector<double, 3>      Vector3d;
typedef TVector<double, 4>      Vector4d;
typedef TVector<int, 1>         Vector1i;
typedef TVector<int, 2>         Vector2i;
typedef TVector<int, 3>         Vector3i;
typedef TVector<int, 4>         Vector4i;
typedef TPoint<float, 1>        Point1f;
typedef TPoint<float, 2>        Point2f;
typedef TPoint<float, 3>        Point3f;
typedef TPoint<float, 4>        Point4f;
typedef TPoint<double, 1>       Point1d;
typedef TPoint<double, 2>       Point2d;
typedef TPoint<double, 3>       Point3d;
typedef TPoint<double, 4>       Point4d;
typedef TPoint<int, 1>          Point1i;
typedef TPoint<int, 2>          Point2i;
typedef TPoint<int, 3>          Point3i;
typedef TPoint<int, 4>          Point4i;
typedef TBoundingBox<Point1f>   BoundingBox1f;
typedef TBoundingBox<Point2f>   BoundingBox2f;
typedef TBoundingBox<Point3f>   BoundingBox3f;
typedef TBoundingBox<Point4f>   BoundingBox4f;
typedef TBoundingBox<Point1d>   BoundingBox1d;
typedef TBoundingBox<Point2d>   BoundingBox2d;
typedef TBoundingBox<Point3d>   BoundingBox3d;
typedef TBoundingBox<Point4d>   BoundingBox4d;
typedef TBoundingBox<Point1i>   BoundingBox1i;
typedef TBoundingBox<Point2i>   BoundingBox2i;
typedef TBoundingBox<Point3i>   BoundingBox3i;
typedef TBoundingBox<Point4i>   BoundingBox4i;
typedef TRay<Point2f, Vector2f> Ray2f;
typedef TRay<Point3f, Vector3f> Ray3f;
typedef TRay<Point2d, Vector2d> Ray2d;
typedef TRay<Point3d, Vector3d> Ray3d;


typedef Vector2i V2i;
typedef Vector2d V2d;

typedef Point2f P2f;
typedef Point2d P2d;
typedef Point2i P2i;

typedef Vector3f V3f;
typedef Vector3d V3d;
typedef Vector3i V3i;

typedef Vector4f V4f;
typedef Vector4d V4d;
typedef Vector4i V4i;

typedef Point3f P3f;
typedef Point3d P3d;

typedef Point4d P4d;


typedef BoundingBox2f Box2f;
typedef BoundingBox2d Box2d;
typedef BoundingBox2i Box2i;

typedef BoundingBox3f Box3f;
typedef BoundingBox3d Box3d;
typedef BoundingBox3i Box3i;

typedef Eigen::Matrix<float, 3, 3> M33f;
typedef Eigen::Matrix<double, 3, 3> M33d;

typedef Eigen::Matrix<float, 4, 4> M44f;
typedef Eigen::Matrix<double, 4, 4> M44d;

/// Import cout, cerr, endl for debugging purposes
using std::cout;
using std::cerr;
using std::endl;

// identifies what we are integrating
// this is important for non-symmetric bsdfs
enum ETransportMode
{
	/// Radiance transport (paths start at camera and we want to integrate radiance)
	ERadiance = 0,
	/// Importance transport (paths start at light source and we want to integrate importance)
	EImportance = 1,
};


/// Measures associated with probability distributions
enum EMeasure {
	EUnknownMeasure = 0,
	ESolidAngle,
	EDiscrete
};



template<typename T>
void sincos(T theta, T *_sin, T *_cos)
{
	*_sin = std::sin(theta);
	*_cos = std::cos(theta);
}

template<typename T>
inline double safe_acos( const T& in )
{
	return std::acos( std::max( T(-1.0), std::min(T(1.0), in)) );
}

template<typename T>
inline double safe_asin( const T& in )
{
	return std::asin( std::max( T(-1.0), std::min(T(1.0), in)) );
}

template<typename T>
T mean( const std::vector<T>& samples )
{
	T m = 0.0;

	int numSamples = int(samples.size());
	for( int i=0;i<numSamples;++i )
		m += (samples[i]-m)/T(i+1);

	return m;
}

// returns sample variance
template<typename T>
T variance( const std::vector<T>& samples, T _mean )
{
	int n = int(samples.size());
	T sum = T(0.0);
	for( int i=0;i<n;++i )
	{
		T e = samples[i] - _mean;
		sum += e*e;
	}
	return T(1.0)/T(n-1)*sum;
}

// returns sample variance
template<typename T>
T variance( const std::vector<T>& samples )
{
	int n = int(samples.size());
	T m = mean<T>( samples );
	T sum = T(0.0);
	for( int i=0;i<n;++i )
	{
		T e = samples[i] - m;
		sum += e*e;
	}
	return T(1.0)/T(n-1)*sum;
}


template<typename T>
T covariance( const std::vector<T>& samples_a, const std::vector<T>& samples_b )
{
	if( samples_a.size() != samples_b.size() )
	{
		throw std::runtime_error("covariance requires equal length sample vectors");
		return 0.0;
	}
	int numSamples = int(samples_a.size());

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat(numSamples, 2);
	for( int i=0;i<numSamples;++i )
	{
		mat(i, 0) = samples_a[i];
		mat(i, 1) = samples_b[i];
	}
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> centered = mat.rowwise() - mat.colwise().mean();
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> cov = (centered.adjoint() * centered) / double(mat.rows() - 1);

	return cov(1, 0);
}

template<typename T>
T correlation( const std::vector<T>& samples_a, const std::vector<T>& samples_b )
{
	T cov_a_b = covariance<T>(samples_a, samples_b);
	T var_a = variance<T>(samples_a);
	T var_b = variance<T>(samples_b);
	return cov_a_b/std::sqrt( var_a*var_b );
}

//// Convert radians to degrees
inline float radToDeg(float value) { return value * (180.0f / float(M_PI)); }

/// Convert degrees to radians
inline float degToRad(float value) { return value * (float(M_PI) / 180.0f); }


template<typename T>
T clamp(T value, T min, T max) {
	if (value < min)
		return min;
	else if (value > max)
		return max;
	else return value;
}

template<typename T, typename R>
inline T lerp( T x0, T x1, R t )
{
	return T(x0*((R)(1.0)-t) + x1*t);
}

template<typename T>
inline T sign( T x )
{
	return (x < T(0)) ? T(-1) : T(1);
}

template<typename T>
inline T signum( T x )
{
	return std::copysign( T(1), x);
}

template<typename T>
void coordinateSystem(const TVector<T, 3> &a, TVector<T, 3> &b, TVector<T, 3> &c)
{
	if (std::abs(a.x()) > std::abs(a.y())) {
		float invLen = 1.0f / std::sqrt(a.x() * a.x() + a.z() * a.z());
		c = TVector<T, 3>(a.z() * invLen, 0.0f, -a.x() * invLen);
	} else {
		float invLen = 1.0f / std::sqrt(a.y() * a.y() + a.z() * a.z());
		c = TVector<T, 3>(0.0f, a.z() * invLen, -a.y() * invLen);
	}
	b = c.cross(a);
}

//typedef coordinateSystem<float> coordinateSystemf;

/// take the absolute of the dot product between two vectors
extern float absDot( const Vector3f &a, const Vector3f &b );

/// Uniformly sample a vector on the unit sphere with respect to solid angles
extern Vector3f squareToUniformSphere(const Point2f &sample);

/// Uniformly sample a vector on the unit hemisphere with respect to solid angles
extern Vector3f squareToUniformHemisphere(const Point2f &sample);

/// Uniformly sample a vector on the unit hemisphere with respect to projected solid angles
extern Vector3f squareToCosineHemisphere(const Point2f &sample);

/// Density of \ref squareToCosineHemiphere with respect to solid angles
extern float squareToCosineHemispherePdf(const Vector3f &d);

/// Uniformly sample a vector on a 2D disk
extern Point2f squareToUniformDisk(const Point2f &sample);



/// Convert an uniformly distributed square sample into barycentric coordinates
extern Point2f squareToUniformTriangle(const Point2f &sample);

/// Compute a direction for the given coordinates in acos()l coordinates
template<typename T>
TVector<T, 3> sphericalDirection(T theta, T phi)
{
	T sinTheta, cosTheta, sinPhi, cosPhi;

	sincos(theta, &sinTheta, &cosTheta);
	sincos(phi, &sinPhi, &cosPhi);

	return TVector<T, 3>(
	   sinTheta * cosPhi,
	   sinTheta * sinPhi,
	   cosTheta
	);
}


/// Compute a direction for the given coordinates in spherical coordinates
template<typename T>
TPoint<T, 2> sphericalCoordinates(const TVector<T, 3> &v)
{
	TPoint<T, 2> result(
		std::acos(v.z()),
		std::atan2(v.y(), v.x())
	);
	if (result.y() < 0)
		result.y() += 2*M_PI;
	return result;
}

// up needs to be normalized
template<typename T>
Eigen::Matrix<T, 4, 4> lookAt( const TPoint<T, 3>& origin, const TPoint<T, 3>& target, const TVector<T, 3>& up = TVector<T, 3>(0.0, 1.0, 0.0) )
{
	TVector<T, 3> dir = (target - origin).normalized();
	TVector<T, 3> left = up.normalized().cross(dir);
	TVector<T, 3> newUp = dir.cross(left);

	Eigen::Matrix<T, 4, 4> trafo;
	trafo << left, newUp, dir, origin,
			  0, 0, 0, 1;

	return trafo;
	//m_transform = Eigen::Affine3f(trafo) * m_transform;
}

template<typename T>
T remap( const T &sourceRangeMin, const T &sourceRangeMax, const T &value )
{
	return (value-sourceRangeMin) / (sourceRangeMax - sourceRangeMin);
}

/// Allocate an aligned region of memory
extern void *allocAligned(size_t size);

/// Free an aligned region of memory
extern void freeAligned(void *ptr);

/// Return the number of cores (real and virtual)
extern int getCoreCount();

