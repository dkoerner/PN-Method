#pragma once

#include <math/common.h>
#include <Eigen/Dense>


/* ===================================================================
    This file contains a few templates and specializations that
    provide 2/3D points, vectors, and normals over different
    underlying data types. Points, vectors, and normals are distinct
    in Nori, because they transform differently under homogeneous
    coordinate transformations.
 * =================================================================== */

/**
 * \brief Generic N-dimensional vector data structure based on Eigen::Matrix
 */
template <typename _Scalar, int _Dimension> struct TVector : public Eigen::Matrix<_Scalar, _Dimension, 1> {
public:
	enum {
		Dimension = _Dimension
	};

	typedef _Scalar                             Scalar;
	typedef Eigen::Matrix<Scalar, Dimension, 1> Base;
	typedef TVector<Scalar, Dimension>          VectorType;
	typedef TPoint<Scalar, Dimension>           PointType;

	/// Create a new vector with constant component vlaues
	inline TVector(Scalar value = (Scalar) 0) { this->setConstant(value); }

	/// Create a new 2D vector (type error if \c Dimension != 2)
	inline TVector(Scalar x, Scalar y) : Base(x, y) { }
	
	/// Create a new 3D vector (type error if \c Dimension != 3)
	inline TVector(Scalar x, Scalar y, Scalar z) : Base(x, y, z) { }

	/// Create a new 4D vector (type error if \c Dimension != 4)
	inline TVector(Scalar x, Scalar y, Scalar z, Scalar w) : Base(x, y, z, w) { }

	/// Construct a vector from MatrixBase (needed to play nice with Eigen)
	template <typename Derived> inline TVector(const Eigen::MatrixBase<Derived>& p) 
		: Base(p) { }

	/// Assign a vector from MatrixBase (needed to play nice with Eigen)
    template <typename Derived> TVector &operator=(const Eigen::MatrixBase<Derived>& p) {
		this->Base::operator=(p);
		return *this;
    }

//	/// Return a human-readable string summary
//	inline QString toString() const {
//		QString result;
//		for (size_t i=0; i<Dimension; ++i) {
//			result += QString("%1").arg(this->coeff(i));
//			if (i+1 < Dimension)
//				result += ", ";
//		}
//		return "[" + result + "]";
//	}
	std::string toString()const
	{
		std::string result = "(";
		for (size_t i=0; i<Dimension; ++i)
		{
			result += ::toString(this->coeff(i));
			if (i+1 < Dimension)
				result += ", ";
		}
		return result + ")";
	}
};

/**
 * \brief Generic N-dimensional point data structure based on Eigen::Matrix
 */
template <typename _Scalar, int _Dimension> struct TPoint : public Eigen::Matrix<_Scalar, _Dimension, 1> {
public:
	enum {
		Dimension = _Dimension
	};

	typedef _Scalar                             Scalar;
	typedef Eigen::Matrix<Scalar, Dimension, 1> Base;
	typedef TVector<Scalar, Dimension>          VectorType;
	typedef TPoint<Scalar, Dimension>           PointType;

	/// Create a new point with constant component vlaues
	inline TPoint(Scalar value = (Scalar) 0) { this->setConstant(value); }

	/// Create a new 2D point (type error if \c Dimension != 2)
	inline TPoint(Scalar x, Scalar y) : Base(x, y) { }
	
	/// Create a new 3D point (type error if \c Dimension != 3)
	inline TPoint(Scalar x, Scalar y, Scalar z) : Base(x, y, z) { }

	/// Create a new 4D point (type error if \c Dimension != 4)
	inline TPoint(Scalar x, Scalar y, Scalar z, Scalar w) : Base(x, y, z, w) { }

	/// Construct a point from MatrixBase (needed to play nice with Eigen)
	template <typename Derived> inline TPoint(const Eigen::MatrixBase<Derived>& p) 
		: Base(p) { }

	/// Assign a point from MatrixBase (needed to play nice with Eigen)
    template <typename Derived> TPoint &operator=(const Eigen::MatrixBase<Derived>& p) {
		this->Base::operator=(p);
		return *this;
    }

	std::string toString()const
	{
		std::string result = "(";
		for (size_t i=0; i<Dimension; ++i)
		{
			result += ::toString(this->coeff(i));
			if (i+1 < Dimension)
				result += ", ";
		}
		return result + ")";
	}
};

/**
 * \brief 3-dimensional surface normal representation
 */
template<typename T>
struct TNormal3 : public Eigen::Matrix<T, 3, 1> {
public:
	enum {
		Dimension = 3
	};

	typedef T                               Scalar;
	typedef Eigen::Matrix<Scalar, Dimension, 1> Base;
	typedef TVector<Scalar, Dimension>          VectorType;
	typedef TPoint<Scalar, Dimension>           PointType;


	/// Create a new normal with constant component vlaues
	inline TNormal3(Scalar value = 0.0f) { this->setConstant(value); }

	/// Create a new 3D normal 
	inline TNormal3(Scalar x, Scalar y, Scalar z) : Base(x, y, z) { }

	/// Construct a normal from MatrixBase (needed to play nice with Eigen)
	template <typename Derived> inline TNormal3(const Eigen::MatrixBase<Derived>& p)
		: Base(p) { }

	/// Assign a normal from MatrixBase (needed to play nice with Eigen)
	template <typename Derived> TNormal3 &operator=(const Eigen::MatrixBase<Derived>& p) {
		this->Base::operator=(p);
		return *this;
    }

//	/// Return a human-readable string summary
//	inline QString toString() const {
//		return QString("[%1, %2, %3]").arg(coeff(0)).arg(coeff(1)).arg(coeff(2));
//	}
	std::string toString()const
	{
		std::string result = "(";
		for (size_t i=0; i<Dimension; ++i)
		{
			result += ::toString(this->coeff(i));
			if (i+1 < Dimension)
				result += ", ";
		}
		return result + ")";
	}
};


typedef TNormal3<float> Normal3f;
typedef TNormal3<float> N3f;
typedef TNormal3<double> N3d;
typedef TNormal3<double> Normal3d;

//template<typename T>
//inline T dot( const TVector<T, 3>& a, const TVector<T, 3>& b )
//{
//	return a.dot(b);
//}

template<typename T>
inline T dotT( const TVector<T, 3>& a, const TVector<T, 3>& b )
{
	return a.dot(b);
}

template<typename T>
inline TVector<T, 3> cross( const TVector<T, 3>& a, const TVector<T, 3>& b )
{
	return a.cross(b);
}


inline double dot( const V3d& a, const V3d& b )
{
	return a.dot(b);
}
inline float dot2( const V3f& a, const V3f& b )
{
	return a.dot(b);
}
inline float dotf( const V3f& a, const V3f& b )
{
	return a.dot(b);
}

template<typename T>
inline TVector<T, 3> normalizeT( const TVector<T, 3>& v )
{
	return v.normalized();
}


inline V3d normalize( const V3d& v )
{
	return v.normalized();
}

inline V3d normalized( const V3d& v, double& distance )
{
	distance = v.norm();
	return v.normalized();
}


// For a given incident vector I and surface normal N reflect returns the reflection direction calculated as I - 2.0 * dot(N, I) * N. N should be normalized in order to achieve the desired result.
template<typename T>
inline TVector<T, 3> reflect( const TVector<T, 3> &i,
							  const TVector<T, 3> &n )
{
	return -i + T(2.0)*dotT<T>(n,i)*n;
}

template<typename T>
inline TVector<T, 3> refract( const TVector<T, 3> &i,
							  const TVector<T, 3> &n,
							  T eta, // ior_inside/ior_outside
							  T cosThetaT)
{
	if (cosThetaT < T(0))
		eta = 1 / eta;
	return n * (dotT<T>(i, n) * eta + cosThetaT) - i * eta;
}

///// Complete the set {a} to an orthonormal base
//extern void coordinateSystem(const Vector3f &a, Vector3f &b, Vector3f &c);
