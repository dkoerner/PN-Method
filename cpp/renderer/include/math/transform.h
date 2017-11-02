#pragma once

#include <math/common.h>
#include <math/ray.h>


/**
 * \brief Homogeneous coordinate transformation
 *
 * This class stores a general homogeneous coordinate tranformation, such as
 * rotation, translation, uniform or non-uniform scaling, and perspective
 * transformations. The inverse of this transformation is also recorded
 * here, since it is required when transforming normal vectors.
 */
template<typename T>
struct Transform
{
public:

	typedef Eigen::Matrix<T, 4, 4> M44;
	typedef TPoint<T, 3> Point;
	typedef TVector<T, 3> Vector;
	typedef TNormal3<T> Normal;
	typedef TRay<Point, Vector> Ray;

	/// Create the identity transform
	Transform() : 
		m_transform(M44::Identity()),
		m_inverse(M44::Identity()) { }

	/// Create a new transform instance for the given matrix
	Transform(const Eigen::Transform<T,3,Eigen::Affine> &trafo)
		: m_transform(trafo.matrix()), m_inverse(trafo.matrix().inverse()) { }

	/// Create a new transform instance for the given matrix 
	Transform(const M44 &trafo)
		: m_transform(trafo), m_inverse(trafo.inverse()) { }

	/// Create a new transform instance for the given matrix and its inverse
	Transform(const M44 &trafo, const M44 &inv)
		: m_transform(trafo), m_inverse(inv) { }

	/// Return the underlying matrix
	inline const M44 &getMatrix() const {
		return m_transform;
	}

	/// Return the inverse of the underlying matrix
	inline const M44 &getInverseMatrix() const {
		return m_inverse;
	}

	/// Return the inverse transformation
	Transform<T> inverse() const {
		return Transform<T>(m_inverse, m_transform);
	}

	/// Concatenate with another transform
	Transform<T> operator*(const Transform<T> &t) const;

	/// Apply the homogeneous transformation to a 3D vector
	inline Vector operator*(const Vector &v) const {
		return m_transform.template topLeftCorner<3,3>() * v;
	}

	/// Apply the homogeneous transformation to a 3D normal
	inline Normal operator*(const Normal &n) const {
		return m_inverse.template topLeftCorner<3, 3>().transpose() * n;
	}

	/// Transform a point by an arbitrary matrix in homogeneous coordinates
	inline Point operator*(const Point &p) const {
		TVector<T, 4> result = m_transform*TVector<T, 4>(p[0], p[1], p[2], 1.0);
		return result.template head<3>() / result.w();
	}

	/// Apply the homogeneous transformation to a ray
	inline Ray operator*(const Ray &r) const
	{
		return Ray(
			operator*(r.o), 
			operator*(r.d), 
			r.mint, r.maxt
		);
	}

	static Transform<T> from_aabb( const TBoundingBox<Point>& bound_dest, const TBoundingBox<Point>& bound_src = TBoundingBox<Point>(Point(0.0), Point(1.0)) )
	{
		V3d extent_src = bound_src.getExtents();
		V3d extent_dest = bound_dest.getExtents();

		Eigen::Transform<T,3,Eigen::Affine> test =
				Eigen::Translation<double, 3>(bound_dest.min-bound_src.min)*
				Eigen::DiagonalMatrix<double, 3>(extent_dest.cwiseQuotient(extent_src));

		return Transformd(test);
	}


private:
	M44 m_transform;
	M44 m_inverse;
};

typedef Transform<float> Transformf;
typedef Transform<double>  Transformd;

