#pragma once
#include "vector.h"



/**
 * \brief Simple n-dimensional ray segment data structure
 * 
 * Along with the ray origin and direction, this data structure additionally
 * stores a ray segment [mint, maxt] (whose entries may include positive/negative
 * infinity), as well as the componentwise reciprocals of the ray direction.
 * That is just done for convenience, as these values are frequently required.
 *
 * \remark Important: be careful when changing the ray direction. You must
 * call \ref update() to compute the componentwise reciprocals as well, or Nori's
 * ray-triangle intersection code will go haywire.
 */
template <typename _PointType, typename _VectorType> struct TRay {
	typedef _PointType                  PointType;
	typedef _VectorType                 VectorType;
	typedef typename PointType::Scalar  Scalar;

	PointType o;     ///< Ray origin
	VectorType d;    ///< Ray direction
	VectorType dRcp; ///< Componentwise reciprocals of the ray direction
	Scalar mint;     ///< Minimum position on the ray segment
	Scalar maxt;     ///< Maximum position on the ray segment

	/// Construct a new ray
	inline TRay() : mint(Epsilon), 
		maxt(std::numeric_limits<Scalar>::infinity()) { }
	
	/// Construct a new ray
	inline TRay(const PointType &o, const VectorType &d) : o(o), d(d), 
			mint(Epsilon), maxt(std::numeric_limits<Scalar>::infinity()) {
		update();
	}

	/// Construct a new ray
	inline TRay(const PointType &o, const VectorType &d, 
		Scalar mint, Scalar maxt) : o(o), d(d), mint(mint), maxt(maxt) {
		update();
	}

	/// Copy constructor
	inline TRay(const TRay &ray) 
	 : o(ray.o), d(ray.d), dRcp(ray.dRcp),
	   mint(ray.mint), maxt(ray.maxt) { }

	/// Copy a ray, but change the covered segment of the copy
	inline TRay(const TRay &ray, Scalar mint, Scalar maxt) 
	 : o(ray.o), d(ray.d), dRcp(ray.dRcp), mint(mint), maxt(maxt) { }

	/// Update the reciprocal ray directions after changing 'd'
	inline void update() {
		dRcp = d.cwiseInverse();
	}

	/// Return the position of a point along the ray
	inline PointType operator() (Scalar t) const { return o + t * d; }

	/// Return a ray that points into the opposite direction
	inline Ray3f reverse() const {
		Ray3f result;
		result.o = o; result.d = -d; result.dRcp = -dRcp;
		result.mint = mint; result.maxt = maxt;
		return result;
	}

	inline std::string toString() const
	{
		std::string result;
		result = o.toString() + " " + d.toString() + " " + ::toString(mint)+ " " + ::toString(maxt);
		return result;
	}

//	/// Return a human-readable string summary of this ray
//	inline QString toString() const {
//		return QString(
//				"Ray[\n"
//				"  o = %1,\n"
//				"  d = %2,\n"
//				"  mint = %3,\n"
//				"  maxt = %4\n"
//				"]")
//			.arg(o.toString())
//			.arg(d.toString())
//			.arg(mint)
//			.arg(maxt);
//	}
};


template<typename T>
inline void closestPoint(const TRay<TPoint<T, 3>, TVector<T, 3>> &r1,
						 const TRay<TPoint<T, 3>, TVector<T, 3>> &r2,
						 T &t1,
						 T &t2)
{
	TVector<T, 3> p13 = r1.o - r2.o;
	TVector<T, 3> p43 = r2.d;
	TVector<T, 3> p21 = r1.d;

	T d1343 = p13.x() * p43.x() + p13.y() * p43.y() + p13.z() * p43.z();
	T d4321 = p43.x() * p21.x() + p43.y() * p21.y() + p43.z() * p21.z();
	T d1321 = p13.x() * p21.x() + p13.y() * p21.y() + p13.z() * p21.z();
	T d4343 = p43.x() * p43.x() + p43.y() * p43.y() + p43.z() * p43.z();
	T d2121 = p21.x() * p21.x() + p21.y() * p21.y() + p21.z() * p21.z();

	T denom = d2121 * d4343 - d4321 * d4321;
	T numer = d1343 * d4321 - d1321 * d4343;

	t1 = numer / denom;
	t2 = (d1343 + d4321 * t1) / d4343;
}



template<typename T>
inline TPoint<T, 3> closestPoint(const TRay<TPoint<T, 3>, TVector<T, 3>> &r,
								 const TPoint<T, 3>& p )
{
	return r(r.d.dot( p-r.o ));
}
