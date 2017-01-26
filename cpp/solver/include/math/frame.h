#pragma once

#include <svpt/math/vector.h>


/**
 * \brief Stores a three-dimensional orthonormal coordinate frame
 *
 * This class is mostly used to quickly convert between different
 * cartesian coordinate systems and to efficiently compute certain
 * quantities (e.g. \ref cosTheta(), \ref tanTheta, ..).
 */
struct Frame
{
	Vector3f s, t;
	Normal3f n;

	/// Default constructor -- performs no initialization!
	inline Frame() { }

	/// Given a normal and tangent vectors, construct a new coordinate frame
	inline Frame(const Vector3f &s, const Vector3f &t, const Normal3f &n)
	 : s(s), t(t), n(n) { }

	/// Construct a frame from the given orthonormal vectors
	inline Frame(const Vector3f &x, const Vector3f &y, const Vector3f &z)
	 : s(x), t(y), n(z) { }

	/// Construct a new coordinate frame from a single vector
	inline Frame(const Vector3f &n) : n(n) {
		coordinateSystem(n, s, t);
	}

	/// Convert from world coordinates to local coordinates
	inline Vector3f toLocal(const Vector3f &v) const {
		return Vector3f(
			v.dot(s), v.dot(t), v.dot(n)
		);
	}

	/// Convert from local coordinates to world coordinates
	inline Vector3f toWorld(const Vector3f &v) const {
		return s * v.x() + t * v.y() + n * v.z();
	}

	/** \brief Assuming that the given direction is in the local coordinate 
	 * system, return the cosine of the angle between the normal and v */
	inline static float cosTheta(const Vector3f &v) {
		return v.z();
	}


	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the squared cosine of the angle between the normal and v */
	inline static double cosTheta2(const Vector3f &v)
	{
		return v.z() * v.z();
	}

	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the sine of the angle between the normal and v */
	inline static float sinTheta(const Vector3f &v) {
		float temp = sinTheta2(v);
		if (temp <= 0.0f)
			return 0.0f;
		return std::sqrt(temp);
	}

	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the tangent of the angle between the normal and v */
	inline static float tanTheta(const Vector3f &v) {
		float temp = 1 - v.z()*v.z();
		if (temp <= 0.0f)
			return 0.0f;
		return std::sqrt(temp) / v.z();
	}

	inline static float tanTheta2(const Vector3f &v) {
		float temp = 1 - v.z()*v.z();
		if (temp <= 0.0f)
			return 0.0f;
		return temp / (v.z() * v.z());

	}

	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the squared sine of the angle between the normal and v */
	inline static float sinTheta2(const Vector3f &v) {
		return 1.0f - v.z() * v.z();
	}

	/** \brief Assuming that the given direction is in the local coordinate 
	 * system, return the sine of the phi parameter in spherical coordinates */
	inline static float sinPhi(const Vector3f &v) {
		float sinTheta = Frame::sinTheta(v);
		if (sinTheta == 0.0f)
			return 1.0f;
		return clamp(v.y() / sinTheta, -1.0f, 1.0f);
	}

	/** \brief Assuming that the given direction is in the local coordinate 
	 * system, return the cosine of the phi parameter in spherical coordinates */
	inline static float cosPhi(const Vector3f &v) {
		float sinTheta = Frame::sinTheta(v);
		if (sinTheta == 0.0f)
			return 1.0f;
		return clamp(v.x() / sinTheta, -1.0f, 1.0f);
	}

	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the squared sine of the phi parameter in  spherical
	 * coordinates */
	inline static float sinPhi2(const Vector3f &v) {
		return clamp(v.y() * v.y() / sinTheta2(v), 0.0f, 1.0f);
	}

	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the squared cosine of the phi parameter in  spherical
	 * coordinates */
	inline static float cosPhi2(const Vector3f &v) {
		return clamp(v.x() * v.x() / sinTheta2(v), 0.0f, 1.0f);
	}

	/// \brief return reflected vector in local space
	inline static Vector3f reflect(const Vector3f &wi)
	{
		return Vector3f( -wi.x(), -wi.y(), wi.z() );
	}

	/// Equality test
	inline bool operator==(const Frame &frame) const {
		return frame.s == s && frame.t == t && frame.n == n;
	}

	/// Inequality test
	inline bool operator!=(const Frame &frame) const {
		return !operator==(frame);
	}

//	/// Return a human-readable string summary of this frame
//	inline QString toString() const {
//		return QString(
//				"Frame[\n"
//				"  s = %1,\n"
//				"  t = %2,\n"
//				"  n = %3\n"
//				"]")
//			.arg(s.toString()).arg(t.toString()).arg(n.toString());
//	}
};




struct Framed
{
	Vector3d s, t;
	Normal3d n;

	/// Default constructor -- performs no initialization!
	inline Framed() { }

	/// Given a normal and tangent vectors, construct a new coordinate frame
	inline Framed(const Vector3d &s, const Vector3d &t, const Normal3d &n)
	 : s(s), t(t), n(n) { }

	/// Construct a frame from the given orthonormal vectors
	inline Framed(const Vector3d &x, const Vector3d &y, const Vector3d &z)
	 : s(x), t(y), n(z) { }

	/// Construct a new coordinate frame from a single vector
	inline Framed(const Vector3d &n) : n(n) {
		coordinateSystem<double>(n, s, t);
	}

	/// Convert from world coordinates to local coordinates
	inline Vector3d toLocal(const Vector3d &v) const {
		return Vector3d(
			v.dot(s), v.dot(t), v.dot(n)
		);
	}

	/// Convert from local coordinates to world coordinates
	inline Vector3d toWorld(const Vector3d &v) const {
		return s * v.x() + t * v.y() + n * v.z();
	}

	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the cosine of the angle between the normal and v */
	inline static double cosTheta(const Vector3d &v) {
		return v.z();
	}

	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the squared cosine of the angle between the normal and v */
	inline static double cosTheta2(const Vector3d &v) {
		return v.z() * v.z();
	}

	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the sine of the angle between the normal and v */
	inline static double sinTheta(const Vector3d &v) {
		double temp = sinTheta2(v);
		if (temp <= 0.0)
			return 0.0;
		return std::sqrt(temp);
	}

	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the tangent of the angle between the normal and v */
	inline static double tanTheta(const Vector3d &v) {
		double temp = 1 - v.z()*v.z();
		if (temp <= 0.0)
			return 0.0;
		return std::sqrt(temp) / v.z();
	}

	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the squared tangent of the angle between the normal and v */
	inline static double tanTheta2(const Vector3d &v) {
		double temp = 1 - v.z()*v.z();
		if (temp <= 0.0f)
			return 0.0f;
		return temp / (v.z() * v.z());
	}


	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the squared sine of the angle between the normal and v */
	inline static double sinTheta2(const Vector3d &v) {
		return 1.0f - v.z() * v.z();
	}

	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the sine of the phi parameter in spherical coordinates */
	inline static double sinPhi(const Vector3d &v) {
		double sinTheta = Framed::sinTheta(v);
		if (sinTheta == 0.0)
			return 1.0;
		return clamp(v.y() / sinTheta, -1.0, 1.0);
	}

	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the cosine of the phi parameter in spherical coordinates */
	inline static double cosPhi(const Vector3d &v) {
		double sinTheta = Framed::sinTheta(v);
		if (sinTheta == 0.0)
			return 1.0;
		return clamp(v.x() / sinTheta, -1.0, 1.0);
	}

	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the squared sine of the phi parameter in  spherical
	 * coordinates */
	inline static double sinPhi2(const Vector3d &v) {
		return clamp(v.y() * v.y() / sinTheta2(v), 0.0, 1.0);
	}

	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the squared cosine of the phi parameter in  spherical
	 * coordinates */
	inline static double cosPhi2(const Vector3d &v) {
		return clamp(v.x() * v.x() / sinTheta2(v), 0.0, 1.0);
	}

	/// \brief return reflected vector in local space
	inline static Vector3d reflect(const Vector3d &wi)
	{
		return Vector3d( -wi.x(), -wi.y(), wi.z() );
	}

	/// Equality test
	inline bool operator==(const Framed &frame) const {
		return frame.s == s && frame.t == t && frame.n == n;
	}

	/// Inequality test
	inline bool operator!=(const Framed &frame) const {
		return !operator==(frame);
	}


	inline std::string toString()const
	{
		return std::string("frame: right=")+s.toString() +"  up=" + t.toString() + "  forward=" + n.toString() + "\n";
	}

//	/// Return a human-readable string summary of this frame
//	inline QString toString() const {
//		return QString(
//				"Frame[\n"
//				"  s = %1,\n"
//				"  t = %2,\n"
//				"  n = %3\n"
//				"]")
//			.arg(s.toString()).arg(t.toString()).arg(n.toString());
//	}
};
