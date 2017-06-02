#pragma once
#include <math/common.h>

/**
 * \brief Represents a linear RGB color value
 */
struct Color3f : public Eigen::Array3f {
public:
	typedef Eigen::Array3f Base;

	/// Initialize the color vector with a uniform value
	inline Color3f(float value = 0) : Base(value, value, value) { }

	/// Initialize the color vector with specific per-channel values
	inline Color3f(float r, float g, float b) : Base(r, g, b) { }

	/// Construct a color vector from ArrayBase (needed to play nice with Eigen)
	template <typename Derived> inline Color3f(const Eigen::ArrayBase<Derived>& p) 
		: Base(p) { }

	/// Assign a color vector from ArrayBase (needed to play nice with Eigen)
    template <typename Derived> Color3f &operator=(const Eigen::ArrayBase<Derived>& p) {
		this->Base::operator=(p);
		return *this;
    }

	/// Return a reference to the red channel
	inline float &r() { return x(); }
	/// Return a reference to the red channel (const version)
	inline const float &r() const { return x(); }
	/// Return a reference to the green channel
	inline float &g() { return y(); }
	/// Return a reference to the green channel (const version)
	inline const float &g() const { return y(); }
	/// Return a reference to the blue channel
	inline float &b() { return z(); }
	/// Return a reference to the blue channel (const version)
	inline const float &b() const { return z(); }

	/// Clamp to the positive range
	inline Color3f clamp() const { return Color3f(std::max(r(), 0.0f),
		std::max(g(), 0.0f), std::max(b(), 0.0f)); }

	/// Check if the color vector contains a NaN/Inf/negative value
	bool isValid() const {
		for (int i=0; i<3; ++i) {
			float value = coeff(i);
			int cl = std::fpclassify(value);
			if (value < 0 || cl == FP_INFINITE || cl == FP_NAN)
				return false;
		}
		return true;
	}

	/// Convert from sRGB to linear RGB
	Color3f toLinearRGB() const {
		Color3f result;

		for (int i=0; i<3; ++i) {
			float value = coeff(i);

			if (value <= 0.04045f)
				result[i] = value * (1.0f / 12.92f);
			else
				result[i] = std::pow((value + 0.055f)
					* (1.0f / 1.055f), 2.4f);
		}

		return result;
	}

	/// Convert from linear RGB to sRGB
	Color3f toSRGB() const {
		Color3f result;

		for (int i=0; i<3; ++i) {
			float value = coeff(i);

			if (value <= 0.0031308f)
				result[i] = 12.92f * value;
			else
				result[i] = (1.0f + 0.055f)
					* std::pow(value, 1.0f/2.4f) -  0.055f;
		}

		return result;
	}

	/// Return the associated luminance
	float getLuminance() const
	{
		return x() * 0.212671f + y() * 0.715160f + z() * 0.072169f;
	}

//	/// Return a human-readable string summary
//	inline QString toString() const {
//		return QString("[%1, %2, %3]").arg(coeff(0)).arg(coeff(1)).arg(coeff(2));
//	}

	std::string toString()const
	{
		std::string result = "(";
		for (size_t i=0; i<3; ++i)
		{
			result += ::toString(this->coeff(i));
			if (i+1 < 3)
				result += ", ";
		}
		return result + ")";
	}
};

/**
 * \brief Represents a linear RGB color and a weight
 *
 * This is used by Nori's image reconstruction filter code
 */
/*
struct Color4f : public Eigen::Array4f {
public:
	typedef Eigen::Array4f Base;

	/// Create an zero value
	inline Color4f() : Base(0.0f, 0.0f, 0.0f, 0.0f) { }

	/// Create from a 3-channel color
	inline Color4f(const Color3f &c) : Base(c.r(), c.g(), c.b(), 1.0f) { }

	/// Initialize the color vector with specific per-channel values
	inline Color4f(float r, float g, float b, float w) : Base(r, g, b, w) { }

	/// Construct a color vector from ArrayBase (needed to play nice with Eigen)
	template <typename Derived> inline Color4f(const Eigen::ArrayBase<Derived>& p) 
		: Base(p) { }

	/// Assign a color vector from ArrayBase (needed to play nice with Eigen)
    template <typename Derived> Color4f &operator=(const Eigen::ArrayBase<Derived>& p) {
		this->Base::operator=(p);
		return *this;
    }

	/// Normalize and convert into a \ref Color3f value 
	inline Color3f normalized() const {
		if (EXPECT_TAKEN(w() != 0))
			return head<3>() / w();
		else
			return Color3f(0.0f);
	}

	/// Return a human-readable string summary
	inline QString toString() const {
		return QString("[%1, %2, %3, %4]").arg(coeff(0)).arg(coeff(1)).arg(coeff(2)).arg(coeff(4));
	}
};
*/

