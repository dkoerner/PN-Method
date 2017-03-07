#include <util/moexp.h>


#include <map>
#include <util/cas.h>
#include <houio/Geometry.h>
#include <houio/HouGeoIO.h>



namespace sh {

namespace {

// Number of precomputed factorials and double-factorials that can be
// returned in constant time.
const int kCacheSize = 16;

const int kHardCodedOrderLimit = 4;

const int kIrradianceOrder = 2;
const int kIrradianceCoeffCount = GetCoefficientCount(kIrradianceOrder);

// For efficiency, the cosine lobe for normal = (0, 0, 1) as the first 9
// spherical harmonic coefficients are hardcoded below. This was computed by
// evaluating:
//   ProjectFunction(kIrradianceOrder, [] (double phi, double theta) {
//     return Clamp(Eigen::Vector3d::UnitZ().dot(ToVector(phi, theta)),
//                  0.0, 1.0);
//   }, 10000000);
const std::vector<double> cosine_lobe = { 0.886227, 0.0, 1.02333, 0.0, 0.0, 0.0,
										  0.495416, 0.0, 0.0 };

// A zero template is required for EvalSHSum to handle its template
// instantiations and a type's default constructor does not necessarily
// initialize to zero.
template<typename T> T Zero();
template<> double Zero() { return 0.0; }
template<> float Zero() { return 0.0; }
template<> Eigen::Array3f Zero() { return Eigen::Array3f::Zero(); }

template <class T>
using VectorX = Eigen::Matrix<T, Eigen::Dynamic, 1>;

// Usage: CHECK(bool, string message);
// Note that it must end a semi-colon, making it look like a
// valid C++ statement (hence the awkward do() while(false)).
#ifndef NDEBUG
# define CHECK(condition, message) \
  do { \
	if (!(condition)) { \
	  std::cerr << "Check failed (" #condition ") in " << __FILE__ \
		<< ":" << __LINE__ << ", message: " << message << std::endl; \
	  std::exit(EXIT_FAILURE); \
	} \
  } while(false)
#else
# define ASSERT(condition, message) do {} while(false)
#endif

// Clamp the first argument to be greater than or equal to the second
// and less than or equal to the third.
double Clamp(double val, double min, double max) {
  if (val < min) {
	val = min;
  }
  if (val > max) {
	val = max;
  }
  return val;
}

// Return true if the first value is within epsilon of the second value.
bool NearByMargin(double actual, double expected) {
  double diff = actual - expected;
  if (diff < 0.0) {
	diff = -diff;
  }
  // 5 bits of error in mantissa (source of '32 *')
  return diff < 32 * std::numeric_limits<double>::epsilon();
}

// Return floating mod x % m.
double FastFMod(double x, double m) {
  return x - (m * floor(x / m));
}

// Hardcoded spherical harmonic functions for low orders (l is first number
// and m is second number (sign encoded as preceeding 'p' or 'n')).
//
// As polynomials they are evaluated more efficiently in cartesian coordinates,
// assuming that @d is unit. This is not verified for efficiency.
double HardcodedSH00(const Eigen::Vector3d& d) {
  // 0.5 * sqrt(1/pi)
  return 0.282095;
}

double HardcodedSH1n1(const Eigen::Vector3d& d) {
  // -sqrt(3/(4pi)) * y
  return -0.488603 * d.y();
}

double HardcodedSH10(const Eigen::Vector3d& d) {
  // sqrt(3/(4pi)) * z
  return 0.488603 * d.z();
}

double HardcodedSH1p1(const Eigen::Vector3d& d) {
  // -sqrt(3/(4pi)) * x
  return -0.488603 * d.x();
}

double HardcodedSH2n2(const Eigen::Vector3d& d) {
  // 0.5 * sqrt(15/pi) * x * y
  return 1.092548 * d.x() * d.y();
}

double HardcodedSH2n1(const Eigen::Vector3d& d) {
  // -0.5 * sqrt(15/pi) * y * z
  return -1.092548 * d.y() * d.z();
}

double HardcodedSH20(const Eigen::Vector3d& d) {
  // 0.25 * sqrt(5/pi) * (-x^2-y^2+2z^2)
  return 0.315392 * (-d.x() * d.x() - d.y() * d.y() + 2.0 * d.z() * d.z());
}

double HardcodedSH2p1(const Eigen::Vector3d& d) {
  // -0.5 * sqrt(15/pi) * x * z
  return -1.092548 * d.x() * d.z();
}

double HardcodedSH2p2(const Eigen::Vector3d& d) {
  // 0.25 * sqrt(15/pi) * (x^2 - y^2)
  return 0.546274 * (d.x() * d.x() - d.y() * d.y());
}

double HardcodedSH3n3(const Eigen::Vector3d& d) {
  // -0.25 * sqrt(35/(2pi)) * y * (3x^2 - y^2)
  return -0.590044 * d.y() * (3.0 * d.x() * d.x() - d.y() * d.y());
}

double HardcodedSH3n2(const Eigen::Vector3d& d) {
  // 0.5 * sqrt(105/pi) * x * y * z
  return 2.890611 * d.x() * d.y() * d.z();
}

double HardcodedSH3n1(const Eigen::Vector3d& d) {
  // -0.25 * sqrt(21/(2pi)) * y * (4z^2-x^2-y^2)
  return -0.457046 * d.y() * (4.0 * d.z() * d.z() - d.x() * d.x()
							 - d.y() * d.y());
}

double HardcodedSH30(const Eigen::Vector3d& d) {
  // 0.25 * sqrt(7/pi) * z * (2z^2 - 3x^2 - 3y^2)
  return 0.373176 * d.z() * (2.0 * d.z() * d.z() - 3.0 * d.x() * d.x()
							 - 3.0 * d.y() * d.y());
}

double HardcodedSH3p1(const Eigen::Vector3d& d) {
  // -0.25 * sqrt(21/(2pi)) * x * (4z^2-x^2-y^2)
  return -0.457046 * d.x() * (4.0 * d.z() * d.z() - d.x() * d.x()
							 - d.y() * d.y());
}

double HardcodedSH3p2(const Eigen::Vector3d& d) {
  // 0.25 * sqrt(105/pi) * z * (x^2 - y^2)
  return 1.445306 * d.z() * (d.x() * d.x() - d.y() * d.y());
}

double HardcodedSH3p3(const Eigen::Vector3d& d) {
  // -0.25 * sqrt(35/(2pi)) * x * (x^2-3y^2)
  return -0.590044 * d.x() * (d.x() * d.x() - 3.0 * d.y() * d.y());
}

double HardcodedSH4n4(const Eigen::Vector3d& d) {
  // 0.75 * sqrt(35/pi) * x * y * (x^2-y^2)
  return 2.503343 * d.x() * d.y() * (d.x() * d.x() - d.y() * d.y());
}

double HardcodedSH4n3(const Eigen::Vector3d& d) {
  // -0.75 * sqrt(35/(2pi)) * y * z * (3x^2-y^2)
  return -1.770131 * d.y() * d.z() * (3.0 * d.x() * d.x() - d.y() * d.y());
}

double HardcodedSH4n2(const Eigen::Vector3d& d) {
  // 0.75 * sqrt(5/pi) * x * y * (7z^2-1)
  return 0.946175 * d.x() * d.y() * (7.0 * d.z() * d.z() - 1.0);
}

double HardcodedSH4n1(const Eigen::Vector3d& d) {
  // -0.75 * sqrt(5/(2pi)) * y * z * (7z^2-3)
  return -0.669047 * d.y() * d.z() * (7.0 * d.z() * d.z() - 3.0);
}

double HardcodedSH40(const Eigen::Vector3d& d) {
  // 3/16 * sqrt(1/pi) * (35z^4-30z^2+3)
  double z2 = d.z() * d.z();
  return 0.105786 * (35.0 * z2 * z2 - 30.0 * z2 + 3.0);
}

double HardcodedSH4p1(const Eigen::Vector3d& d) {
  // -0.75 * sqrt(5/(2pi)) * x * z * (7z^2-3)
  return -0.669047 * d.x() * d.z() * (7.0 * d.z() * d.z() - 3.0);
}

double HardcodedSH4p2(const Eigen::Vector3d& d) {
  // 3/8 * sqrt(5/pi) * (x^2 - y^2) * (7z^2 - 1)
  return 0.473087 * (d.x() * d.x() - d.y() * d.y())
	  * (7.0 * d.z() * d.z() - 1.0);
}

double HardcodedSH4p3(const Eigen::Vector3d& d) {
  // -0.75 * sqrt(35/(2pi)) * x * z * (x^2 - 3y^2)
  return -1.770131 * d.x() * d.z() * (d.x() * d.x() - 3.0 * d.y() * d.y());
}

double HardcodedSH4p4(const Eigen::Vector3d& d) {
  // 3/16*sqrt(35/pi) * (x^2 * (x^2 - 3y^2) - y^2 * (3x^2 - y^2))
  double x2 = d.x() * d.x();
  double y2 = d.y() * d.y();
  return 0.625836 * (x2 * (x2 - 3.0 * y2) - y2 * (3.0 * x2 - y2));
}

// Compute the factorial for an integer @x. It is assumed x is at least 0.
// This implementation precomputes the results for low values of x, in which
// case this is a constant time lookup.
//
// The vast majority of SH evaluations will hit these precomputed values.
double Factorial(int x) {
  const double factorial_cache[kCacheSize] = {1, 1, 2, 6, 24, 120, 720, 5040,
											  40320, 362880, 3628800, 39916800,
											  479001600, 6227020800,
											  87178291200, 1307674368000};

  if (x < kCacheSize) {
	return factorial_cache[x];
  } else {
	double s = 1.0;
	for (int n = 2; n <= x; n++) {
	  s *= n;
	}
	return s;
  }
}

// Compute the double factorial for an integer @x. This assumes x is at least
// 0.  This implementation precomputes the results for low values of x, in
// which case this is a constant time lookup.
//
// The vast majority of SH evaluations will hit these precomputed values.
// See http://mathworld.wolfram.com/DoubleFactorial.html
double DoubleFactorial(int x) {
  const double dbl_factorial_cache[kCacheSize] = {1, 1, 2, 3, 8, 15, 48, 105,
												  384, 945, 3840, 10395, 46080,
												  135135, 645120, 2027025};

  if (x < kCacheSize) {
	return dbl_factorial_cache[x];
  } else {
	double s = 1.0;
	double n = x;
	while (n > 1.0) {
	  s *= n;
	  n -= 2.0;
	}
	return s;
  }
}

// Evaluate the associated Legendre polynomial of degree @l and order @m at
// coordinate @x. The inputs must satisfy:
// 1. l >= 0
// 2. 0 <= m <= l
// 3. -1 <= x <= 1
// See http://en.wikipedia.org/wiki/Associated_Legendre_polynomials
//
// This implementation is based off the approach described in [1],
// instead of computing Pml(x) directly, Pmm(x) is computed. Pmm can be
// lifted to Pmm+1 recursively until Pml is found
double EvalLegendrePolynomial(int l, int m, double x) {
  // Compute Pmm(x) = (-1)^m(2m - 1)!!(1 - x^2)^(m/2), where !! is the double
  // factorial.
  double pmm = 1.0;
  // P00 is defined as 1.0, do don't evaluate Pmm unless we know m > 0
  if (m > 0) {
	double sign = (m % 2 == 0 ? 1 : -1);
	pmm = sign * DoubleFactorial(2 * m - 1) * pow(1 - x * x, m / 2.0);
  }

  if (l == m) {
	// Pml is the same as Pmm so there's no lifting to higher bands needed
	return pmm;
  }

  // Compute Pmm+1(x) = x(2m + 1)Pmm(x)
  double pmm1 = x * (2 * m + 1) * pmm;
  if (l == m + 1) {
	// Pml is the same as Pmm+1 so we are done as well
	return pmm1;
  }

  // Use the last two computed bands to lift up to the next band until l is
  // reached, using the recurrence relationship:
  // Pml(x) = (x(2l - 1)Pml-1 - (l + m - 1)Pml-2) / (l - m)
  for (int n = m + 2; n <= l; n++) {
	double pmn = (x * (2 * n - 1) * pmm1 - (n + m - 1) * pmm) / (n - m);
	pmm = pmm1;
	pmm1 = pmn;
  }
  // Pmm1 at the end of the above loop is equal to Pml
  return pmm1;
}

// ---- The following functions are used to implement SH rotation computations
//      based on the recursive approach described in [1, 4]. The names of the
//      functions correspond with the notation used in [1, 4].

// See http://en.wikipedia.org/wiki/Kronecker_delta
double KroneckerDelta(int i, int j) {
  if (i == j) {
	return 1.0;
  } else {
	return 0.0;
  }
}

// [4] uses an odd convention of referring to the rows and columns using
// centered indices, so the middle row and column are (0, 0) and the upper
// left would have negative coordinates.
//
// This is a convenience function to allow us to access an Eigen::MatrixXd
// in the same manner, assuming r is a (2l+1)x(2l+1) matrix.
double GetCenteredElement(const Eigen::MatrixXd& r, int i, int j) {
  // The shift to go from [-l, l] to [0, 2l] is (rows - 1) / 2 = l,
  // (since the matrix is assumed to be square, rows == cols).
  int offset = (r.rows() - 1) / 2;
  return r(i + offset, j + offset);
}

// P is a helper function defined in [4] that is used by the functions U, V, W.
// This should not be called on its own, as U, V, and W (and their coefficients)
// select the appropriate matrix elements to access (arguments @a and @b).
double P(int i, int a, int b, int l, const std::vector<Eigen::MatrixXd>& r) {
  if (b == l) {
	return GetCenteredElement(r[1], i, 1) *
		GetCenteredElement(r[l - 1], a, l - 1) -
		GetCenteredElement(r[1], i, -1) *
		GetCenteredElement(r[l - 1], a, -l + 1);
  } else if (b == -l) {
	return GetCenteredElement(r[1], i, 1) *
		GetCenteredElement(r[l - 1], a, -l + 1) +
		GetCenteredElement(r[1], i, -1) *
		GetCenteredElement(r[l - 1], a, l - 1);
  } else {
	return GetCenteredElement(r[1], i, 0) * GetCenteredElement(r[l - 1], a, b);
  }
}

// The functions U, V, and W should only be called if the correspondingly
// named coefficient u, v, w from the function ComputeUVWCoeff() is non-zero.
// When the coefficient is 0, these would attempt to access matrix elements that
// are out of bounds. The list of rotations, @r, must have the @l - 1
// previously completed band rotations. These functions are valid for l >= 2.

double U(int m, int n, int l, const std::vector<Eigen::MatrixXd>& r) {
  // Although [1, 4] split U into three cases for m == 0, m < 0, m > 0
  // the actual values are the same for all three cases
  return P(0, m, n, l, r);
}

double V(int m, int n, int l, const std::vector<Eigen::MatrixXd>& r) {
  if (m == 0) {
	return P(1, 1, n, l, r) + P(-1, -1, n, l, r);
  } else if (m > 0) {
	return P(1, m - 1, n, l, r) * sqrt(1 + KroneckerDelta(m, 1)) -
		P(-1, -m + 1, n, l, r) * (1 - KroneckerDelta(m, 1));
  } else {
	// Note there is apparent errata in [1,4,4b] dealing with this particular
	// case. [4b] writes it should be P*(1-d)+P*(1-d)^0.5
	// [1] writes it as P*(1+d)+P*(1-d)^0.5, but going through the math by hand,
	// you must have it as P*(1-d)+P*(1+d)^0.5 to form a 2^.5 term, which
	// parallels the case where m > 0.
	return P(1, m + 1, n, l, r) * (1 - KroneckerDelta(m, -1)) +
		P(-1, -m - 1, n, l, r) * sqrt(1 + KroneckerDelta(m, -1));
  }
}

double W(int m, int n, int l, const std::vector<Eigen::MatrixXd>& r) {
  if (m == 0) {
	// whenever this happens, w is also 0 so W can be anything
	return 0.0;
  } else if (m > 0) {
	return P(1, m + 1, n, l, r) + P(-1, -m - 1, n, l, r);
  } else {
	return P(1, m - 1, n, l, r) - P(-1, -m + 1, n, l, r);
  }
}

// Calculate the coefficients applied to the U, V, and W functions. Because
// their equations share many common terms they are computed simultaneously.
void ComputeUVWCoeff(int m, int n, int l, double* u, double* v, double* w) {
  double d = KroneckerDelta(m, 0);
  double denom = (abs(n) == l ? 2.0 * l * (2.0 * l - 1) : (l + n) * (l - n));

  *u = sqrt((l + m) * (l - m) / denom);
  *v = 0.5 * sqrt((1 + d) * (l + abs(m) - 1.0) * (l + abs(m)) / denom)
	  * (1 - 2 * d);
  *w = -0.5 * sqrt((l - abs(m) - 1) * (l - abs(m)) / denom) * (1 - d);
}

// Calculate the (2l+1)x(2l+1) rotation matrix for the band @l.
// This uses the matrices computed for band 1 and band l-1 to compute the
// matrix for band l. @rotations must contain the previously computed l-1
// rotation matrices, and the new matrix for band l will be appended to it.
//
// This implementation comes from p. 5 (6346), Table 1 and 2 in [4] taking
// into account the corrections from [4b].
void ComputeBandRotation(int l, std::vector<Eigen::MatrixXd>* rotations) {
  // The band's rotation matrix has rows and columns equal to the number of
  // coefficients within that band (-l <= m <= l implies 2l + 1 coefficients).
  Eigen::MatrixXd rotation(2 * l + 1, 2 * l + 1);
  for (int m = -l; m <= l; m++) {
	for (int n = -l; n <= l; n++) {
	  double u, v, w;
	  ComputeUVWCoeff(m, n, l, &u, &v, &w);

	  // The functions U, V, W are only safe to call if the coefficients
	  // u, v, w are not zero
	  if (!NearByMargin(u, 0.0))
		  u *= U(m, n, l, *rotations);
	  if (!NearByMargin(v, 0.0))
		  v *= V(m, n, l, *rotations);
	  if (!NearByMargin(w, 0.0))
		  w *= W(m, n, l, *rotations);

	  rotation(m + l, n + l) = (u + v + w);
	}
  }

  rotations->push_back(rotation);
}

}  // namespace

Eigen::Vector3d ToVector(double phi, double theta) {
  double r = sin(theta);
  return Eigen::Vector3d(r * cos(phi), r * sin(phi), cos(theta));
}

void ToSphericalCoords(const Eigen::Vector3d& dir, double* phi, double* theta) {
  //CHECK(NearByMargin(dir.squaredNorm(), 1.0), "dir is not unit");
  // Explicitly clamp the z coordinate so that numeric errors don't cause it
  // to fall just outside of acos' domain.
  *theta = acos(Clamp(dir.z(), -1.0, 1.0));
  // We don't need to divide dir.y() or dir.x() by sin(theta) since they are
  // both scaled by it and atan2 will handle it appropriately.
  *phi = atan2(dir.y(), dir.x());
}

double ImageXToPhi(int x, int width) {
  // The directions are measured from the center of the pixel, so add 0.5
  // to convert from integer pixel indices to float pixel coordinates.
  return 2.0 * M_PI * (x + 0.5) / width;
}

double ImageYToTheta(int y, int height) {
  return M_PI * (y + 0.5) / height;
}

Eigen::Vector2d ToImageCoords(double phi, double theta, int width, int height) {
  // Allow theta to repeat and map to 0 to pi. However, to account for cases
  // where y goes beyond the normal 0 to pi range, phi may need to be adjusted.
  theta = Clamp(FastFMod(theta, 2.0 * M_PI), 0.0, 2.0 * M_PI);
  if (theta > M_PI) {
	// theta is out of bounds. Effectively, theta has rotated past the pole
	// so after adjusting theta to be in range, rotating phi by pi forms an
	// equivalent direction.
	theta = 2.0 * M_PI - theta;  // now theta is between 0 and pi
	phi += M_PI;
  }
  // Allow phi to repeat and map to the normal 0 to 2pi range.
  // Clamp and map after adjusting theta in case theta was forced to update phi.
  phi = Clamp(FastFMod(phi, 2.0 * M_PI), 0.0, 2.0 * M_PI);

  // Now phi is in [0, 2pi] and theta is in [0, pi] so it's simple to inverse
  // the linear equations in ImageCoordsToSphericalCoords, although there's no
  // -0.5 because we're returning floating point coordinates and so don't need
  // to center the pixel.
  return Eigen::Vector2d(width * phi / (2.0 * M_PI), height * theta / M_PI);
}

double EvalSHSlow(int l, int m, double phi, double theta) {
  //CHECK(l >= 0, "l must be at least 0.");
  //CHECK(-l <= m && m <= l, "m must be between -l and l.");

  double kml = sqrt((2.0 * l + 1) * Factorial(l - abs(m)) /
					(4.0 * M_PI * Factorial(l + abs(m))));
  if (m > 0) {
	return sqrt(2.0) * kml * cos(m * phi) *
		EvalLegendrePolynomial(l, m, cos(theta));
  } else if (m < 0) {
	return sqrt(2.0) * kml * sin(-m * phi) *
		EvalLegendrePolynomial(l, -m, cos(theta));
  } else {
	return kml * EvalLegendrePolynomial(l, 0, cos(theta));
  }
}

double EvalSHSlow(int l, int m, const Eigen::Vector3d& dir) {
  double phi, theta;
  ToSphericalCoords(dir, &phi, &theta);
  return EvalSH(l, m, phi, theta);
}

double EvalSH(int l, int m, double phi, double theta) {
  // If using the hardcoded functions, switch to cartesian
  if (l <= kHardCodedOrderLimit) {
	return EvalSH(l, m, ToVector(phi, theta));
  } else {
	// Stay in spherical coordinates since that's what the recurrence
	// version is implemented in
	return EvalSHSlow(l, m, phi, theta);
  }
}

double EvalSH(int l, int m, const Eigen::Vector3d& dir) {
  if (l <= kHardCodedOrderLimit) {
	// Validate l and m here (don't do it generally since EvalSHSlow also
	// checks it if we delegate to that function).
	//CHECK(l >= 0, "l must be at least 0.");
	//CHECK(-l <= m && m <= l, "m must be between -l and l.");
	//CHECK(NearByMargin(dir.squaredNorm(), 1.0), "dir is not unit.");

	switch (l) {
	  case 0:
		return HardcodedSH00(dir);
	  case 1:
		switch (m) {
		  case -1:
			return HardcodedSH1n1(dir);
		  case 0:
			return HardcodedSH10(dir);
		  case 1:
			return HardcodedSH1p1(dir);
		}
	  case 2:
		switch (m) {
		  case -2:
			return HardcodedSH2n2(dir);
		  case -1:
			return HardcodedSH2n1(dir);
		  case 0:
			return HardcodedSH20(dir);
		  case 1:
			return HardcodedSH2p1(dir);
		  case 2:
			return HardcodedSH2p2(dir);
		}
	  case 3:
		switch (m) {
		  case -3:
			return HardcodedSH3n3(dir);
		  case -2:
			return HardcodedSH3n2(dir);
		  case -1:
			return HardcodedSH3n1(dir);
		  case 0:
			return HardcodedSH30(dir);
		  case 1:
			return HardcodedSH3p1(dir);
		  case 2:
			return HardcodedSH3p2(dir);
		  case 3:
			return HardcodedSH3p3(dir);
		}
	  case 4:
		switch (m) {
		  case -4:
			return HardcodedSH4n4(dir);
		  case -3:
			return HardcodedSH4n3(dir);
		  case -2:
			return HardcodedSH4n2(dir);
		  case -1:
			return HardcodedSH4n1(dir);
		  case 0:
			return HardcodedSH40(dir);
		  case 1:
			return HardcodedSH4p1(dir);
		  case 2:
			return HardcodedSH4p2(dir);
		  case 3:
			return HardcodedSH4p3(dir);
		  case 4:
			return HardcodedSH4p4(dir);
		}
	}

	// This is unreachable given the CHECK's above but the compiler can't tell.
	return 0.0;
  } else {
	// Not hard-coded so use the recurrence relation (which will convert this
	// to spherical coordinates).
	return EvalSHSlow(l, m, dir);
  }
}

std::unique_ptr<std::vector<double>> ProjectFunction(
	int order, const SphericalFunction& func, int sample_count) {
  //CHECK(order >= 0, "Order must be at least zero.");
  //CHECK(sample_count > 0, "Sample count must be at least one.");

  // This is the approach demonstrated in [1] and is useful for arbitrary
  // functions on the sphere that are represented analytically.
  const int sample_side = static_cast<int>(floor(sqrt(sample_count)));
  std::unique_ptr<std::vector<double>> coeffs(new std::vector<double>());
  coeffs->assign(GetCoefficientCount(order), 0.0);

  // generate sample_side^2 uniformly and stratified samples over the sphere
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> rng(0.0, 1.0);
  for (int t = 0; t < sample_side; t++) {
	for (int p = 0; p < sample_side; p++) {
	  double alpha = (t + rng(gen)) / sample_side;
	  double beta = (p + rng(gen)) / sample_side;
	  // See http://www.bogotobogo.com/Algorithms/uniform_distribution_sphere.php
	  double phi = 2.0 * M_PI * beta;
	  double theta = acos(2.0 * alpha - 1.0);

	  // evaluate the analytic function for the current spherical coords
	  double func_value = func(phi, theta);

	  // evaluate the SH basis functions up to band O, scale them by the
	  // function's value and accumulate them over all generated samples
	  for (int l = 0; l <= order; l++) {
		for (int m = -l; m <= l; m++) {
		  double sh = EvalSH(l, m, phi, theta);
		  (*coeffs)[GetIndex(l, m)] += func_value * sh;
		}
	  }
	}
  }

  // scale by the probability of a particular sample, which is
  // 4pi/sample_side^2. 4pi for the surface area of a unit sphere, and
  // 1/sample_side^2 for the number of samples drawn uniformly.
  double weight = 4.0 * M_PI / (sample_side * sample_side);
  for (unsigned int i = 0; i < coeffs->size(); i++) {
	 (*coeffs)[i] *= weight;
  }

  return coeffs;
}


template <typename T>
T EvalSHSum(int order, const std::vector<T>& coeffs, double phi, double theta) {
  if (order <= kHardCodedOrderLimit) {
	// It is faster to compute the cartesian coordinates once
	return EvalSHSum(order, coeffs, ToVector(phi, theta));
  }

  //CHECK(GetCoefficientCount(order) == coeffs.size(),
  //  "Incorrect number of coefficients provided.");
  T sum = Zero<T>();
  for (int l = 0; l <= order; l++) {
	for (int m = -l; m <= l; m++) {
	  sum += EvalSH(l, m, phi, theta) * coeffs[GetIndex(l, m)];
	}
  }
  return sum;
}

template <typename T>
T EvalSHSum(int order, const std::vector<T>& coeffs,
			const Eigen::Vector3d& dir) {
  if (order > kHardCodedOrderLimit) {
	// It is faster to switch to spherical coordinates
	double phi, theta;
	ToSphericalCoords(dir, &phi, &theta);
	return EvalSHSum(order, coeffs, phi, theta);
  }

  //CHECK(GetCoefficientCount(order) == coeffs.size(), "Incorrect number of coefficients provided.");
  //CHECK(NearByMargin(dir.squaredNorm(), 1.0), "dir is not unit.");

  T sum = Zero<T>();
  for (int l = 0; l <= order; l++) {
	for (int m = -l; m <= l; m++) {
	  sum += EvalSH(l, m, dir) * coeffs[GetIndex(l, m)];
	}
  }
  return sum;
}


// ---- Template specializations -----------------------------------------------

template double EvalSHSum<double>(int order, const std::vector<double>& coeffs,
								  double phi, double theta);
template double EvalSHSum<double>(int order, const std::vector<double>& coeffs,
								  const Eigen::Vector3d& dir);

template float EvalSHSum<float>(int order, const std::vector<float>& coeffs,
								double phi, double theta);
template float EvalSHSum<float>(int order, const std::vector<float>& coeffs,
								const Eigen::Vector3d& dir);

template Eigen::Array3f EvalSHSum<Eigen::Array3f>(
	int order,  const std::vector<Eigen::Array3f>& coeffs,
	double phi, double theta);
template Eigen::Array3f EvalSHSum<Eigen::Array3f>(
	int order,  const std::vector<Eigen::Array3f>& coeffs,
	const Eigen::Vector3d& dir);

}  // namespace sh




namespace moexp
{

	int ipow(int base, int exp)
	{
		int result = 1;
		while (exp)
		{
			if (exp & 1)
				result *= base;
			exp >>= 1;
			base *= base;
		}

		return result;
	}

	//const double M_PI = 3.14159265358979323846;

	// Number of precomputed factorials and double-factorials that can be
	// returned in constant time.
	const int kCacheSize = 16;

	// Compute the factorial for an integer @x. It is assumed x is at least 0.
	// This implementation precomputes the results for low values of x, in which
	// case this is a constant time lookup.
	//
	// The vast majority of SH evaluations will hit these precomputed values.
	double factorial(int x)
	{
		const double factorial_cache[kCacheSize] = {1, 1, 2, 6, 24, 120, 720, 5040,
												  40320, 362880, 3628800, 39916800,
												  479001600, 6227020800,
												  87178291200, 1307674368000};

		if (x < kCacheSize)
		{
			return factorial_cache[x];
		}else
		{
			double s = 1.0;
			for (int n = 2; n <= x; n++)
			{
				s *= n;
			}
			return s;
		}
	}

	// Compute the double factorial for an integer @x. This assumes x is at least
	// 0.  This implementation precomputes the results for low values of x, in
	// which case this is a constant time lookup.
	//
	// The vast majority of SH evaluations will hit these precomputed values.
	// See http://mathworld.wolfram.com/DoubleFactorial.html
	double doubleFactorial(int x)
	{
	  const double dbl_factorial_cache[kCacheSize] = {1, 1, 2, 3, 8, 15, 48, 105,
													  384, 945, 3840, 10395, 46080,
													  135135, 645120, 2027025};

		if (x < kCacheSize)
		{
			return dbl_factorial_cache[x];
		}else
		{
			double s = 1.0;
			double n = x;
			while (n > 1.0)
			{
				s *= n;
				n -= 2.0;
			}
			return s;
		}
	}

	int numTensorComponents( int order, int dimension )
	{
		return ipow(dimension, order);
	}

	// condon-shortley phase
	double csp( int m )
	{
		return (m % 2 == 0 ? 1.0 : -1.0);
	}

	double C(int l, int m)
	{
		double a = csp(m);
		double b1 = (2*l+1)*INV_FOURPI;
		double b2 = factorial(l-m)/factorial(l+m);
		return a*std::sqrt(b1*b2);
		return a;
	}

	double a( int l, int m, int j )
	{
		double a = csp(j);
		double b = ipow(2, l)*factorial(j)*factorial(l-j);
		double c = factorial(2*l-2*j);
		double frac1 = a/b;
		double frac2;
		double d_fac = l-m-2*j;
		// it appears that fractions which contain negative factorials are considered zero by convention
		// see http://mathoverflow.net/questions/10124/the-factorial-of-1-2-3#comment14723_10129
		if( d_fac < 0)
			frac2 = 0.0;
		else
			frac2 = c/factorial(d_fac);
		return frac1*frac2;
	}


	// Evaluate the associated Legendre polynomial of degree @l and order @m at
	// coordinate @x. The inputs must satisfy:
	// 1. l >= 0
	// 2. 0 <= m <= l
	// 3. -1 <= x <= 1
	// See http://en.wikipedia.org/wiki/Associated_Legendre_polynomials
	//
	// This implementation is based off the approach described in [1],
	// instead of computing Pml(x) directly, Pmm(x) is computed. Pmm can be
	// lifted to Pmm+1 recursively until Pml is found
	//
	// note that the Condon-Shorteley Phase is included...
	double P(int l, int m, double x)
	{
		// Compute Pmm(x) = (-1)^m(2m - 1)!!(1 - x^2)^(m/2), where !! is the double
		// factorial.
		double pmm = 1.0;
		// P00 is defined as 1.0, do don't evaluate Pmm unless we know m > 0
		if (m > 0)
		{
			//double sign = (m % 2 == 0 ? 1 : -1);
			double sign = 1.0;
			pmm = sign * doubleFactorial(2 * m - 1) * pow(1 - x * x, m / 2.0);
		}

		if (l == m)
		{
			// Pml is the same as Pmm so there's no lifting to higher bands needed
			return pmm;
		}

		// Compute Pmm+1(x) = x(2m + 1)Pmm(x)
		double pmm1 = x * (2 * m + 1) * pmm;
		if (l == m + 1)
		{
			// Pml is the same as Pmm+1 so we are done as well
			return pmm1;
		}

		// Use the last two computed bands to lift up to the next band until l is
		// reached, using the recurrence relationship:
		// Pml(x) = (x(2l - 1)Pml-1 - (l + m - 1)Pml-2) / (l - m)
		for (int n = m + 2; n <= l; n++)
		{
			double pmn = (x * (2 * n - 1) * pmm1 - (n + m - 1) * pmm) / (n - m);
			pmm = pmm1;
			pmm1 = pmn;
		}
		// Pmm1 at the end of the above loop is equal to Pml
		return pmm1;
	}





	moexp::complex Y_gg( int l, int m, double theta, double phi )
	{
		return C(l,m)*P(l, m, std::cos(theta))*complex(std::cos(m*phi), std::sin(m*phi));
	}

	moexp::complex Y_cc( int l, int m, double theta, double phi )
	{
		return C(l,m)*P(l, m, std::cos(theta))*complex(std::cos(m*phi), -std::sin(m*phi));
	}

	complex Y( int l, int m, double theta, double phi )
	{
		if(m>=0)
		{
			return Y_gg(l, m, theta, phi);
		}else
		{
			return csp(m)*Y_cc(l, std::abs(m), theta, phi);
		}
	}

	// real spherical harmonics...
	// l must be >= 0
	// m must be >= -l and <= l
	double Y_real(int l, int m, double theta, double phi)
	{
		//CHECK(l >= 0, "l must be at least 0.");
		//CHECK(-l <= m && m <= l, "m must be between -l and l.");
		double kml = sqrt((2.0 * l + 1) * factorial(l - abs(m)) / (4.0 * M_PI * factorial(l + abs(m))));
		if (m > 0)
		{
			return sqrt(2.0)*kml*cos(m * phi)*csp(m)*P(l, m, cos(theta));
		}else
		if(m < 0)
		{
			return sqrt(2.0)*kml*sin(-m * phi)*csp(-m)*P(l, -m, cos(theta));
		}else
		{
			return kml*P(l, 0, cos(theta));
		}
	}

	std::unique_ptr<std::vector<complex>> project_Y(int order, const SphericalFunction<double>& func, int sample_count)
	{
		std::unique_ptr<std::vector<complex>> coeffs(new std::vector<complex>());
		coeffs->assign(numSHCoefficients(order), 0.0);

		RNGd rng;
		// integrate func times sh basis function to get the coefficient associated with that basis
		for( int i=0;i<sample_count;++i )
		{
			V3d d = sampleSphere<double>(rng);
			P2d theta_phi = sphericalCoordinates(d);
			double theta = theta_phi.x();
			double phi = theta_phi.y();
			double d_pdf = sampleSpherePDF();

			double f = func( theta_phi.x(), theta_phi.y() );

			// integrate and accumulate for each sh basis function
			for( int l=0;l<=order;++l )
				for( int m=-l;m<=l;++m )
				{
					complex sh = Y(l, m, theta, phi);
					complex sample = f*sh/d_pdf;
					complex& sh_coeff = (*coeffs)[shIndex(l, m)];
					sh_coeff += (sample - sh_coeff)/double(i+1);
				}
		}

		return coeffs;
	}

	complex Y_sum(int order, const std::vector<complex>& coeffs, double theta, double phi)
	{
		//CHECK(GetCoefficientCount(order) == coeffs.size(), "Incorrect number of coefficients provided.");
		complex sum = complex(0.0);//Zero<T>();
		for (int l = 0; l <= order; l++)
		{
			for (int m = -l; m <= l; m++)
			{
				sum += Y(l, m, theta, phi) * coeffs[shIndex(l, m)];
			}
		}
		return sum;
	}






	struct YTensorBuilder
	{

		std::vector<Tensor<complex>> compute_Y_tensors( int l, int m )
		{
			buildComponentCounts(l);

			// this will hold all the tensors by their rank
			std::map<int, Tensor<complex>> tensors;

			cas::Expression::Ptr expanded;
			cas::Scope scope;
			scope.m_variables["l"] = cas::num(l);

			// here we instantiate the variables l and m and expand the equation 2.7 in the Thorne paper
			if( m>= 0 )
			{
				scope.m_variables["m"] = cas::num(m);
				expanded = expand(expr_Y->deep_copy(), scope);
				if(l==3 && m==0)
					std::cout << expanded->toLatex() << std::endl;
			}else
			{
				scope.m_variables["m"] = cas::num(std::abs(m));
				expanded = expand(expr_Y_negative_m->deep_copy(), scope);
			}

			// now we extract the coefficients by analysing the resulting terms of the expanded equation
			cas::Addition::Ptr add = std::dynamic_pointer_cast<cas::Addition>(expanded);
			if(!add)
			{
				std::cout << "YTensorBuilder::compute_Y_tensors: error: equation expansion failed\n";
				throw std::runtime_error("hmghmhgm");
			}

			int numTerms = add->getNumOperands();
			for( int i=0;i<numTerms;++i )
			{
				cas::Multiplication::Ptr mul = std::dynamic_pointer_cast<cas::Multiplication>(add->getOperand(i));
				if(!mul)
				{
					std::cout << "YTensorBuilder::compute_Y_tensors: error: equation expansion failed\n";
					throw std::runtime_error("hmghmhgm");
				}
				//std::cout << "term " << i << "= " << mul->toLatex() << std::endl;

				// each term contributes to one tensor component
				double component_contribution = 1.0;

				// in each term, we expect l number of variables to occur, from
				// which we can figure out the component...this is called the component code
				V3i code(0,0,0);

				bool flag = false;


				int numFactors = mul->getNumOperands();
				for( int j=0;j<numFactors;++j )
				{
					cas::Expression::Ptr factor = mul->getOperand(j);

					cas::Variable::Ptr var = std::dynamic_pointer_cast<cas::Variable>(factor);
					cas::Power::Ptr pow = std::dynamic_pointer_cast<cas::Power>(factor);

					//std::cout << "\tfactor " << j << "= " << factor->toLatex() << std::endl;

					if( std::dynamic_pointer_cast<cas::Number>(factor) )
					{
						cas::Number::Ptr n = std::dynamic_pointer_cast<cas::Number>(factor);
						double value = 1.0;
						switch(n->getType())
						{
							case cas::Number::EInteger:value = n->get_int();break;
							case cas::Number::EReal:value = n->get_real();break;
							default:
							{
								std::cout << "YTensorBuilder::compute_Y_tensors: error: unable to handle number type\n";
								throw std::runtime_error("hrhhrhth");
							}break;
						};

						component_contribution *= value;
						//std::cout << "\tgot number  value=" << value << std::endl;
					}else
					// check for Clm or almj
					if( std::dynamic_pointer_cast<cas::Index>(factor) )
					{
						double index_contribution = 1.0;
						cas::Index::Ptr index = std::dynamic_pointer_cast<cas::Index>(factor);
						if(index->getBaseName() == "C")
							index_contribution *= C(index->getExponent(0), index->getExponent(1));
						else
						if(index->getBaseName() == "a")
							index_contribution *= a(index->getExponent(0), index->getExponent(1), index->getExponent(2));

						component_contribution *= index_contribution;
						//std::cout << "\tgot index " << index->getBaseName() << " value=" << index_contribution << std::endl;
					}else
					if(var)
					{
						if(var->getName() == "n_x")
							code[0] += 1;
						else
						if(var->getName() == "in_y")
							code[1] += 1;
						else
						if(var->getName() == "n_z")
							code[2] += 1;
					}else
					if(!var && pow)
					{
						var = std::dynamic_pointer_cast<cas::Variable>(pow->getOperand(0));
						cas::Number::Ptr exp = std::dynamic_pointer_cast<cas::Number>(pow->getOperand(1));
						if(!var)
						{
							std::cout << "YTensorBuilder::compute_Y_tensors: error: power of non variable encountered\n";
							throw std::runtime_error("dsfhdfhdfh");
						}
						if(!exp)
						{
							std::cout << "YTensorBuilder::compute_Y_tensors: error: power with non number exponent encountered\n";
							throw std::runtime_error("dsfhdfhdfh");
						}
						if(exp->getType() != cas::Number::EInteger)
						{
							std::cout << "YTensorBuilder::compute_Y_tensors: error: power with non integer exponent encountered\n";
							throw std::runtime_error("dsfhdfhdfh");
						}

						int exp_number = exp->get_int();
						if(exp_number >=0)
						{
							if(var->getName() == "n_x")
								code[0] += exp_number;
							else
							if(var->getName() == "in_y")
								code[1] += exp_number;
							else
							if(var->getName() == "n_z")
								code[2] += exp_number;
						}else
							flag = true;
					}else
					{
						std::cout << "YTensorBuilder::compute_Y_tensors: error: unable to handle factor type\n";std::flush(std::cout);
						std::cout << factor->toLatex() << std::endl;
						throw std::runtime_error("hrhhrhth");
					}

				} // for each factor of current term


				if(flag)
				{
					//std::cout << "YTensorBuilder::compute_Y_tensors: warning: variable with negative exponent found contribution=" << component_contribution << std::endl;
					continue;
				}

				// complete the final contribution associated with the current term
				complex final_contribution = std::pow( complex(0.0, 1.0), code[1] )*component_contribution;

				// now we need to find out which tensor component this coefficient is associated with
				// find the rank of the tensor to which the current component belongs
				int component_l = code[0] + code[1] + code[2];

				// look for the tensor in our map and create one if needed
				if( tensors.find(component_l) == tensors.end() )
					tensors[component_l] = Tensor<complex>(component_l);


				// we use the tensor associated with the rank which is implied by the current component code
				Tensor<complex>* tensor = &tensors[component_l];
				for( auto it = tensor->begin(), end = tensor->end();it!=end;++it )
				{
					if( it.code() == code )
					{
						it.value() += final_contribution/double(m_component_count[code]);
					}
				}


			} // for each term in the addition


			std::vector<Tensor<complex>> tensor_list;
			for( auto& it:tensors )
				tensor_list.push_back(it.second);
			return tensor_list;
		}

		static YTensorBuilder* instance()
		{
			if(!g_instance)
				g_instance = new YTensorBuilder();
			return g_instance;
		}


	private:
		YTensorBuilder():
			m_built_level(-1)
		{
			{
				using namespace cas;
				expr_Y = mul( index(var("C"), var("l"), var("m")), pow(add(var("n_x"), var("in_y")), var("m")),sum( "j", largest_integer(mul(pow(num(2), num(-1)) , add(var("l"), mul(num(-1), var("m"))))), mul( index( var("a"), var("l"), var("m"), var("j")), pow( var("n_z"), add(var("l"), mul( num(-1), var("m")), mul( num(-2), var("j")))))));
				expr_Y_negative_m = mul( pow(num(-1), var("m")), index(var("C"), var("l"), var("m")), pow(add(var("n_x"), mul(num(-1), var("in_y"))), var("m")),sum( "j", largest_integer(mul(pow(num(2), num(-1)) , add(var("l"), mul(num(-1), var("m"))))), mul( index( var("a"), var("l"), var("m"), var("j")), pow( var("n_z"), add(var("l"), mul( num(-1), var("m")), mul( num(-2), var("j")))))));
			}
		}

		struct compare_V3i
		{
			bool operator()(const V3i& a, const V3i& b) const
			{
				return std::make_tuple(a.x(), a.y(), a.z()) < std::make_tuple(b.x(), b.y(), b.z());
			}
		};

		void buildComponentCounts( int order )
		{
			while( m_built_level < order )
			{
				++m_built_level;
				int l = m_built_level;
				{
					Tensor<double> t(l);
					for( auto it = t.begin(), end = t.end(); it!=end;++it )
					{
						V3i code = it.code();
						if( m_component_count.find(code) == m_component_count.end() )
							m_component_count[code] = 0;
						++m_component_count[code];
					}
				}
			}
		}


		// this is used for extracting the tensor components from an expanded expression for Y
		cas::Expression::Ptr expr_Y;
		cas::Expression::Ptr expr_Y_negative_m;

		// this contains for each component the number of times it appears in its tensor
		// this is used for our way of dealing with symmetries for tensors of arbitrary rank
		std::map<V3i, int, compare_V3i> m_component_count;
		int m_built_level; // just a way to remember which component counts we already have built

		static YTensorBuilder* g_instance;
	};

	YTensorBuilder* YTensorBuilder::g_instance = 0;


	std::vector<Tensor<complex>> Y_tensors( int l, int m )
	{
		YTensorBuilder* builder = YTensorBuilder::instance();
		return builder->compute_Y_tensors( l, m );
	}
} // namespace moexp





void rasterizeSphericalFunctionSphere(const std::string& filename, moexp::SphericalFunction<Color3f> func, double exposure  )
{
	double scale = std::pow(2.0, exposure);
	houio::Geometry::Ptr geo = houio::Geometry::createSphere(120, 120, 1.0);
	houio::Attribute::Ptr pAttr = geo->getAttr("P");
	houio::Attribute::Ptr cdAttr = houio::Attribute::createV3f(pAttr->numElements());
	for( int i=0;i<pAttr->numElements();++i )
	{
		houio::math::V3f p = pAttr->get<houio::math::V3f>(i);
		P2d theta_phi = sphericalCoordinates<double>(V3d(p.x, p.y, p.z));
		double theta = theta_phi.x();
		double phi = theta_phi.y();
		Color3f col = func(theta, phi)*scale;
		cdAttr->set<houio::math::V3f>( i, houio::math::V3f(col.r(), col.g(), col.b()) );
	}
	geo->setAttr("Cd", cdAttr);
	houio::HouGeoIO::xport( filename, geo);
}

void rasterizeSphericalFunctionMap( EnvMap& envmap, moexp::SphericalFunction<Color3f> func )
{
	//std::ofstream f( "test_pixels_coords.txt", std::ios::binary | std::ios::trunc );
	int xres = envmap.bitmap().cols();
	int yres = envmap.bitmap().rows();
	for( int j=0;j<yres;++j )
		for( int i=0;i<xres;++i )
		{
			P2d xy( i+0.5f, j+0.5f );
			V3d d = envmap.xyToDirection(xy);
			P2d theta_phi = sphericalCoordinates(d);
			double theta = theta_phi.x();
			double phi = theta_phi.y();
			envmap.bitmap().coeffRef(j, i) = func(theta, phi);
			//f << theta << " "  << phi << std::endl;
		}
}

void rasterizeSphericalFunctionMap( const std::string& filename, moexp::SphericalFunction<Color3f> func )
{
	EnvMap envmap;
	//std::ofstream f( "test_pixels_coords.txt", std::ios::binary | std::ios::trunc );
	int xres = envmap.bitmap().cols();
	int yres = envmap.bitmap().rows();
	for( int j=0;j<yres;++j )
		for( int i=0;i<xres;++i )
		{
			P2d xy( i+0.5f, j+0.5f );
			V3d d = envmap.xyToDirection(xy);
			P2d theta_phi = sphericalCoordinates(d);
			double theta = theta_phi.x();
			double phi = theta_phi.y();
			envmap.bitmap().coeffRef(j, i) = func(theta, phi);
			//f << theta << " "  << phi << std::endl;
		}
	envmap.bitmap().saveEXR(filename);
}


// compute sh coefficients
void validate_moexp( int order )
{
	EnvMap map("envmap.exr");


	moexp::SphericalFunction<double> fun = [&](double theta, double phi) -> double
	{
		return map.eval(theta, phi).getLuminance();
	};

	rasterizeSphericalFunctionSphere( "groundtruth.bgeo", fun );

	// complex sh ---
	{
		// project
		int order = 30;
		int numSamples = 50000;
		std::unique_ptr<std::vector<moexp::complex>> sh_coeffs;
		sh_coeffs = moexp::project_Y(order, fun, numSamples);

		// reconstruct
		for( int l=0;l<=order;++l )
		{
			moexp::SphericalFunction<double> reconstruction = [&](double theta, double phi) -> double
			{
				// for some reason our complexSH reconstruction is flipped in y
				// when compared to real SH reconstruction
				V3d d = sphericalDirection(theta, phi);
				d.y() = -d.y();
				double theta2, phi2;
				P2d theta_phi2 = sphericalCoordinates(d);
				theta2 = theta_phi2.x();
				phi2 = theta_phi2.y();

				moexp::complex value = moexp::Y_sum(l, *sh_coeffs.get(), theta2, phi2);
				//return std::abs(value);
				return value.real();
			};

			std::string filename("testfit_complexsh_reconstruction_$0.bgeo");
			filename = replace(filename, "$0", toString(l));

			rasterizeSphericalFunctionSphere( filename, reconstruction );
		}
	}

	// real sh
	{
		// project
		int order = 30;
		int numSamples = 50000;
		std::unique_ptr<std::vector<double>> sh_coeffs;
		sh_coeffs = moexp::project_Y_real(order, fun, numSamples);

		// reconstruct
		for( int l=0;l<=order;++l )
		{
			moexp::SphericalFunction<double> reconstruction = [&](double theta, double phi) -> double
			{
				return moexp::Y_real_sum<double>(l, *sh_coeffs.get(), theta, phi);
			};
			std::string filename("testfit_realsh_reconstruction_$0.bgeo");
			filename = replace(filename, "$0", toString(l));
			rasterizeSphericalFunctionSphere( filename, reconstruction );
		}
	}

	// moment expansion ---
	{
		// produce tensors...
		std::vector<std::vector<moexp::Tensor<moexp::complex>>> Y_tensors_all;
		for( int l=0;l<order;++l )
		{
			std::cout << "building Y_tensor for order l=" << l << std::endl;
			for( int m=-l;m<=l;++m )
			{
				std::cout << "\tm=" << m << std::endl;
				std::vector<moexp::Tensor<moexp::complex>> tensors = moexp::Y_tensors(l, m);
				Y_tensors_all.push_back(tensors);
			}
		}

		// validate that comples SH basis functions are identical to our ylm tensors ---
		/*
		std::vector<Bitmap::Ptr> error_bitmaps;
		int Ylm_index = 0;
		for( int l=0;l<order;++l )
		{
			for( int m=-l;m<=l;++m, ++Ylm_index )
			{
				//if( (l != 1) && (m!=1))
				//	continue;
				std::vector<moexp::Tensor<moexp::complex>>& Y_tensors = Y_tensors_all[Ylm_index];

				EnvMap map_error;


				moexp::SphericalFunction<Color3f> error = [&](double theta, double phi) -> Color3f
				{
					V3d n = sphericalDirection(theta, phi);
					moexp::complex csh = moexp::Y(l, m, theta, phi);
					moexp::complex csh_from_tensor_contraction;

					for( auto& tensor:Y_tensors )
						csh_from_tensor_contraction += moexp::contract(tensor, n);

					double error = std::abs(csh-csh_from_tensor_contraction);

					return Color3f(error);
				};

				rasterizeSphericalFunction(map_error, error);

				{
					std::string filename("error_$0_$1.exr");
					filename = replace(filename, "$0", toString(l));
					filename = replace(filename, "$1", toString(m));
					map_error.bitmap().saveEXR(filename);
					error_bitmaps.push_back(std::make_shared<Bitmap>(filename));
				}
			}
		}
		compositeImages(error_bitmaps, "error_all.exr");
		*/


		// project
		int numSamples = 50000;
		std::unique_ptr<std::vector<moexp::complex>> sh_coeffs;
		sh_coeffs = moexp::project_Y(order, fun, numSamples);

		// we wedge this for every order
		for( int k=1;k<order;++k )
		{
			// compute F_kl tensors according to eq. 2.13a in [thorn80]
			std::vector<moexp::Tensor<moexp::complex>> F_list;

			for( int l=0;l<k;++l )
				F_list.push_back(moexp::Tensor<moexp::complex>(l));

			for( int l=0;l<k;++l )
				for( int m=-l;m<=l;++m )
				{
					int shindex = moexp::shIndex(l,m);
					std::vector<moexp::Tensor<moexp::complex>> Y_tensors = Y_tensors_all[shindex];
					for( auto& tensor : Y_tensors )
						F_list[tensor.getOrder()].multiplyAdd((*sh_coeffs)[shindex], tensor);
				}

			// reconstruct
			moexp::SphericalFunction<double> reconstruction = [&](double theta, double phi) -> double
			{
				V3d n = sphericalDirection(theta, phi);

				// for some reason our complexSH reconstruction is flipped in y
				// when compared to real SH reconstruction
				n.y() = -n.y();

				moexp::complex value(0.0, 0.0);
				for( auto& F:F_list )
					value += moexp::contract(F, n);

				return value.real();
			};

			std::string filename("testfit_moexp_reconstruction_$0.bgeo");
			filename = replace(filename, "$0", toString(k-1));
			rasterizeSphericalFunctionSphere( filename, reconstruction );
		}
	}
}


void EnvMap::saveGeo( const std::string& filename, double exposure )
{
	double scale = std::pow(2.0,exposure );
	houio::Geometry::Ptr geo = houio::Geometry::createSphere(120, 120, 1.0);
	houio::Attribute::Ptr pAttr = geo->getAttr("P");
	houio::Attribute::Ptr cdAttr = houio::Attribute::createV3f(pAttr->numElements());
	for( int i=0;i<pAttr->numElements();++i )
	{
		houio::math::V3f p = pAttr->get<houio::math::V3f>(i);
		P2d theta_phi = sphericalCoordinates<double>(V3d(p.x, p.y, p.z));
		double theta = theta_phi.x();
		double phi = theta_phi.y();
		Color3f col = this->eval(theta, phi)*scale;
		cdAttr->set<houio::math::V3f>( i, houio::math::V3f(col.r(), col.g(), col.b()) );
	}
	geo->setAttr("Cd", cdAttr);
	houio::HouGeoIO::xport( filename, geo);
}

