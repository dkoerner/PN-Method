#pragma once
#include <vector>
#include <functional>
#include <Eigen/Core>
#include <math/vector.h>



// Precompute normalization coefficients for the first 10 bands
#define SH_NORMTBL_SIZE 10
typedef V3d Vector;


struct SHVector;

double integrate_spherical_function(std::function<double(double, double)>& fun , double theta_start, double theta_end, double phi_start, double phi_end);

//
// * \brief Stores the diagonal blocks of a spherical harmonic
// * rotation matrix
// *
// * \ingroup libcore
// * \ingroup libpython
//
struct SHRotation
{
	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;

    std::vector<Matrix> blocks;

    /// Construct a new rotation storage for the given number of bands
	inline SHRotation(int bands) : blocks(bands)
	{
        for (int i=0; i<bands; ++i) {
            int dim = 2*i+1;
            blocks[i] = Matrix(dim, dim);
        }
    }

	//
	// * \brief Transform a coefficient vector and store the result into
	// * the given target vector.
	// *
	// * The source and target must have the same number of bands.
	//
    void operator()(const SHVector &source, SHVector &target) const;
};

// *
// * \brief Stores a truncated real spherical harmonics representation of
// * an L2-integrable function.
// *
// * Also provides some other useful functionality, such as evaluation,
// * projection and rotation.
// *
// * The Mathematica equivalent of the basis functions implemented here is:
// *
// * \code
// * SphericalHarmonicQ[l_, m_, \[Theta]_, \[Phi]_] :=
// *   Piecewise[{
// *      {SphericalHarmonicY[l, m, \[Theta], \[Phi]], m == 0},
// *      {Sqrt[2]*Re[SphericalHarmonicY[l, m, \[Theta], \[Phi]]], m > 0},
// *      {Sqrt[2]*Im[SphericalHarmonicY[l, -m, \[Theta], \[Phi]]], m < 0}
// *  }]
// * \endcode
// *
// * \ingroup libcore
// * \ingroup libpython
//
struct SHVector
{
public:
    /// Construct an invalid SH vector
    inline SHVector()
	    : m_bands(0)
	{
    }

    /// Construct a new SH vector (initialized to zero)
    inline SHVector(int bands)
	    : m_bands(bands), m_coeffs(bands*bands)
	{
        clear();
    }

    /// Unserialize a SH vector to a binary data stream
	//SHVector(Stream *stream);

    /// Copy constructor
    inline SHVector(const SHVector &v) : m_bands(v.m_bands),
	    m_coeffs(v.m_coeffs)
	{
    }

    /// Return the number of stored SH coefficient bands
	inline int getBands() const
	{
        return m_bands;
    }

    /// Serialize a SH vector to a binary data stream
	//void serialize(Stream *stream) const;

    /// Get the energy per band
	inline double energy(int band) const
	{
		double result = 0;
        for (int m=-band; m<=band; ++m)
            result += std::abs(operator()(band,m));
        return result;
    }

    /// Assignment
	inline SHVector &operator=(const SHVector &v)
	{
        m_bands = v.m_bands;
        m_coeffs = v.m_coeffs;
        return *this;
    }

    /// Set all coefficients to zero
	inline void clear()
	{
        m_coeffs.setZero();
    }

    /// Component-wise addition
	inline SHVector& operator+=(const SHVector &v)
	{
        ptrdiff_t extendBy = v.m_coeffs.size() - m_coeffs.size();
        if (extendBy > 0) {
            m_coeffs.conservativeResize(v.m_coeffs.rows());
            m_coeffs.tail(extendBy).setZero();
            m_bands = v.m_bands;
        }
        m_coeffs.head(m_coeffs.size()) += v.m_coeffs.head(m_coeffs.size());
        return *this;
    }

    /// Component-wise addition
	inline SHVector operator+(const SHVector &v) const
	{
        SHVector vec(std::max(m_bands, v.m_bands));
        if (m_bands > v.m_bands) {
            vec.m_coeffs = m_coeffs;
            vec.m_coeffs.head(v.m_coeffs.size()) += v.m_coeffs;
        } else {
            vec.m_coeffs = v.m_coeffs;
            vec.m_coeffs.head(m_coeffs.size()) += m_coeffs;
        }
        return vec;
    }

    /// Component-wise subtraction
	inline SHVector& operator-=(const SHVector &v)
	{
        ptrdiff_t extendBy = v.m_coeffs.size() - m_coeffs.size();
        if (extendBy > 0) {
            m_coeffs.conservativeResize(v.m_coeffs.rows());
            m_coeffs.tail(extendBy).setZero();
            m_bands = v.m_bands;
        }
        m_coeffs.head(m_coeffs.size()) -= v.m_coeffs.head(m_coeffs.size());
        return *this;
    }

    /// Component-wise subtraction
	inline SHVector operator-(const SHVector &v) const
	{
        SHVector vec(std::max(m_bands, v.m_bands));
        if (m_bands > v.m_bands) {
            vec.m_coeffs = m_coeffs;
            vec.m_coeffs.head(v.m_coeffs.size()) -= v.m_coeffs;
        } else {
            vec.m_coeffs = -v.m_coeffs;
            vec.m_coeffs.head(m_coeffs.size()) += m_coeffs;
        }
        return vec;
    }

    /// Add a scalar multiple of another vector
	inline SHVector& madd(double f, const SHVector &v)
	{
        ptrdiff_t extendBy = v.m_coeffs.size() - m_coeffs.size();
        if (extendBy > 0) {
            m_coeffs.conservativeResize(v.m_coeffs.rows());
            m_coeffs.tail(extendBy).setZero();
            m_bands = v.m_bands;
        }
        m_coeffs.head(m_coeffs.size()) += v.m_coeffs.head(m_coeffs.size()) * f;

        return *this;
    }

    /// Scalar multiplication
	inline SHVector &operator*=(double f)
	{
        m_coeffs *= f;
        return *this;
    }

    /// Scalar multiplication
	inline SHVector operator*(double f) const
	{
        SHVector vec(m_bands);
        vec.m_coeffs = m_coeffs * f;
        return vec;
    }

    /// Scalar division
	inline SHVector &operator/=(double f)
	{
		m_coeffs *= (double) 1 / f;
        return *this;
    }

    /// Scalar division
	inline SHVector operator/(double f) const
	{
        SHVector vec(m_bands);
        vec.m_coeffs = m_coeffs * (1/f);
        return vec;
    }

    /// Negation operator
	inline SHVector operator-() const
	{
        SHVector vec(m_bands);
        vec.m_coeffs = -m_coeffs;
        return vec;
    }

    /// Access coefficient m (in {-l, ..., l}) on band l
	inline double &operator()(int l, int m)
	{
        return m_coeffs[l*(l+1) + m];
    }

    /// Access coefficient m (in {-l, ..., l}) on band l
	inline const double &operator()(int l, int m) const
	{
        return m_coeffs[l*(l+1) + m];
    }

    /// Evaluate for a direction given in spherical coordinates
	double eval(double theta, double phi) const;

    /// Evaluate for a direction given in Cartesian coordinates
	double eval(const Vector &v) const;

	//
	// * \brief Evaluate for a direction given in spherical coordinates.
	// *
	// * This function is much faster but only works for azimuthally
	// * invariant functions
	//
	double evalAzimuthallyInvariant(double theta, double phi) const;

	//
	// * \brief Evaluate for a direction given in cartesian coordinates.
	// *
	// * This function is much faster but only works for azimuthally
	// * invariant functions
	//
	double evalAzimuthallyInvariant(const Vector &v) const;

    /// Check if this function is azumuthally invariant
    bool isAzimuthallyInvariant() const;

    /// Equality comparison operator
	inline bool operator==(const SHVector &v) const
	{
        return m_bands == v.m_bands && m_coeffs == v.m_coeffs;
    }

    /// Equality comparison operator
	inline bool operator!=(const SHVector &v) const
	{
        return !operator==(v);
    }

    /// Dot product
	inline friend double dot(const SHVector &v1, const SHVector &v2);

    /// Normalize so that the represented function becomes a valid distribution
    void normalize();
/*
    /// Compute the second spherical moment (analytic)
	//Matrix3x3 mu2() const;
*/

    /// Brute-force search for the minimum value over the sphere
	double findMinimum(int res) const;

    /// Add a constant value
	void addOffset(double value);
/*
	//
	// * \brief Convolve the SH representation with the supplied kernel.
	// *
	// * Based on the Funk-Hecke theorem -- the kernel must be rotationally
	// * symmetric around the Z-axis.
	//
    void convolve(const SHVector &kernel);

    /// Project the given function onto a SH basis (using a 2D composite Simpson's rule)
	template<typename Functor> void project(const Functor &f, int res = 32)
	{
        SAssert(res % 2 == 0);
		// Nested composite Simpson's rule
		double hExt = M_PI / res,
              hInt = (2*M_PI)/(res*2);

        for (int l=0; l<m_bands; ++l)
            for (int m=-l; m<=l; ++m)
                operator()(l,m) = 0;

		double *sinPhi = (double *) alloca(sizeof(double)*m_bands),
			  *cosPhi = (double *) alloca(sizeof(double)*m_bands);

        for (int i=0; i<=res; ++i) {
			double theta = hExt*i, cosTheta = std::cos(theta);
			double weightExt = (i & 1) ? 4.0f : 2.0f;
            if (i == 0 || i == res)
                weightExt = 1.0f;

            for (int j=0; j<=res*2; ++j) {
				double phi = hInt*j;
				double weightInt = (j & 1) ? 4.0f : 2.0f;
                if (j == 0 || j == 2*res)
                    weightInt = 1.0f;

                for (int m=0; m<m_bands; ++m) {
                    sinPhi[m] = std::sin((m+1)*phi);
                    cosPhi[m] = std::cos((m+1)*phi);
                }

				double value = f(sphericalDirection(theta, phi))*std::sin(theta)
                    * weightExt*weightInt;

                for (int l=0; l<m_bands; ++l) {
                    for (int m=1; m<=l; ++m) {
						double L = legendreP(l, m, cosTheta) * normalization(l, m);
                        operator()(l, -m) += value * SQRT_TWO * sinPhi[m-1] * L;
                        operator()(l, m) += value * SQRT_TWO * cosPhi[m-1] * L;
                    }

                    operator()(l, 0) += value * legendreP(l, 0, cosTheta) * normalization(l, 0);
                }
            }
        }

        for (int l=0; l<m_bands; ++l)
            for (int m=-l; m<=l; ++m)
                operator()(l,m) *= hExt*hInt/9;
    }

    /// Compute the relative L2 error
	template<typename Functor> double l2Error(const Functor &f, int res = 32) const
	{
        SAssert(res % 2 == 0);
		// Nested composite Simpson's rule
		double hExt = M_PI / res,
              hInt = (2*M_PI)/(res*2);
		double error = 0.0f, denom=0.0f;

        for (int i=0; i<=res; ++i) {
			double theta = hExt*i;
			double weightExt = (i & 1) ? 4.0f : 2.0f;
            if (i == 0 || i == res)
                weightExt = 1.0f;

            for (int j=0; j<=res*2; ++j) {
				double phi = hInt*j;
				double weightInt = (j & 1) ? 4.0f : 2.0f;
                if (j == 0 || j == 2*res)
                    weightInt = 1.0f;

				double value1 = f(sphericalDirection(theta, phi));
				double value2 = eval(theta, phi);
				double diff = value1-value2;
				double weight = std::sin(theta)*weightInt*weightExt;

                error += diff*diff*weight;
                denom += value1*value1*weight;
            }
        }

        return error/denom;
    }

    /// Turn into a string representation
    std::string toString() const;
	*/

    /// Return a normalization coefficient
	inline static double normalization(int l, int m)
	{
        if (l < SH_NORMTBL_SIZE)
            return m_normalization[l*(l+1)/2 + m];
        else
            return computeNormalization(l, m);
    }

	/*

	//
	// * \brief Recursively computes rotation matrices for each band of SH coefficients.
	// *
	// * Based on 'Rotation Matrices for Real Spherical Harmonics. Direct Determination by Recursion'
	// * by Ivanic and Ruedenberg. The implemented tables follow the notation in
	// * 'Spherical Harmonic Lighting: The Gritty Details' by Robin Green.
	//
    static void rotation(const Transform &t, SHRotation &rot);
	*/

    /// Precomputes normalization coefficients for the first few bands
    static void staticInitialization();

    /// Free the memory taken up by staticInitialization()
    static void staticShutdown();
	/*
protected:
    /// Helper function for rotation() -- computes a diagonal block based on the previous level
    static void rotationBlock(const SHRotation::Matrix &M1, const SHRotation::Matrix &Mp, SHRotation::Matrix &Mn);
	*/
    /// Compute a normalization coefficient
	static double computeNormalization(int l, int m);
public:
//private:

	int m_bands;
	Eigen::Matrix<double, Eigen::Dynamic, 1> m_coeffs;
	static double *m_normalization;
};

/*
inline double dot(const SHVector &v1, const SHVector &v2)
{
    const size_t size = std::min(v1.m_coeffs.size(), v2.m_coeffs.size());
    return v1.m_coeffs.head(size).dot(v2.m_coeffs.head(size));
}
*/


//
// * \brief Implementation of 'Importance Sampling Spherical Harmonics'
// * by W. Jarsz, N. Carr and H. W. Jensen (EUROGRAPHICS 2009)
// *
// * \ingroup libcore
// * \ingroup libpython
//
class SHSampler
{
public:

	//
	// * \brief Precompute a spherical harmonics sampler object
	// *
	// * \param bands Number of SH coefficient bands to support
	// * \param depth Number of recursive sample warping steps.
	//
    SHSampler(int bands, int depth);

	//
	// * \brief Warp a uniform sample in [0,1]^2 to one that is
	// * approximately proportional to the specified function.
	// *
	// * The resulting sample will have spherical coordinates
	// * [0,pi]x[0,2pi] and its actual PDF (which might be
	// * slightly different from the function evaluated at the
	// * sample, even if $f$ is a distribution) will be returned.
	//
	double warp(const SHVector &f, P2d &sample) const;

    /// Return information on the size of the precomputed tables
    std::string toString() const;

	Eigen::MatrixXd getBlocks(int depth, const SHVector &f, bool use_org)const;
	//MTS_DECLARE_CLASS()
protected:
    /// Virtual destructor
    virtual ~SHSampler();

	// Index into the assoc. legendre polynomial table
    inline int I(int l, int m) const { return l*(l+1)/2 + m; }

	// Index into the phi table
    inline int P(int m) const { return m + m_bands; }

	inline double lookupIntegral(int depth, int zBlock, int phiBlock, int l, int m) const
	{
		return -m_phiMap[depth][phiBlock][P(m)] * m_legendreMap[depth][zBlock][I(l, std::abs(m))];
    }
    /// Recursively compute assoc. legendre & phi integrals
	double *legendreIntegrals(double a, double b);
	double *phiIntegrals(double a, double b);

	// integrates SH function over given block by using precomputed integration tables
	double integrate(int depth, int zBlock, int phiBlock, const SHVector &f) const;
	// integrates SH function over given block by summing the integrals of all its child blocks
	double integrateChilds(int depth, int i, int j, const SHVector &f) const;
protected:
    int m_bands;
    int m_depth;
	double ***m_phiMap;
	double ***m_legendreMap;
    int m_dataSize;
	double *m_normalization;
};
