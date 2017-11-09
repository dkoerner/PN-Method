#pragma once

#include <memory>

#include <math/common.h>
#include <math/vector.h>
#include <math/bbox.h>
#include <math/rng.h>








// datastructure for working with the solution from a PNSolver
// basically a voxelgrid with a constant number of elements per voxel
// TODO: switch everything to real only
struct PNSolution
{
	typedef std::shared_ptr<PNSolution> Ptr;


	struct SHVector
	{
		static double eval(double theta, double phi, double* coeffs, int order);
		static double& get( int l, int m, double* coeffs );
		static const double& get( int l, int m, const double* coeffs );
		static void staticInit();
		static void staticShutdown();
		static double legendreP(int l, int m, double x);
	private:
		static double factorial(int x);
		inline static double normalization(int l, int m)
		{
			if (l < m_normTableSize)
				return m_normalization[l*(l+1)/2 + m];
			else
				return computeNormalization(l, m);
		}
		static double computeNormalization(int l, int m);
		static double *m_normalization;
		static int m_normTableSize;
	};

	struct SHSampler
	{
		SHSampler(int bands, int depth=5);

		//
		// Warp a uniform sample in [0,1]^2 to one that is
		// approximately proportional to the specified function.
		//
		// The resulting sample will have spherical coordinates
		// [0,pi]x[0,2pi] and its actual PDF (which might be
		// slightly different from the function evaluated at the
		// sample, even if $f$ is a distribution) will be returned.
		//
		double warp(const double *coeffs, P2d &sample) const;

	//protected:
		~SHSampler();

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
		double integrate(int depth, int zBlock, int phiBlock, const double *coeffs) const;
		// integrates SH function over given block by summing the integrals of all its child blocks
		double integrateChilds(int depth, int i, int j, const double *coeffs) const;

		int indexofSmallestElement(double array[], int size);
	protected:
		int m_bands;
		int m_depth;
		double ***m_phiMap;
		double ***m_legendreMap;
		int m_dataSize;
		double *m_normalization;
	};

	PNSolution(const std::string& filename);
	PNSolution(int order, const V3i& resolution, const Box3d& bound, const std::complex<double> *data);

	// evaluation ---
	// evaluate <- evalute SH function for given direction at given worldspace position (does interpolation of SH coefficients)
	double eval(const P3d &pWS, const V3d &direction)const;
	std::complex<double> evalCoefficient(const P3d &pWS, int coeff_index)const;
	Eigen::VectorXd evalCoefficients(const P3d &pWS)const;

	// sampling ---
	// samples a direction at the given worldspace position
	V3d sample( const P3d& pWS, double& pdf, RNGd& rng );
	V3d sample( const P3d& pWS, double& pdf, const P2d& sample );


	// io ---
	void save( const std::string& filename );


	// space transforms ---
	P3d voxelToLocal(const P3d& pVS)const;
	P3d localToVoxel(const P3d& pLS)const;
	P3d localToWorld( const P3d& pLS )const;
	P3d worldToLocal( const P3d& pWS )const;
	P3d voxelToWorld( const P3d& pVS)const;
	P3d worldToVoxel( const P3d& pWS)const;
	bool contains_voxel( const P3i& voxel )const;

	// info ---
	int getOrder()const;
	int getNumCoeffs()const;
	V3i getResolution()const;
	P3d getBoundMin()const;
	P3d getBoundMax()const;
	std::complex<double>* data();
private:

	int getIndex( const V3i& voxel_coord ) const;
	int m_order; // SH order
	int m_numCoeffs;
	V3i m_resolution;
	Box3d m_bound;
	V3d m_extend;
	std::vector<std::complex<double>> m_data; // number of voxels*number of coefficients
	SHSampler m_shsampler;
};

