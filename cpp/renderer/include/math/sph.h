#pragma once

#include <complex>
#include <functional>



namespace sph
{
/*
    typedef std::complex<double> complex;
	inline int numSHCoeffs( int order )
	{
		return (order+1)*(order+1);
	}

	inline int shIndex( int l, int m )
	{
		return l * (l + 1) + m;
	}

	double P(int l, int m, double x);

    complex sph_basis( int l, int m, double theta, double phi );
*/







	int numCoeffs(int order);
	int index( int l, int m );
	double eval(double theta, double phi, const double *coeffs, int order);
	double& get( int l, int m, double* coeffs );
	const double& get( int l, int m, const double* coeffs );
	void staticInit();
	void staticShutdown();
	double legendreP(int l, int m, double x);
	//void convolve( const double* coeffs, const double* coeffs_filter, int order, double* coeffs_result );
	void convolve( double* coeffs, const double* coeffs_filter, int order );
	void project(std::function<double(double, double)> fun, int order, int res, double* coeffs_result );
	double csp( int m );
	double lambda( int l );
	Eigen::MatrixXcd buildComplexToRealConversionMatrix( int order );
	double basis( int l, int m, double theta, double phi );
	std::complex<double> complex_basis( int l, int m, double theta, double phi );
	double factorial(int x);
	double computeNormalization(int l, int m);
	extern double *m_normalization;
	extern int m_normTableSize;
	inline static double normalization(int l, int m)
	{
		if (l < m_normTableSize)
			return m_normalization[l*(l+1)/2 + m];
		else
			return computeNormalization(l, m);
	}
}

