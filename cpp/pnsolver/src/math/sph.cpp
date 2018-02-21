#include <math/common.h>
#include <math/sph.h>









namespace sph
{

	std::complex<double> Y_gg( int l, int m, double theta, double phi )
	{
		return  csp(m)*
				normalization(l, std::abs(m))*
				legendreP(l, m, std::cos(theta))*std::complex<double>(std::cos(m*phi), std::sin(m*phi));
	}

	std::complex<double> Y_cc( int l, int m, double theta, double phi )
	{
		return  csp(m)*
				normalization(l, std::abs(m))*
				legendreP(l, m, std::cos(theta))*std::complex<double>(std::cos(m*phi), -std::sin(m*phi));
	}



	std::complex<double> complex_basis( int l, int m, double theta, double phi )
	{
		if(m>=0)
		{
			return Y_gg(l, m, theta, phi);
		}else
		{
			return csp(m)*Y_cc(l, std::abs(m), theta, phi);
		}
	}



	int numCoeffs(int order)
	{
		return (order+1)*(order+1);
	}

	int index( int l, int m )
	{
		return l * (l + 1) + m;
	}



	double eval(double theta, double phi, const double* coeffs, int order)
	{
		int numBands = order+1;
		double result = 0;
		double cosTheta = std::cos(theta);
		double *sinPhi = (double *) alloca(sizeof(double)*numBands),
			  *cosPhi = (double *) alloca(sizeof(double)*numBands);

		for (int m=0; m<numBands; ++m)
		{
			sinPhi[m] = std::sin((m+1) * phi);
			cosPhi[m] = std::cos((m+1) * phi);
		}

		for (int l=0; l<numBands; ++l)
		{
			for (int m=1; m<=l; ++m)
			{
				double L = legendreP(l, m, cosTheta) * normalization(l, m);
				result += get(l, -m, coeffs) * SQRT_TWO * sinPhi[m-1] * L;
				result += get(l, m, coeffs)  * SQRT_TWO * cosPhi[m-1] * L;
			}

			result += get(l, 0, coeffs) * legendreP(l, 0, cosTheta) * normalization(l, 0);
		}
		return result;
	}

	/*
	double basis( double theta, double phi, int l, int m )
	{
		double cosTheta = std::cos(theta);
		if(m<0)
		{
			return SQRT_TWO * std::sin(-m * phi) * legendreP(l, -m, cosTheta) * normalization(l, -m);
		}else
		if(m==0)
		{
			return legendreP(l, 0, cosTheta) * normalization(l, 0);
		}else
		if(m>0)
		{
			return SQRT_TWO * std::cos(m * phi) * legendreP(l, m, cosTheta) * normalization(l, m);
		}
		throw std::runtime_error("unexpected");
	}
	*/

	double& get( int l, int m, double* coeffs )
	{
		return coeffs[l*(l+1) + m];
	}

	const double& get( int l, int m, const double* coeffs )
	{
		return coeffs[l*(l+1) + m];
	}

	double legendreP(int l, int m, double x)
	{
		double p_mm = 1;

		if (m > 0) {
			double somx2 = std::sqrt((1 - x) * (1 + x));
			double fact = 1;
			for (int i=1; i<=m; i++) {
				p_mm *= (-fact) * somx2;
				fact += 2;
			}
		}

		if (l == m)
			return p_mm;

		double p_mmp1 = x * (2*m + 1) * p_mm;
		if (l == m+1)
			return p_mmp1;

		double p_ll = 0;
		for (int ll=m+2; ll <= l; ++ll) {
			p_ll = ((2*ll-1)*x*p_mmp1 - (ll+m-1) * p_mm) / (ll-m);
			p_mm = p_mmp1;
			p_mmp1 = p_ll;
		}

		return p_ll;
	}

	/*
	void convolve(const double *coeffs, const double *coeffs_filter, int order, double *coeffs_result)
	{
		for (int l=0; l<=order; ++l)
		{
			double alpha = std::sqrt(4 * (double) M_PI / (2*l + 1));
			for (int m=-l; m<=l; ++m)
				get(l, m, coeffs_result) *= alpha * get(l, 0, coeffs_filter);
		}
	}
	*/

	void convolve(double *coeffs, const double *coeffs_filter, int order)
	{
		for (int l=0; l<=order; ++l)
		{
			double alpha = std::sqrt(4 * (double) M_PI / (2*l + 1));
			for (int m=-l; m<=l; ++m)
				get(l, m, coeffs) *= alpha * get(l, 0, coeffs_filter);
		}
	}

	void project( std::function<double(double, double)> fun, int order, int res, double* coeffs_result)
	{
		int numBands = order+1;

		//SAssert(res % 2 == 0);
		// Nested composite Simpson's rule
		double  hExt = M_PI / res,
				hInt = (2*M_PI)/(res*2);

		for (int l=0; l<numBands; ++l)
			for (int m=-l; m<=l; ++m)
				get(l, m, coeffs_result) = 0.0;

		double *sinPhi = (double *) alloca(sizeof(double)*numBands),
			   *cosPhi = (double *) alloca(sizeof(double)*numBands);

		for (int i=0; i<=res; ++i)
		{
			double theta = hExt*i, cosTheta = std::cos(theta);
			double weightExt = (i & 1) ? 4.0f : 2.0f;
			if (i == 0 || i == res)
				weightExt = 1.0f;

			for (int j=0; j<=res*2; ++j)
			{
				double phi = hInt*j;
				double weightInt = (j & 1) ? 4.0f : 2.0f;
				if (j == 0 || j == 2*res)
					weightInt = 1.0f;

				for (int m=0; m<numBands; ++m)
				{
					sinPhi[m] = std::sin((m+1)*phi);
					cosPhi[m] = std::cos((m+1)*phi);
				}

				double value = fun(theta, phi)*std::sin(theta)*weightExt*weightInt;

				for (int l=0; l<numBands; ++l)
				{
					for (int m=1; m<=l; ++m)
					{
						double L = legendreP(l, m, cosTheta) * normalization(l, m);
						get(l, -m, coeffs_result) += value * SQRT_TWO * sinPhi[m-1] * L;
						get(l, m, coeffs_result) += value * SQRT_TWO * cosPhi[m-1] * L;
					}

					get(l, 0, coeffs_result) += value * legendreP(l, 0, cosTheta) * normalization(l, 0);
				}
			}
		}

		for (int l=0; l<numBands; ++l)
			for (int m=-l; m<=l; ++m)
			{
				get(l,m, coeffs_result) *= hExt*hInt/9;
			}
	}


	double computeNormalization(int l, int m)
	{
		//SAssert(m>=0);
		return std::sqrt( ((2*l+1) * factorial(l-m))/
						  (4 * (double) M_PI * factorial(l+m)));
	}




	// Compute the factorial for an integer @x. It is assumed x is at least 0.
	// This implementation precomputes the results for low values of x, in which
	// case this is a constant time lookup.
	//
	// The vast majority of SH evaluations will hit these precomputed values.
	double factorial(int x)
	{
		// Number of precomputed factorials and double-factorials that can be
		// returned in constant time.
		const int kCacheSize = 16;

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

	// condon-shortley phase
	double csp( int m )
	{
		return (m % 2 == 0 ? 1.0 : -1.0);
	}

	// lambda_l for the convolution thing. the eigenvalue associated with the eigenfunction Y_lm
	double lambda(int l)
	{
		return std::sqrt(4.0*M_PI/(2*l+1));
	}

	Eigen::MatrixXcd buildComplexToRealConversionMatrix( int order )
	{
		int nCoeffs = numCoeffs(order);

		Eigen::MatrixXcd S( nCoeffs, nCoeffs );
		S.fill(std::complex<double>(0.0, 0.0));

		for( int l=0;l<=order;++l )
		{
			for( int m=-l;m<=l;++m )
			{
				int i = index(l, m);
				if( m < 0 )
				{
					S.coeffRef(i, index(l, m)) = std::complex<double>(0.0, 1.0)/std::sqrt(2.0);
					S.coeffRef(i, index(l, -m)) = -std::complex<double>(0.0, 1.0)/std::sqrt(2.0)*csp(m);
				}else
				if( m == 0 )
				{
					S.coeffRef(i, index(l, 0)) = 1.0;
				}else
				if( m > 0 )
				{
					S.coeffRef(i, index(l, -m)) = 1.0/std::sqrt(2.0);
					S.coeffRef(i, index(l, m)) = 1.0/std::sqrt(2.0)*csp(m);
				}
			}
		}

		return S;
	}


	double basis( int l, int m, double theta, double phi )
	{
		double cosTheta = std::cos(theta);
		double L = legendreP(l, std::abs(m), cosTheta) * normalization(l, std::abs(m));

		if(m<0)
			return SQRT_TWO * std::sin(-m * phi) * L;
		else
		if(m==0)
			return L;
		else
		if(m>0)
			return SQRT_TWO * std::cos(m * phi) * L;
		throw std::runtime_error("asdasdsd");
	}


	void staticInit()
	{
		m_normalization = new double[m_normTableSize*(m_normTableSize+1)/2];
		for (int l=0; l<m_normTableSize; ++l)
			for (int m=0; m<=l; ++m)
				m_normalization[l*(l+1)/2 + m] = computeNormalization(l, m);

	}

	void staticShutdown()
	{
		delete[] m_normalization;
		m_normalization = NULL;
	}

	double *m_normalization = NULL;
	int m_normTableSize = 10;







}
