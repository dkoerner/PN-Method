#include <util/sh.h>


#include <math.h>






namespace sh
{

	// Get the total number of coefficients for a function represented by
	// all spherical harmonic basis of degree <= @order (it is a point of
	// confusion that the order of an SH refers to its degree and not the order).
	int numCoefficients(int order)
	{
		return (order + 1) * (order + 1);
	}

	// Get the one dimensional index associated with a particular degree @l
	// and order @m. This is the index that can be used to access the Coeffs
	// returned by SHSolver.
	int index(int l, int m)
	{
		return l * (l + 1) + m;
	}


	const double M_PI = 3.14159265358979323846;

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
	double EvalLegendrePolynomial(int l, int m, double x)
	{
		// Compute Pmm(x) = (-1)^m(2m - 1)!!(1 - x^2)^(m/2), where !! is the double
		// factorial.
		double pmm = 1.0;
		// P00 is defined as 1.0, do don't evaluate Pmm unless we know m > 0
		if (m > 0)
		{
			double sign = (m % 2 == 0 ? 1 : -1);
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

	// l must be >= 0
	// m must be >= -l and <= l
	double eval(int l, int m, double theta, double phi)
	{
		//CHECK(l >= 0, "l must be at least 0.");
		//CHECK(-l <= m && m <= l, "m must be between -l and l.");
		double kml = sqrt((2.0 * l + 1) * factorial(l - abs(m)) /
						(4.0 * M_PI * factorial(l + abs(m))));
		if (m > 0)
		{
			return sqrt(2.0) * kml * cos(m * phi) *
				EvalLegendrePolynomial(l, m, cos(theta));
		}else
		if(m < 0)
		{
			return sqrt(2.0) * kml * sin(-m * phi) *
				EvalLegendrePolynomial(l, -m, cos(theta));
		}else
		{
			return kml * EvalLegendrePolynomial(l, 0, cos(theta));
		}
	}





} // namespace sh
