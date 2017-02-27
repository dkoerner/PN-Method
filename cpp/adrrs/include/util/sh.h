#pragma once

#include <functional>
#include <vector>
#include <memory>





namespace sh
{

	template<typename T>
	using SphericalFunction = std::function<T (double, double)>;

	// Get the total number of coefficients for a function represented by
	// all spherical harmonic basis of degree <= @order (it is a point of
	// confusion that the order of an SH refers to its degree and not the order).
	int numCoefficients(int order);

	// Get the one dimensional index associated with a particular degree @l
	// and order @m. This is the index that can be used to access the Coeffs
	// returned by SHSolver.
	int index(int l, int m);

	double factorial(int x);
	double doubleFactorial(int x);

	// includes Condon-Shorteley Phase...
	double P(int l, int m, double x);

	// l must be >= 0
	// m must be >= -l and <= l
	double eval(int l, int m, double theta, double phi);

	template <typename T>
	T evalSum(int order, const std::vector<T>& coeffs, double theta, double phi)
	{
		//CHECK(GetCoefficientCount(order) == coeffs.size(), "Incorrect number of coefficients provided.");
		T sum = T(0.0);//Zero<T>();
		for (int l = 0; l <= order; l++)
		{
			for (int m = -l; m <= l; m++)
			{
				sum += eval(l, m, theta, phi) * coeffs[index(l, m)];
			}
		}
		return sum;
	}

	template<typename T>
	std::unique_ptr<std::vector<T>> project(int order, const SphericalFunction<T>& func, int sample_count)
	{
		std::unique_ptr<std::vector<T>> coeffs(new std::vector<T>());
		coeffs->assign(numCoefficients(order), 0.0);

		RNGd rng;
		// integrate func times sh basis function to get the coefficient associated with that basis
		for( int i=0;i<sample_count;++i )
		{
			V3d d = sampleSphere<double>(rng);
			P2d theta_phi = sphericalCoordinates(d);
			double theta = theta_phi.x();
			double phi = theta_phi.y();
			double d_pdf = sampleSpherePDF();

			T f = func( theta_phi.x(), theta_phi.y() );

			// integrate and accumulate for each sh basis function
			for( int l=0;l<=order;++l )
				for( int m=-l;m<=l;++m )
				{
					double sh = eval(l, m, theta, phi);
					T sample = f*sh/d_pdf;
					T& sh_coeff = (*coeffs)[index(l, m)];
					sh_coeff += (sample - sh_coeff)/float(i+1);
				}
		}

		return coeffs;
	}
} // namespace sh
