#include <util/moexp2d.h>











namespace moexp
{
namespace _2d
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

	int numFourierCoefficients( int order )
	{
		return 1+2*order;
	}

	int numTensorComponents( int order )
	{
		int dimension = 2;
		return ipow(dimension, order);
	}



	double integrate( const Function<double>& func, int numSamples )
	{
		double tmin = 0.0;
		double tmax = 2.0*M_PI;
		double P = tmax-tmin;
		double dt = P/double(numSamples);
		double sum = 0.0;
		for( int i=0;i<numSamples;++i )
		{
			double t = dt*i;
			sum += func(t)*dt;
		}
		return sum;
	}


	std::unique_ptr<std::vector<double>> projectFourier( int order, const Function<double>& func )
	{
		std::unique_ptr<std::vector<double>> coeffs(new std::vector<double>());
		coeffs->assign(numFourierCoefficients(order), 0.0);

		int numSamples = 10000;

		(*coeffs)[0] = INV_TWOPI*integrate( func, numSamples );
		for( int i=1;i<=order;++i )
		{
			(*coeffs)[1 + (i-1)*2+0] = INV_PI*integrate([&](double t)->double{return func(t)*std::cos(i*t);}, numSamples);
			(*coeffs)[1 + (i-1)*2+1] = INV_PI*integrate([&](double t)->double{return func(t)*std::sin(i*t);}, numSamples);
		}

		return coeffs;
	}

	double fourierSum( int order, double* coeffs, double t )
	{
		double result = coeffs[0];

		for( int i=1;i<=order;++i )
		{
			result += coeffs[1 + (i-1)*2+0]*std::cos(i*t);
			result += coeffs[1 + (i-1)*2+1]*std::sin(i*t);
		}

		return result;
	}

	double fourierBasis(int order, double* coeffs, double t)
	{
		if(order==0)
			return coeffs[0];
		return  coeffs[1 + (order-1)*2+0]*std::cos(order*t)+
				coeffs[1 + (order-1)*2+1]*std::sin(order*t);
	}


	std::vector<Tensor<double>> convertFourierCoefficientsToMoments(double* coeffs)
	{
		std::vector<Tensor<double>> result;

		// hardcoded for now, later we can do this more automated using the same approach as in
		// the 3d case

		// zero moment, m0 == a0
		Tensor<double> m0(0);
		m0.componentRef(0) = coeffs[0];
		result.push_back(m0);

		// first moment, m1 == (a1, b1)
		Tensor<double> m1(1);
		m1.componentRef(0) = coeffs[1];
		m1.componentRef(1) = coeffs[2];
		result.push_back(m1);

		// second moment, m2 == [[a2, b2], [b2, -a2]]
		std::cout << "coeffs[3]=" << coeffs[3] << std::endl;
		std::cout << "coeffs[4]=" << coeffs[4] << std::endl;
		Tensor<double> m2(2);
		m2.componentRef(0) = coeffs[3];
		m2.componentRef(1) = coeffs[4];
		m2.componentRef(2) = coeffs[4];
		m2.componentRef(3) = -coeffs[3];
		result.push_back(m2);

		return result;
	}
}
}
