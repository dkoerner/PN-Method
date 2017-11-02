#pragma once

#include <complex>



namespace sph
{
    typedef std::complex<double> complex;
	inline int numSHCoeffs( int order )
	{
		return (order+1)*(order+1);
	}

	inline int shIndex( int l, int m )
	{
		return l * (l + 1) + m;
	}

    complex sph_basis( int l, int m, double theta, double phi );
}

