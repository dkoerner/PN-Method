#pragma once

#include <complex>



namespace sph
{
    typedef std::complex<double> complex;
    complex sph_basis( int l, int m, double theta, double phi );
}

