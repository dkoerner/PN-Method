#pragma once

#include <complex>
#include <iostream>

#include <math/common.h>
#include <math/vector.h>
#include <math/ray.h>
#include <math/sph.h>

#include <field/Field.h>
#include <field/Constant.h>








// This class implements the RadianceField interface and uses
// an underlying SH representation, where each coefficient is
// given as a (potentially spatially varying) scalarfield
struct SHEXP : public RadianceField
{
	typedef std::shared_ptr<SHEXP> Ptr;


	SHEXP( int order ):
	    RadianceField(),
		m_coeffs(order)
	{
	}


	void setCoefficientField( int l, int m, Field::Ptr field )
	{
		m_coeffs.setField(l, m, field);
	}

	Field::Ptr getCoefficientField( int l, int m, int index )
	{
		return m_coeffs.getField(l,m);
	}

	virtual double eval( const P3d& pWS, const V3d& omega )const override
	{
		/*
		P2d theta_phi = sphericalCoordinates(omega);
		double theta = theta_phi[0];
		double phi = theta_phi[1];
		std::complex<double> result = 0.0;
		for(int l=0;l<m_coeffs.getOrder()+1;++l)
			for( int m=-l;m<=l;++m )
				result+=m_coeffs.eval(l,m,pWS)*sph::basis(l,m,theta, phi);
		return result.real();
		*/
		throw std::runtime_error("not implemented");
		return 0.0;
	}

private:
	SHCoefficientFieldArray m_coeffs;
};
