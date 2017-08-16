#pragma once

#include <complex>
#include <iostream>

#include <math/common.h>
#include <math/vector.h>
#include <math/ray.h>
#include <math/sph.h>

#include <field/Field.h>
#include <field/Constant.h>





struct SHEXP : public RadianceField
{
	typedef std::shared_ptr<SHEXP> Ptr;


	static int num_sh_coeffs( int order )
	{
		return (order + 1) * (order + 1);
	}
	static int sh_index(int l, int m)
	{
		//if l<0 or m < -l or m > l:
		//		return None
		return l * (l + 1) + m;
	}



	SHEXP( int order ):
	    RadianceField(),
	    m_order(order)
	{
		int num_coeffs = num_sh_coeffs(order);
		m_coefficient_fields.resize(num_coeffs);
		for( int i=0;i<num_coeffs;++i )
			set_coeff(i, std::make_shared<Constant>(0.0));
	}

	void set_coeff( int index, Field::Ptr field )
	{
		m_coefficient_fields[index] = field;
	}

	Field::Ptr get_coeff( int index )
	{
		return m_coefficient_fields[index];
	}

	virtual std::complex<double> eval( const P2d& pWS, const V3d& omega )const override
	{
		P2d theta_phi = sphericalCoordinates(omega);
		double theta = theta_phi[0];
		double phi = theta_phi[1];
		std::complex<double> result = 0.0;
		for(int l=0;l<m_order+1;++l)
			for( int m=-l;m<=l;++m )
				result+=m_coefficient_fields[sh_index(l,m)]->eval(pWS)*sph::sph_basis(l,m,theta, phi);
		return result;
	}
	virtual std::complex<double> dx(const P2d& pWS, const V3d& omega)const override
	{
		P2d theta_phi = sphericalCoordinates(omega);
		double theta = theta_phi[0];
		double phi = theta_phi[1];
		std::complex<double> result = 0.0;
		for(int l=0;l<m_order+1;++l)
			for( int m=-l;m<=l;++m )
				result+=m_coefficient_fields[sh_index(l,m)]->dx(pWS)*sph::sph_basis(l,m,theta, phi);
		return result;
	}
	virtual std::complex<double> dxdx(const P2d& pWS, const V3d& omega)const override
	{
		P2d theta_phi = sphericalCoordinates(omega);
		double theta = theta_phi[0];
		double phi = theta_phi[1];
		std::complex<double> result = 0.0;
		for(int l=0;l<m_order+1;++l)
			for( int m=-l;m<=l;++m )
				result+=m_coefficient_fields[sh_index(l,m)]->dxdx(pWS)*sph::sph_basis(l,m,theta, phi);
		return result;
	}
	virtual std::complex<double> dxdy(const P2d& pWS, const V3d& omega)const override
	{
		P2d theta_phi = sphericalCoordinates(omega);
		double theta = theta_phi[0];
		double phi = theta_phi[1];
		std::complex<double> result = 0.0;
		for(int l=0;l<m_order+1;++l)
			for( int m=-l;m<=l;++m )
				result+=m_coefficient_fields[sh_index(l,m)]->dxdy(pWS)*sph::sph_basis(l,m,theta, phi);
		return result;
	}
	virtual std::complex<double> dy(const P2d& pWS, const V3d& omega)const override
	{
		P2d theta_phi = sphericalCoordinates(omega);
		double theta = theta_phi[0];
		double phi = theta_phi[1];
		std::complex<double> result = 0.0;
		for(int l=0;l<m_order+1;++l)
			for( int m=-l;m<=l;++m )
				result+=m_coefficient_fields[sh_index(l,m)]->dy(pWS)*sph::sph_basis(l,m,theta, phi);
		return result;
	}
	virtual std::complex<double> dydy(const P2d& pWS, const V3d& omega)const override
	{
		P2d theta_phi = sphericalCoordinates(omega);
		double theta = theta_phi[0];
		double phi = theta_phi[1];
		std::complex<double> result = 0.0;
		for(int l=0;l<m_order+1;++l)
			for( int m=-l;m<=l;++m )
				result+=m_coefficient_fields[sh_index(l,m)]->dydy(pWS)*sph::sph_basis(l,m,theta, phi);
		return result;
	}
	virtual std::complex<double> dydx(const P2d& pWS, const V3d& omega)const override
	{
		P2d theta_phi = sphericalCoordinates(omega);
		double theta = theta_phi[0];
		double phi = theta_phi[1];
		std::complex<double> result = 0.0;
		for(int l=0;l<m_order+1;++l)
			for( int m=-l;m<=l;++m )
				result+=m_coefficient_fields[sh_index(l,m)]->dydx(pWS)*sph::sph_basis(l,m,theta, phi);
		return result;
	}
	virtual std::complex<double> dz(const P2d& pWS, const V3d& omega)const override
	{
		P2d theta_phi = sphericalCoordinates(omega);
		double theta = theta_phi[0];
		double phi = theta_phi[1];
		std::complex<double> result = 0.0;
		for(int l=0;l<m_order+1;++l)
			for( int m=-l;m<=l;++m )
				result+=m_coefficient_fields[sh_index(l,m)]->dz(pWS)*sph::sph_basis(l,m,theta, phi);
		return result;
	}
	virtual std::complex<double> integral_over_solid_angle(const P2d& pWS)const override
	{
		return m_coefficient_fields[0]->eval(pWS)*std::sqrt(4.0*M_PI);
	}


private:
	int m_order;
	std::vector<Field::Ptr> m_coefficient_fields;
};
