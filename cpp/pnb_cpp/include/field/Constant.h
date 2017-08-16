#pragma once

#include <complex>
#include <iostream>

#include <math/common.h>
#include <math/vector.h>
#include <math/ray.h>

#include <field/Field.h>





struct Constant : public Field
{
	typedef std::shared_ptr<Constant> Ptr;

	Constant( std::complex<double> value ):
		m_value(value)
	{
	}


	virtual std::complex<double> eval( const P2d& pWS )const override
	{
		return m_value;
	}
	virtual std::complex<double> dx(const P2d& pWS)const override
	{
		return 0.0;
	}
	virtual std::complex<double> dxdx(const P2d& pWS)const override
	{
		return 0.0;
	}
	virtual std::complex<double> dxdy(const P2d& pWS)const override
	{
		return 0.0;
	}
	virtual std::complex<double> dy(const P2d& pWS)const override
	{
		return 0.0;
	}
	virtual std::complex<double> dydy(const P2d& pWS)const override
	{
		return 0.0;
	}
	virtual std::complex<double> dydx(const P2d& pWS)const override
	{
		return 0.0;
	}
	virtual std::complex<double> dz(const P2d& pWS)const override
	{
		return 0.0;
	}


private:
	std::complex<double> m_value;
};
