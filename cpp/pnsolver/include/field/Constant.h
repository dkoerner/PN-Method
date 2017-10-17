#pragma once

#include <complex>
#include <iostream>

#include <math/common.h>
#include <math/vector.h>
#include <math/ray.h>

#include <field/Field.h>



// A scalarfield which is constant throughout
struct Constant : public Field
{
	typedef std::shared_ptr<Constant> Ptr;

	Constant( std::complex<double> value ):
		m_value(value)
	{
	}


	virtual std::complex<double> eval( const P3d& pWS )const override
	{
		return m_value;
	}


	virtual Field::Ptr createRestricted()const override
	{
		return std::make_shared<Constant>(m_value);
	}

private:
	std::complex<double> m_value;
};

