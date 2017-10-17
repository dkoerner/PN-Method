#pragma once

#include <complex>
#include <iostream>

#include <math/common.h>
#include <math/vector.h>
#include <math/ray.h>

#include <field/Field.h>



// Adapter for lambda functions to scalar fields
struct Function : public Field
{
	typedef std::shared_ptr<Function> Ptr;
	typedef std::complex<double> ValueType;
	typedef std::function<ValueType (const P3d&)> FunctionType;


	Function( FunctionType fun ):
		m_fun(fun)
	{
	}


	virtual std::complex<double> eval( const P3d& pWS )const override
	{
		return m_fun(pWS);
	}

	virtual Field::Ptr createRestricted()const override
	{
		throw std::runtime_error("restriction step not implemented for Function");
		return Ptr();
	}


private:
	FunctionType m_fun;
};

