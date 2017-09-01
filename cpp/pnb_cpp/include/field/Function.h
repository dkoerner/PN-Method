#pragma once

#include <complex>
#include <iostream>

#include <math/common.h>
#include <math/vector.h>
#include <math/ray.h>

#include <field/Field.h>





struct Function : public Field
{
	typedef std::shared_ptr<Function> Ptr;
	typedef std::complex<double> ValueType;
	typedef std::function<ValueType (const P2d&)> FunctionType;


	Function( FunctionType fun ):
		m_fun(fun)
	{
	}


	virtual std::complex<double> eval( const P2d& pWS )const override
	{
		return m_fun(pWS);
	}


private:
	FunctionType m_fun;
};
