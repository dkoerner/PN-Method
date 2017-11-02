#pragma once
#include <memory>
#include <math/vector.h>

// Abstract interface to a scalarfield which
// can be evaluated at worldspace positions
template<typename T>
struct Field
{
	typedef std::shared_ptr<Field<T>> Ptr;
	virtual T eval( const P3d& pLS )const=0; // NB: here field uses local space

	// this is required for delta tracking where we need to know the maximum
	virtual std::pair<T, T> getValueRange()const=0;
};

typedef Field<double> Fieldd;
typedef Field<V3d> Field3d;



// A scalarfield which is constant throughout
template<typename T>
struct ConstantField : public Field<T>
{
	typedef std::shared_ptr<ConstantField<T>> Ptr;

	ConstantField( T value ):
	    m_value(value)
	{
	}


	virtual T eval( const P3d& pLS )const override
	{
		return m_value;
	}
	virtual std::pair<T, T> getValueRange()const
	{
		return std::make_pair(m_value, m_value);
	}

private:
	T m_value;
};

typedef ConstantField<double> ConstantFieldd;
typedef ConstantField<V3d> ConstantField3d;



