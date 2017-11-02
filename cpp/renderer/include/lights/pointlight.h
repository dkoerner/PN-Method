#pragma once
#include <light.h>
#include <math/bbox.h>

struct PointLight : public Light
{
	typedef std::shared_ptr<PointLight> Ptr;

	PointLight():
	    m_position(0.0, 0.0, 0.0),
	    m_power(1.0, 1.0, 1.0)
	{
	}

	PointLight( const P3d& position, const V3d& power )
	    :m_position(position),
	     m_power(power)
	{
	}

	virtual V3d sample( LightSample& ls )const override
	{
		ls.p =m_position;
		ls.d = normalized( ls.refP-ls.p, ls.distance );
		V3d radiance( m_power[0]*INV_FOURPI, m_power[1]*INV_FOURPI, m_power[2]*INV_FOURPI );
		double G = 1.0/(ls.distance*ls.distance);
		//TODO: take surface normal at reference location into account once we do surfaces too
		return radiance*G;
	}

	virtual std::string toString()const override
	{
		std::ostringstream ss;
		ss << "PointLight position=" << m_position.toString() << std::endl;
		return ss.str();
	}

private:
	P3d m_position;
	V3d m_power;
};

