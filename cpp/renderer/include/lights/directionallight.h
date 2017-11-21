#pragma once
#include <light.h>









struct DirectionalLight : public Light
{
	typedef std::shared_ptr<DirectionalLight> Ptr;

	DirectionalLight():
	    m_direction(0.0, 0.0, 1.0),
	    m_radiance(1.0, 1.0, 1.0),
	    m_distance(1.0)
	{
	}

	DirectionalLight( const V3d& direction, const V3d& radiance )
	    :m_direction(direction),
	     m_radiance(radiance),
	     m_distance(1.0)
	{
	}

	virtual V3d sample( LightSample& ls )const override
	{
		ls.p = ls.refP - m_direction*m_distance;
		ls.d = m_direction;
		ls.distance = std::numeric_limits<double>::infinity();
		//NB: there is no inverse squared falloff due to the directional lightsource
		// being a delta distribution in solid angle domain
		//TODO: take surface normal at reference location into account once we do surfaces too
		return m_radiance;
	}

	virtual std::string toString()const override
	{
		std::ostringstream ss;
		ss << "DirectionalLight direction=" << m_direction.toString() << std::endl;
		return ss.str();
	}

private:
	V3d m_direction;
	V3d m_radiance;
	double m_distance;
};

