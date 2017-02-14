#include <light.h>

#include <math/bbox.h>







namespace lights
{
	struct PointLight : public Light
	{
		//typedef std::shared_ptr<PointLight> Ptr;

		PointLight():
			m_position(0.0, 0.0, 0.0),
			m_power(1.0f, 1.0f, 1.0f)
		{
		}

		PointLight( const P3d& position )
			:m_position(position),
			 m_power(1.0f, 1.0f, 1.0f)
		{

		}

		virtual Color3f sample( LightSample& ls )const override
		{
			ls.p =m_position;
			ls.d = normalized( ls.refP-ls.p, ls.distance );
			Color3f radiance( m_power.r()*INV_FOURPI, m_power.g()*INV_FOURPI, m_power.b()*INV_FOURPI );
			double G = 1.0/(ls.distance*ls.distance);
			//TODO: take surface normal at reference location into account once we do surfaces too
			return radiance*G;
		}

	private:
		P3d m_position;
		Color3f m_power;
	};

	struct DirectionalLight : public Light
	{
		//typedef std::shared_ptr<PointLight> Ptr;

		DirectionalLight():
			m_direction(0.0, 0.0, 0.0),
			m_radiance(1.0f, 1.0f, 1.0f),
			m_distance(1.0)
		{
		}

		DirectionalLight( const V3d& direction, const Box3d& scene_bound )
			:m_direction(direction),
			 m_radiance(1.0f, 1.0f, 1.0f)
		{
			m_distance = scene_bound.getExtents().norm();
		}

		virtual Color3f sample( LightSample& ls )const override
		{
			ls.p =ls.refP - m_direction*m_distance;
			ls.d = m_direction;
			ls.distance = std::numeric_limits<double>::infinity();
			//NB: there is no inverse squared falloff due to the directional lightsource
			// being a delta distribution in solid angle domain
			//TODO: take surface normal at reference location into account once we do surfaces too
			return m_radiance;
		}

	private:
		P3d m_direction;
		double m_distance;
		Color3f m_radiance;
	};


	Light::Ptr point( const P3d& position )
	{
		return std::make_shared<PointLight>(position);
	}

	Light::Ptr directional( const P3d& direction, const Box3d& scene_bound )
	{
		return std::make_shared<DirectionalLight>(direction, scene_bound);
	}
}

