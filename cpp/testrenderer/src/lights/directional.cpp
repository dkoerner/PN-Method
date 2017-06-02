#include <light.h>

#include <math/bbox.h>




namespace lights
{
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

	Light::Ptr directional( const P3d& direction, const Box3d& scene_bound )
	{
		return std::make_shared<DirectionalLight>(direction, scene_bound);
	}
} // namespace light
