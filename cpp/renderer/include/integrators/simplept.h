#pragma once
#include <integrator.h>

#include <scene.h>
#include <math/frame.h>

#include <volume.h>



// intersects the volume bound and starts volume path tracing within the volume
// will stop once the path exits the volume
struct SimplePT : public Integrator
{
	typedef std::shared_ptr<SimplePT> Ptr;

	SimplePT( int maxDepth = std::numeric_limits<int>::max(), bool doSingleScattering = true ):
	    Integrator(),
	    m_maxDepth(maxDepth),
	    m_doSingleScattering(doSingleScattering)
	{
	}


	V3d trace( TraceInfo& ti, RNGd& rng )const
	{
		V3d L(0.0f, 0.0f, 0.0f);

		if(ti.debug)
			std::cout <<"VolumePathTracer::trace\n";

		//L += ti.nee(rng, ti.debug);
		//return L;

		while( ti.depth < m_maxDepth )
		{
			if(ti.debug)
				std::cout <<"VolumePathTracer::trace depth=" << ti.depth << std::endl;

			// get next intersection ---
			if( !ti.intersect() )
				// intersection test failed
				return V3d(0.0, 0.0, 0.0);

			// get next interaction vertex ---
			if( !ti.propagate(rng) )
				// distance sampling failed
				return V3d(0.0, 0.0, 0.0);

			if(ti.current_vertex.getType() == Vertex::ESurface)
				// we intersected the volume boundary...lets quit
				return L;

			// perform next event estimation ---
			if( ((ti.depth == 0) && m_doSingleScattering ) || (ti.depth > 0) )
				L += nee(ti, rng);

			// handle scattering and find outgoing direction at next vertex ---
			if( !ti.scatter(rng) )
				// directional sampling failed
				return V3d(0.0, 0.0, 0.0);

			++ ti.depth;
		}

		return L;
	}

	virtual V3d Li( const Scene* scene, RadianceQuery& rq, RNGd& rng )const override;
	virtual V3d sample_transmittance( const Scene* scene, const Ray3d& ray, double maxt, RNGd& rng )const override
	{
		V3d sigma_t;
		double distance = delta_tracking( scene, ray, maxt, 0, rng, sigma_t );

		if( distance < maxt )
			// since we return 0, we dont need to produce the pdf in the denominator
			return V3d(0.0, 0.0, 0.0);

		//NB: transmittance sampling pdf cancels out with the transmittance term

		return V3d(1.0, 1.0, 1.0);
	}

	virtual std::string toString()const override
	{
		std::ostringstream ss;
		ss << "SimplePT singleScattering=" << m_doSingleScattering << " maxDepth=" << m_maxDepth << std::endl;
		return ss.str();
	}

private:
	bool m_doSingleScattering;
	int m_maxDepth;
};

