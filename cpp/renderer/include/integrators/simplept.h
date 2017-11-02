#include <integrator.h>

#include <scene.h>
#include <math/frame.h>

#include <volume.h>


// returns sigma_t at sampled position (is invalid when we exceeded maxt)
double delta_tracking( const Scene* scene, const Ray3d& ray, double maxt, int component, RNGd& rng, V3d& sigma_t, bool debug = false );

struct Vertex
{
	enum EVertexType
	{
		EInvalid = 0,
		ESurface = 1,
		EVolume = 2
	};

	void setPosition( const P3d& position, V3d sigma_t, V3d albedo )
	{
		m_p = position;
		m_type = EVolume;
		m_sigma_t = sigma_t;
		m_albedo = albedo;
		m_sigma_s = sigma_t.cwiseProduct(albedo);
	}
	void setPosition( const P3d& position, const Framed& frame)
	{
		m_p = position;
		m_type = ESurface;
	}

	const P3d& getPosition()const
	{
		return m_p;
	}

	EVertexType getType()const
	{
		return m_type;
	}

//		private:
	P3d m_p; // position in worldspace
	EVertexType m_type; // volume or surface scattering point
	V3d m_sigma_t;
	V3d m_sigma_s;
	V3d m_albedo;
};



struct TraceInfo
{
	Vertex current_vertex;
	V3d current_direction;
	double current_distance;
	V3d throughput_over_pdf;
	Intersection its;
	int depth;
	const Scene* scene;
	bool debug;

	TraceInfo():
	    throughput_over_pdf(1.0, 1.0, 1.0),
	    debug(false)
	{
	}



	bool intersect()
	{
		return scene->intersect(Ray3d( current_vertex.getPosition(), current_direction ), its);
	}



	bool propagate(RNGd& rng)
	{
		Ray3d ray( current_vertex.getPosition(), current_direction );
		V3d sigma_t;
		double distance = delta_tracking( scene, ray, its.t, 0, rng, sigma_t, debug);

		if( distance < its.t )
		{
			// volume scattering event ---
			P3d pWS = ray(distance);
			current_vertex.setPosition( pWS, sigma_t, scene->volume->evalAlbedo(pWS));

			// apply pdf and volume attenuation
			throughput_over_pdf[0] *= 1.0/sigma_t.x();
			throughput_over_pdf[1] *= 1.0/sigma_t.y();
			throughput_over_pdf[2] *= 1.0/sigma_t.z();

			if( debug )
			{
				std::cout << "TraceInfo::propagate: got scattering event!\n";
			}
		}else
		{
			// surface scattering event ---
			current_vertex.setPosition( its.p, Framed() );
		}

		// note that the transmittance term in the delta tracking pdf cancels
		// out with the transmittance term of our volume attenuation

		// further note that we assume the direction pdf is in solid angle units and we
		// need to convert the pdf of our vertex from solid angle times distance to volume.
		// This conversion cancels out with the geometry factor in throughput.

		// TODO: apply cosine terms and handle bounded volume (where no cosine is required)

		return true;
	}

	bool scatter(RNGd& rng)
	{
		if( current_vertex.getType() == Vertex::EVolume )
		{
			V3d wi = current_direction;
			double pdf;
			double phase_over_pdf = scene->volume->samplePhase(current_vertex.getPosition(), wi, current_direction, pdf, rng);
			throughput_over_pdf = phase_over_pdf*throughput_over_pdf.cwiseProduct(current_vertex.m_sigma_s);
			return true;
		}else
		{
			if( debug )
			{
				std::cout << "TraceInfo::scatter: surface sampling not implemented!\n";
			}
			throw std::runtime_error("surface sampling not implemented");
		}
		return false;
	}



	V3d nee( RNGd& rng, bool debug )
	{
		// sample light source and transmittance ---
		LightSample ls;
		ls.refP = current_vertex.getPosition();
		V3d attenuated_light_over_pdf = scene->sample_attenuated_directlight( ls, rng, debug );
		//std::cout << "attenuated_light_over_pdf=" <<  attenuated_light_over_pdf.toString() << std::endl;
		//return attenuated_light_over_pdf;

		// apply scattering ---
		V3d phase_times_sigma_s = scene->volume->evalPhase( current_vertex.getPosition(), -ls.d, -current_direction )*current_vertex.m_sigma_s;

		//std::cout << "phase_times_sigma_s=" << phase_times_sigma_s.toString() << std::endl;
		//std::cout << "throughput_over_pdf=" << throughput_over_pdf.toString() << std::endl;
		//std::cout << "current_vertex.m_sigma_s=" << current_vertex.m_sigma_s.toString() << std::endl;

		return attenuated_light_over_pdf.cwiseProduct(phase_times_sigma_s).cwiseProduct(throughput_over_pdf);
	}
};

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
				L += ti.nee(rng, ti.debug);

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

