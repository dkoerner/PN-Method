#include <integrator.h>

#include <scene.h>
#include <math/frame.h>


namespace integrators
{
	// returns sigma_t at sampled position (is invalid when we exceeded maxt)
	double delta_tracking( const Scene* scene, const Ray3d& ray, double maxt, int component, RNGd& rng, Color3f& sigma_t )
	{
		double sigma_t_max = scene->volume->getMaxExtinction()[component];

		double t = 0.0;
		while(true)
		{
			double step = -log( 1.0-rng.next1D() )/sigma_t_max;
			t += step;

			if(t>= maxt)
				break;

			sigma_t = scene->volume->evalExtinction(ray(t));

			// russian roulette
			if(rng.next1D()<sigma_t[component]/sigma_t_max)
				break;
		}

		return t;
	}

	struct Vertex
	{
		enum EVertexType
		{
			EInvalid = 0,
			ESurface = 1,
			EVolume = 2
		};

		void setPosition( const P3d& position, Color3f sigma_t, Color3f albedo )
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
		Color3f m_sigma_t;
		Color3f m_sigma_s;
		Color3f m_albedo;
	};



	struct TraceInfo
	{
		Vertex current_vertex;
		V3d current_direction;
		double current_distance;
		Color3f throughput_over_pdf;
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
			Color3f sigma_t;
			double distance = delta_tracking( scene, ray, its.t, 0, rng, sigma_t);
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
					std::cout << "got scattering event!\n";
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
				Color3f phase_over_pdf = scene->volume->getPhaseFunction()->sample(wi, current_direction, pdf, rng);
				throughput_over_pdf = throughput_over_pdf.cwiseProduct(phase_over_pdf).cwiseProduct(current_vertex.m_sigma_s);
				return true;
			}else
			{
				throw std::runtime_error("surface sampling not implemented");
			}
			return false;
		}



		Color3f nee( RNGd& rng, bool debug )
		{
			// sample light source and transmittance ---
			LightSample ls;
			ls.refP = current_vertex.getPosition();
			Color3f attenuated_light_over_pdf = scene->sample_attenuated_directlight( ls, rng );

			// apply scattering ---
			Color3f phase = scene->volume->getPhaseFunction()->eval( -ls.d, -current_direction ).cwiseProduct(current_vertex.m_sigma_s);

			return attenuated_light_over_pdf.cwiseProduct(phase).cwiseProduct(throughput_over_pdf);
		}
	};

	// intersects the volume bound and starts volume path tracing within the volume
	// will stop once the path exits the volume
	struct VolumePathTracer : public Integrator
	{
		typedef std::shared_ptr<VolumePathTracer> Ptr;

		VolumePathTracer( int maxDepth = std::numeric_limits<int>::max(), bool doSingleScattering = true ):
			Integrator(),
			m_maxDepth(maxDepth),
			m_doSingleScattering(doSingleScattering)
		{
		}


		Color3f trace( TraceInfo& ti, RNGd& rng )const
		{
			Color3f L(0.0f, 0.0f, 0.0f);

			if(ti.debug)
				std::cout <<"VolumePathTracer::trace\n";

			while( ti.depth < m_maxDepth )
			{
				// get next intersection ---
				if( !ti.intersect() )
					// intersection test failed
					return Color3f(0.0, 0.0, 0.0);

				// get next interaction vertex ---
				if( !ti.propagate(rng) )
					// distance sampling failed
					return Color3f(0.0, 0.0, 0.0);

				if(ti.current_vertex.getType() == Vertex::ESurface)
					// we intersected the volume boundary...lets quit
					return L;

				// perform next event estimation ---
				if( ((ti.depth == 0) && m_doSingleScattering ) || (ti.depth > 0) )
					L += ti.nee(rng, ti.debug);

				// handle scattering and find outgoing direction at next vertex ---
				if( !ti.scatter(rng) )
					// directional sampling failed
					return Color3f(0.0, 0.0, 0.0);

				++ ti.depth;
			}

			return L;
		}

		virtual Color3f Li( const Scene* scene, RadianceQuery& rq, RNGd& rng )const override
		{
			rq.transmittance = Color3f(1.0f, 1.0f, 1.0f);

			if( rq.debug )
			{
				std::cout << "VolumePathTracer::Li\n";
				std::cout << "VolumePathTracer::Li pix=" << rq.pixel.toString() << std::endl;
				std::cout << "VolumePathTracer::Li ray=" << rq.ray.toString() << std::endl;
			}


			double mint, maxt;
			if( scene->volume->intersectBound(rq.ray, mint, maxt, rq.debug) )
			{
				// start tracing
				TraceInfo ti;
				ti.depth = 0;
				ti.current_vertex = Vertex();
				ti.current_vertex.setPosition(rq.ray(mint+Epsilon), Color3f(0.0, 0.0, 0.0), Color3f(0.0, 0.0, 0.0));
				ti.current_direction = rq.ray.d;
				ti.scene = scene;
				ti.debug = rq.debug;
				return trace( ti, rng );
			}else
			{
				// no intersection with the medium boundary
			}

			return Color3f(0.0, 0.0, 0.0);
		}

		virtual Color3f sample_transmittance( const Scene* scene, const Ray3d& ray, double maxt, RNGd& rng )const override
		{
			Color3f sigma_t;
			double distance = delta_tracking( scene, ray, maxt, 0, rng, sigma_t );

			if( distance < maxt )
				// since we return 0, we dont need to produce the pdf in the denominator
				return Color3f(0.0, 0.0, 0.0);

			//NB: transmittance sampling pdf cancels out with the transmittance term

			return Color3f(1.0, 1.0, 1.0);
		}

		virtual std::string getId()const override
		{
			if(!m_doSingleScattering)
				return "volpt_ms";
			if((m_maxDepth == 1) && m_doSingleScattering)
				return "volpt_ss";
			return "volpt";
		}

	private:
		bool m_doSingleScattering;
		int m_maxDepth;
	};

	Integrator::Ptr volume_path_tracer( int maxDepth, bool doSingleScattering )
	{
		return std::make_shared<VolumePathTracer>(maxDepth, doSingleScattering);
	}

} // namespace integrators
