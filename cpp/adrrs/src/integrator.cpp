#include <integrator.h>

#include <scene.h>
#include <math/frame.h>

RenderTaskInfoGlobal RenderTaskInfo::g;

namespace integrators
{

	// intersects the volume bound and raymarches through the intersection segment
	struct Raymarcher : public Integrator
	{
		typedef std::shared_ptr<Raymarcher> Ptr;

		Raymarcher():m_stepSize(0.1)
		{
		}

		Raymarcher(double stepSize):m_stepSize(stepSize)
		{
		}

		virtual Color3f Li( const Scene* scene, RadianceQuery& rq, RNGd& rng )const override
		{
			Color3f background(0.0, 0.0, 0.0);
			rq.transmittance = Color3f(1.0f, 1.0f, 1.0f);


			double mint, maxt;
			if( scene->volume->intersectBound(rq.ray, mint, maxt, rq.debug) )
			{
				Color3f f(0.0, 0.0, 0.0);

				int numSteps = int( (maxt-mint)/m_stepSize );


				for( int i=0;i<numSteps;++i )
				{
					double t = mint + i*m_stepSize;
					P3d p = rq.ray(t);
					Color3f extinction = scene->volume->evalExtinction(p);
					Color3f emission = scene->volume->evalAlbedo(p);

					f+=rq.transmittance*Color3f(emission.x()*extinction.x()*m_stepSize,
												emission.y()*extinction.x()*m_stepSize,
												emission.z()*extinction.x()*m_stepSize);

					rq.transmittance.x() *= std::exp(-extinction.x()*m_stepSize);
					rq.transmittance.y() *= std::exp(-extinction.y()*m_stepSize);
					rq.transmittance.z() *= std::exp(-extinction.z()*m_stepSize);
				}
				return f+rq.transmittance*background;
			}else
				// no intersection with the medium boundary
				return background;

			return Color3f(0.0, 0.0, 0.0);
		}

		virtual Color3f sample_transmittance( const Scene* scene, const Ray3d& ray, double maxt, RNGd& rng )const override
		{
			return Color3f(0.0, 0.0, 0.0);
		}

	private:
		double m_stepSize;
	};

	Integrator::Ptr raymarcher(double stepSize)
	{
		return std::make_shared<Raymarcher>(stepSize);
	}


	// intersects the volume bound and raymarches through the intersection segment
	// this raymarcher uses a voxelgrid to lookup precomputed values for fluence
	struct CacheRaymarcher : public Integrator
	{
		typedef std::shared_ptr<CacheRaymarcher> Ptr;

		CacheRaymarcher():m_stepSize(0.1)
		{
		}

		CacheRaymarcher(double stepSize, Cache::Ptr cache):
			m_stepSize(stepSize),
			m_cache(cache)
		{
		}

		virtual Color3f Li( const Scene* scene, RadianceQuery& rq, RNGd& rng )const override
		{
			Color3f background(0.0, 0.0, 0.0);
			rq.transmittance = Color3f(1.0f, 1.0f, 1.0f);


			double mint, maxt;
			if( scene->volume->intersectBound(rq.ray, mint, maxt, rq.debug) )
			{
				Color3f f(0.0, 0.0, 0.0);

				int numSteps = int( (maxt-mint)/m_stepSize );

				double offset = rng.next1D()*m_stepSize;
				double offset_pdf = 1.0/m_stepSize;
				for( int i=0;i<numSteps;++i )
				{
					double t = mint + i*m_stepSize + offset;
					P3d pWS = rq.ray(t);
					V3d d = sampleSphere<double>(rng);
					double d_pdf = sampleSpherePDF();
					Color3f sigma_t = scene->volume->evalExtinction(pWS);
					Color3f sigma_s = sigma_t.cwiseProduct(scene->volume->evalAlbedo(pWS));
					Color3f radiance = m_cache->eval(pWS, d);

					Color3f phase = scene->volume->getPhaseFunction()->eval(d, -rq.ray.d);
					Color3f emission = radiance.cwiseProduct(phase).cwiseProduct(sigma_s);

					f+=rq.transmittance.cwiseProduct(emission)/offset_pdf/d_pdf;

					rq.transmittance.x() *= std::exp(-sigma_t.x()*m_stepSize);
					rq.transmittance.y() *= std::exp(-sigma_t.y()*m_stepSize);
					rq.transmittance.z() *= std::exp(-sigma_t.z()*m_stepSize);
				}
				return f+rq.transmittance*background;
			}else
				// no intersection with the medium boundary
				return background;

			return Color3f(0.0, 0.0, 0.0);
		}

		virtual Color3f sample_transmittance( const Scene* scene, const Ray3d& ray, double maxt, RNGd& rng )const override
		{
			return Color3f(0.0, 0.0, 0.0);
		}


	private:
		double m_stepSize;
		Cache::Ptr m_cache;
	};

	Integrator::Ptr cache_raymarcher(double stepSize, Cache::Ptr cache)
	{
		return std::make_shared<CacheRaymarcher>(stepSize, cache);
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

		VolumePathTracer( int maxDepth = std::numeric_limits<int>::max() ):
			Integrator(),
			m_maxDepth(maxDepth)
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

	private:
		int m_maxDepth;
	};

	Integrator::Ptr volume_path_tracer( int maxDepth )
	{
		return std::make_shared<VolumePathTracer>(maxDepth);
	}



	// intersects the volume bound and starts volume path tracing within the volume
	// will stop once the path exits the volume
	struct ADRRSVolumePathTracer : public Integrator
	{
		typedef std::shared_ptr<VolumePathTracer> Ptr;

		ADRRSVolumePathTracer( Bitmap::Ptr pixel_estimates, Field3d::Ptr fluence ):
			Integrator(),
			m_pixel_estimates(pixel_estimates),
			m_fluence_field(fluence)
		{
		}


		Color3f trace( TraceInfo& ti, RNGd& rng )const
		{
			Color3f L(0.0f, 0.0f, 0.0f);

			if(ti.debug)
				std::cout <<"VolumePathTracer::trace\n";

			//while(true)
			while( ti.depth < 2 )
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
				if(ti.depth == 1)
					L += ti.nee(rng, ti.debug);

				++ti.depth;

				// russian roulette ---
				{

				}

				// splitting ---
				{
					int n = 50;
					ti.throughput_over_pdf *= 1.0/double(n);
					for( int i=0;i<n;++i )
					{
						TraceInfo split = ti; // make a copy
						if( split.scatter(rng) )
							L += trace(split, rng);
					}
					//TODO: what we could do is to spawn n-1 traces and continue with the last one
					return L;
				}



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

	private:
		Bitmap::Ptr m_pixel_estimates;
		Field3d::Ptr m_fluence_field;
	};

	Integrator::Ptr adrrs_volume_path_tracer( Bitmap::Ptr pixel_estimates, Field3d::Ptr fluence )
	{
		return std::make_shared<ADRRSVolumePathTracer>(pixel_estimates, fluence);
	}


	struct Dummy :public Integrator
	{
		virtual Color3f Li( const Scene* scene, RadianceQuery& rq, RNGd& rng )const override
		{
			return Color3f(1.0);
		}

		virtual Color3f sample_transmittance( const Scene* scene, const Ray3d& ray, double maxt, RNGd& rng )const override
		{
			return Color3f(1.0);
		}
	};

	Integrator::Ptr dummy()
	{
		return std::make_shared<Dummy>();
	}

} // namespace integrator


