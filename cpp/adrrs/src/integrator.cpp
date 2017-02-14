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

	private:
		double m_stepSize;
	};

	Integrator::Ptr raymarcher(double stepSize)
	{
		return std::make_shared<Raymarcher>(stepSize);
	}

	struct Intersection
	{
		void reset()
		{
			t = std::numeric_limits<double>::infinity();
		}

		P3d p;
		double t;
	};

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

		bool scene_intersect( const Ray3d& ray, Intersection& its )
		{
			its.reset();

			double mint, maxt;
			// TODO: take into account whether vertex is surface or volume (Epsilon)...
			bool result = scene->volume->intersectBound(ray, mint, maxt);
			if( result )
			{
				its.t = maxt;
				its.p = ray(its.t);
				return true;
			}
			return false;
		}

		bool intersect()
		{
			return scene_intersect(Ray3d( current_vertex.getPosition(), current_direction ), its);
		}

		// returns sigma_t at sampled position (is invalid when we exceeded maxt)
		double delta_tracking( const Ray3d& ray, double maxt, int component, RNGd& rng, Color3f& sigma_t )
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

		bool propagate(RNGd& rng)
		{
			Ray3d ray( current_vertex.getPosition(), current_direction );
			Color3f sigma_t;
			double distance = delta_tracking(ray, its.t, 0, rng, sigma_t);
			if( distance < its.t )
			{
				// volume scattering event ---
				current_vertex.setPosition( ray(distance), sigma_t, scene->volume->evalAlbedo(its.p));

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
			// sample light source ---
			LightSample ls;
			ls.refP = current_vertex.getPosition();
			Color3f light_over_pdf = scene->light->sample(ls);

			double maxt =ls.distance;

			// find next intersection point ---
			{
				Intersection its;
				if( scene_intersect( Ray3d( current_vertex.getPosition(), -ls.d ), its ) )
				{
					// there is something between current position and lightsource
					//TODO: we currently assume its just the volume boundary and we are save to ignore it
					maxt = std::min( its.t, ls.distance );
				}else
				{
					// currently we expect to be within a volume
					return Color3f(0.0f, 0.0f, 0.0f);
				}
			}


			// sample transmittance ---
			Color3f sigma_t;
			double distance = delta_tracking( Ray3d( current_vertex.getPosition(), -ls.d ), maxt, 0, rng, sigma_t );

			if( distance < maxt )
				// since we return 0, we dont need to produce the pdf in the denominator
				return Color3f(0.0, 0.0, 0.0);

			//NB: transmittance sampling pdf cancels out with the transmittance term

			//NB: the geometry term between light source position and current vertex is already been dealt with during light sampling

			// apply scattering ---
			Color3f phase = scene->volume->getPhaseFunction()->eval( -ls.d, -current_direction ).cwiseProduct(current_vertex.m_sigma_s);

			return light_over_pdf.cwiseProduct(phase);
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

	private:
		int m_maxDepth;
	};

	Integrator::Ptr volume_path_tracer( int maxDepth )
	{
		return std::make_shared<VolumePathTracer>(maxDepth);
	}


/*
	// =============================================================================
	Color3f raymarching_sample(RenderTaskInfo& tinf, const Ray3d& rayWS, double mint, double maxt, const V2i& pix, Color3f& T, bool debug, int numSteps)
	{
		Color3f f(0.0, 0.0, 0.0);

		//int numSamples = 1000;
		double dt = (maxt-mint)/double(numSteps);

		T = Color3f(1.0f, 1.0f, 1.0f);
		for( int i=0;i<numSteps;++i )
		{
			double t = mint + i*dt;
			P3d p = rayWS(t);
			V3d extinction = tinf.g.volume->evalExtinction(p);
			V3d emission = tinf.g.volume->evalAlbedo(p);

			f+=T*Color3f(emission.x()*extinction.x()*dt,
						 emission.y()*extinction.x()*dt,
						 emission.z()*extinction.x()*dt);

			T.x() *= std::exp(-extinction.x()*dt);
			T.y() *= std::exp(-extinction.y()*dt);
			T.z() *= std::exp(-extinction.z()*dt);
		}
		return f+T*tinf.g.background;
	}

	Integrator raymarching( int numSteps )
	{
		using namespace std::placeholders;

		Integrator i;
		i.id = "groundtruth";
		i.info = "numSteps="+toString(numSteps);
		i.sample_function = std::bind(raymarching_sample, _1, _2, _3, _4, _5, _6, _7, numSteps);
		return i;
	}


	// =============================================================================

	Color3f randomoffset_sample(RenderTaskInfo& tinf, const Ray3d& rayWS, double mint, double maxt, const V2i& xy, Color3f& T, bool debug, int numSteps)
	{
		Color3f f(0.0, 0.0, 0.0);

		double dt = (maxt-mint)/double(numSteps);

		double offset = tinf.rng.next1D();
		//double offset = 0.0;
		for( int i=0;i<numSteps;++i )
		{
			double ti = mint + i*dt;
			double t = ti+offset*dt;

			P3d p = rayWS(t);
			V3d extinction = tinf.g.volume->evalExtinction(p);
			V3d albedo = tinf.g.volume->evalAlbedo(p);
			double pdf = 1.0/dt;

			Color3f contribution = Color3f(float(albedo.x()*extinction.x()),
										   float(albedo.y()*extinction.y()),
										   float(albedo.z()*extinction.z()))*dt;

			f += T*contribution;

			T.x() *= std::exp(-extinction.x()/pdf);
			T.y() *= std::exp(-extinction.y()/pdf);
			T.z() *= std::exp(-extinction.z()/pdf);
		}


		return f+T*tinf.g.background;
	}

	Integrator randomoffset( int numSteps )
	{
		using namespace std::placeholders;
		Integrator i;
		i.id = "randomoffset";
		i.info = "numSteps="+toString(numSteps);
		i.sample_function = std::bind(randomoffset_sample, _1, _2, _3, _4, _5, _6, _7, numSteps);
		return i;
	}


	// =============================================================================
	double deltatracking_sampleDistance( const Ray3d& rayWS, double mint, double maxt, const Volume* volume, RNGd& rng, double t_pdf, bool debug = false )
	{
		// since we render 3 color channels at the same time, we need to decide which channel to use for
		// importance sampling
		int component = 0;
		//double extinction_max = volume->extinction_max[component];
		double extinction_max = volume->getMaxExtinction()[component];


		double t = mint;
		while(true)
		{
			double step = -log( 1.0-rng.next1D() )/extinction_max;
			t += step;

			if(t>= maxt)
				break;

			double extinction = volume->evalExtinction(rayWS(t))[component];

			// russian roulette
			if(rng.next1D()<extinction/extinction_max)
				break;
		}

		return t;
	}

	Color3f deltatracking_sample(RenderTaskInfo& ti, const Ray3d& rayWS, double mint, double maxt, const V2i& pix, Color3f& T, bool debug = false)
	{
		Color3f f(0.0, 0.0, 0.0);

		// single scattering
		double t_pdf;
		double t = deltatracking_sampleDistance( rayWS, mint, maxt, ti.g.volume, ti.rng, t_pdf, debug );

		// if scattering within volume
		if(t<maxt)
		{
			// note: distance sampling pdf cancels out with transmittance
			//		what remains is 1/extinction at p
			P3d p = rayWS(t);

			V3d albedo = ti.g.volume->evalAlbedo(p);
			f = Color3f(albedo.x(), albedo.y(), albedo.z());
			T = Color3f(0.0f, 0.0f, 0.0f);
		}else
		{
			// no scattering within the medium
			// note: distance sampling pdf cancels out with transmittance
			//f = background;
			f = Color3f(1.0, 1.0, 1.0)*ti.g.background;
			T = Color3f(1.0f, 1.0f, 1.0f);
		}
		return f;
	}

	Integrator deltatracking()
	{
		//using namespace std::placeholders;
		Integrator i;
		i.id = "deltatracking";
		//i.sample_function = std::bind(randomoffset_sample, _1, _2, _3, _4, _5, _6, _7, numSteps);
		i.sample_function = deltatracking_sample;
		return i;
	}
*/

}
