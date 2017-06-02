

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

		virtual std::string getId()const override
		{
			return "adrrs";
		}

	private:
		Bitmap::Ptr m_pixel_estimates;
		Field3d::Ptr m_fluence_field;
	};

	Integrator::Ptr adrrs_volume_path_tracer( Bitmap::Ptr pixel_estimates, Field3d::Ptr fluence )
	{
		return std::make_shared<ADRRSVolumePathTracer>(pixel_estimates, fluence);
	}
