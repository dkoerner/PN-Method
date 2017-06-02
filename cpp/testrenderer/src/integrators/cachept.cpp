
	// intersects the volume bound and uses delta_tracking to go through the volume
	// cache is used to lookup indirect illumination
	// NB: originally this was a raymarcher but the introduced bias would cause a
	// noticeable mismatch to the groundtruth which I want to avoid
	struct CachePathTracer : public Integrator
	{
		typedef std::shared_ptr<CachePathTracer> Ptr;

		CachePathTracer(Cache::Ptr cache):
			m_cache(cache)
		{
		}


		/*
		Color3f nee( const Scene* scene, const P3d& pWS, const V3d& current_direction, Color3f sigma_s_at_p, RNGd& rng, bool debug = false )const
		{
			// sample light source and transmittance ---
			LightSample ls;
			ls.refP = pWS;
			Color3f attenuated_light_over_pdf = scene->sample_attenuated_directlight( ls, rng, debug );
			if(debug)
				std::cout << "\tattenuated_light_over_pdf=" << attenuated_light_over_pdf.toString() << std::endl;

			// apply scattering ---
			Color3f phase = scene->volume->getPhaseFunction()->eval( -ls.d, -current_direction ).cwiseProduct(sigma_s_at_p);

			return attenuated_light_over_pdf.cwiseProduct(phase);
		}
		*/


		virtual Color3f Li( const Scene* scene, RadianceQuery& rq, RNGd& rng )const override
		{
			Color3f background(0.0, 0.0, 0.0);
			rq.transmittance = Color3f(1.0f, 1.0f, 1.0f);


			double mint, maxt;
			if( scene->volume->intersectBound(rq.ray, mint, maxt, rq.debug) )
			{
				// find scattering event ---
				Ray3d rayWS( rq.ray( mint+Epsilon ), rq.ray.d );
				double max_distance = maxt-mint;
				Color3f sigma_t;
				double t = delta_tracking(scene, rayWS, max_distance, 0, rng, sigma_t );


				if(t<max_distance)
				{
					// scattering event found
					P3d pWS = rayWS(t);
					Color3f albedo = scene->volume->evalAlbedo(pWS);
					Color3f sigma_s = sigma_t*albedo;

					// NB: transmittance cancels out with delta_tracking pdf
					Color3f throughput_over_pdf = Color3f(1.0f/sigma_t.x(), 1.0f/sigma_t.y(), 1.0f/sigma_t.z());
					Color3f scatter_attenuation = sigma_s*INV_FOURPI; // assuming isotropic scattering here...


					// TODO: direct illumination
					Color3f direct_light_over_pdf(0.0f, 0.0f, 0.0f);

					// use cache for indirect illum
					Color3f indirect_light_over_pdf;
					{
						// sample random direction
						V3d d = sampleSphere<double>(rng);
						double d_pdf = sampleSpherePDF();
						indirect_light_over_pdf = throughput_over_pdf*scatter_attenuation*m_cache->eval( pWS, d )/d_pdf;
					}

					return direct_light_over_pdf+indirect_light_over_pdf;
				}else
					// no scattering event, sample leaves the volume
					// NB: transmittance terms cancels with delta_tracking transmittance
					return background;
			}else
				// no intersection with the medium boundary
				return background;
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
			return "cacherm_"+m_cache->getId();
		}


	private:
		Cache::Ptr m_cache;
	};

	Integrator::Ptr volume_path_tracer_cached(Cache::Ptr cache)
	{
		return std::make_shared<CachePathTracer>(cache);
	}

