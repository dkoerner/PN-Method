#pragma once
#include <map>
#include <integrators/simplept.h>

#include <scene.h>
#include <math/frame.h>

#include <volume.h>

#include <PNSolution.h>











struct TraceInfoPNIS
{
	Vertex current_vertex;
	V3d current_direction;
	double current_distance;
	V3d throughput_over_pdf;
	Intersection its;
	int depth;
	const Scene* scene;
	bool debug;
	const PNSolution* pns;

	TraceInfoPNIS():
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

			double phase_over_pdf;


			// we use one-sample-estimator to decide between phase function and PNS-sampling
			double prob = 0.5;
			if( rng.next1D() < prob )
			{
				// sample according to PN-Solution
				double pdf_shsampling;
				current_direction = -pns->sample( current_vertex.getPosition(), pdf_shsampling, P2d(rng.next1D(), rng.next1D()) );

				if(pdf_shsampling == 0.0)
					return false;

				double pdf_phase = INV_FOURPI;

				double pdf = prob*pdf_shsampling + (1.0-prob)*pdf_phase;

				phase_over_pdf = scene->volume->evalPhase(current_vertex.getPosition(), wi, current_direction)/pdf;

				if( debug )
				{
					dbg()["phase_sampling"] = 0;
					dbg()["phase_pdf"] = pdf;
					dbg()["pdf_shsampling"] = pdf_shsampling;
				}
			}else
			{
				// sample according to phase function
				double pdf_phase;
				scene->volume->samplePhase(current_vertex.getPosition(), wi, current_direction, pdf_phase, rng);

				double pdf_shsampling = pns->pdf(current_vertex.getPosition(), -current_direction);


				double pdf = prob*pdf_shsampling + (1.0-prob)*pdf_phase;
				//double pdf = pdf_phase;


				phase_over_pdf = scene->volume->evalPhase(current_vertex.getPosition(), wi, current_direction)/pdf;

				if( debug )
				{
					dbg()["phase_sampling"] = 1;
					dbg()["phase_pdf"] = pdf;
					dbg()["pdf_shsampling"] = pdf_shsampling;

				}
			}

			if( debug )
			{
				dbg()["phase_over_pdf"] = phase_over_pdf;
				dbg()["pWS_x"] = current_vertex.getPosition().x();
				dbg()["pWS_y"] = current_vertex.getPosition().y();
				dbg()["pWS_z"] = current_vertex.getPosition().z();
				dbg()["dir_x"] = current_direction.x();
				dbg()["dir_y"] = current_direction.y();
				dbg()["dir_z"] = current_direction.z();
			}

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


	std::vector<std::map<std::string, double>> *debug_info;
	std::map<std::string, double>& dbg()const
	{
		return (*debug_info).back();
	}


};







// intersects the volume bound and starts volume path tracing within the volume
// will stop once the path exits the volume
struct PNISPT : public Integrator
{
	typedef std::shared_ptr<PNISPT> Ptr;

	PNISPT( PNSolution::Ptr pns,
			bool doSingleScattering = true,
			int maxDepth = std::numeric_limits<int>::max() ):
		Integrator(),
		m_pns(pns),
		m_maxDepth(maxDepth),
		m_doSingleScattering(doSingleScattering)
	{
	}

	V3d trace( TraceInfoPNIS& ti, RNGd& rng )const
	{
		V3d L(0.0f, 0.0f, 0.0f);

		if(ti.debug)
			std::cout <<"VolumePathTracer::trace\n";

		//L += ti.nee(rng, ti.debug);
		//return L;

		while( ti.depth < m_maxDepth )
		{
			if(ti.debug)
			{
				std::cout <<"VolumePathTracer::trace depth=" << ti.depth << std::endl;
				dbgAdd();
				dbg()["depth"] = ti.depth;
				dbg()["throughput_over_pdf"] = ti.throughput_over_pdf.x();
				dbg()["L"] = L.x();
			}
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
				//L += nee(ti, rng);
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





	// debugging stuff -------------
	mutable std::vector<std::map<std::string, double>> debug_info;

	void dbgAdd()const
	{
		debug_info.push_back(std::map<std::string, double>());
	}

	std::map<std::string, double>& dbg()const
	{
		return debug_info.back();
	}

	Eigen::VectorXd dbgGet( const std::string& id )
	{
		int numElements = debug_info.size();
		Eigen::VectorXd vec(numElements);
		for( int i=0;i<numElements;++i )
		{
			std::map<std::string, double>& map = debug_info[i];
			if( map.find(id) != map.end() )
				vec.coeffRef(i) = map[id];
			else
			{
				vec.coeffRef(i) = std::numeric_limits<double>::quiet_NaN();
			}
		}
		return vec;
	}



private:
	PNSolution::Ptr m_pns;
	bool m_doSingleScattering;
	int m_maxDepth;
};




/*
// intersects the volume bound and starts volume path tracing within the volume
// will stop once the path exits the volume
struct PNISPT : public Integrator
{
	typedef std::shared_ptr<PNISPT> Ptr;

	PNISPT( PNSolution::Ptr pns,
			bool doSingleScattering = true,
			int maxDepth = std::numeric_limits<int>::max() ):
	    Integrator(),
	    m_pns(pns),
	    m_maxDepth(maxDepth),
	    m_doSingleScattering(doSingleScattering)
	{
	}


	//V3d trace( TraceInfoPNIS& ti, RNGd& rng )const
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
			{
				std::cout <<"VolumePathTracer::trace depth=" << ti.depth << std::endl;
				dbgAdd();
				dbg()["depth"] = ti.depth;
				dbg()["throughput_over_pdf"] = ti.throughput_over_pdf.x();
				dbg()["L"] = L.x();
			}

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
				//L += ti.nee(rng, ti.debug);
				L += nee( ti, rng);

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
		ss << "PNISPT singleScattering=" << m_doSingleScattering << " maxDepth=" << m_maxDepth << std::endl;
		return ss.str();
	}






	// debugging stuff -------------
	mutable std::vector<std::map<std::string, double>> debug_info;

	void dbgAdd()const
	{
		debug_info.push_back(std::map<std::string, double>());
	}

	std::map<std::string, double>& dbg()const
	{
		return debug_info.back();
	}

	Eigen::VectorXd dbgGet( const std::string& id )
	{
		int numElements = debug_info.size();
		Eigen::VectorXd vec(numElements);
		for( int i=0;i<numElements;++i )
		{
			std::map<std::string, double>& map = debug_info[i];
			if( map.find(id) != map.end() )
				vec.coeffRef(i) = map[id];
			else
			{
				vec.coeffRef(i) = std::numeric_limits<double>::quiet_NaN();
			}
		}
		return vec;
	}

private:
	PNSolution::Ptr m_pns;
	bool m_doSingleScattering;
	int m_maxDepth;
};
*/
