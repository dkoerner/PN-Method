#pragma once

#include <cuda/CudaVolume.h>
#include <cuda/CudaLight.h>

/*
// the version below takes volume as density grid and extinction cross section
// kcs = extinction crosssection (density*kcs = kt = extinction coefficient)
template<typename T>
__device__ T pt_T_sample( const cumath::RNG &rng, const CudaField<T> *density, T density_max, T kcs, const cumath::Ray3<T> &ray )
{
	cumath::BoundingBox3<T> bound = density->m_bound;

	T t0, t1;
	if( !cumath::intersectionRayBox<T>( ray, bound, t0, t1 ) )
		// no intersection
		return T(1.0);

	T l = t0;
	T kt_max = density_max*kcs;
	while( true )
	{
		T dl = -log(T(rng.randomFloat()))/kt_max;
		l += dl;
		if( l > t1 )
			return T(1.0);
		cumath::Vec3<T> vsP = cumath::transform( ray.getPosition(l), density->m_worldToVoxel );
		T k_t = density->eval(vsP)*kcs;
		if( T(rng.randomFloat()) < k_t/kt_max )
			return T(0.0);
	}


	return T(1.0);
}



// MonteCarlo single+multiple scattering for given ray and voxelgrid ====
// the version below takes volume as density grid and extinction cross section
// kcs = extinction crosssection (density*kcs = kt = extinction coefficient)
template<typename T>
__device__ T pt_ms( const cumath::RNG &rng, const CudaField<T> *density, T density_max, T kcs, T albedo, cumath::Vec3<T> lightDir, T lightIntensity, cumath::Ray3<T> startRay, int maxDepth )
{
	int rrDepth = 3; // indicates when we will start with russian roulette exit tests
	cumath::BoundingBox3<T> bound = density->m_bound;
	T k_extinct_max = density_max*kcs;

	cumath::Ray3<T> ray = startRay;
	T throughput = 1.0f;
	T pathpdf = 1.0f;
	int depth = 0;
	T L = 0.0f;

	while( depth < maxDepth )
	{
		// sample ray (propagation event) ---
		T t0, t1;
		if( !cumath::intersectionRayBox<T>( ray, bound, t0, t1 ) )
			// no intersection
			break;

		T l = t0;
		T k_t = 0.0f;
		bool accepted = false;
		cumath::Vec3<T> p;

		while( true )
		{
			T dl = -log(T(rng.randomFloat()))/k_extinct_max;
			l += dl;
			if( l > t1 )
				break;
			p = ray.getPosition(l);

			cumath::Vec3<T> vsP = transform(p, density->m_worldToVoxel );
			k_t = density->eval(vsP)*kcs;
			if( T(rng.randomFloat()) < k_t/k_extinct_max )
			{
				accepted = true;
				break;
			}
		}


		// if scatter event within slab ---
		if(accepted)
		{
			pathpdf *= k_t; // *T

			// direct lighting contribution at point p
			cumath::Ray3<T> lightRay;
			lightRay.o = p;
			lightRay.d = lightDir;
			if( (pathpdf>(T)0.0) && (pt_T_sample<T>( rng, density, density_max, kcs, lightRay ) > (T)0.0) )
			{
				// phase function
				T phase = ((T)1.0/((T)4.0*MATH_PI));
				T L_direct = lightIntensity*phase*albedo*k_t*throughput;
				L += L_direct / pathpdf;
			}

			// scatter event (isotropic phase function)
			T phase = ((T)1.0/((T)4.0*MATH_PI));
			T wnPDF = T(0.0);
			cumath::Vec3<T> wn = cumath::sample_sphere( rng, &wnPDF );
			throughput *= phase*albedo*k_t;
			pathpdf *= wnPDF;

			// proceed
			ray.o = p;
			ray.d = wn;

			++depth;
		}else
			break;

		if( depth > rrDepth )
		{
			T q = T(0.95);
			if( rng.randomFloat() > q )
				break;
			pathpdf *= q;
		}
	};


	return L;
}



// version below takes volume as extinction coefficient grid directly
template<typename T>
__device__ T pt_T_sample( const cumath::RNG &rng, const CudaField<T> *kt, T kt_max, const cumath::Ray3<T> &ray )
{
	cumath::BoundingBox3<T> bound = kt->m_bound;

	T t0, t1;
	if( !cumath::intersectionRayBox<T>( ray, bound, t0, t1 ) )
		// no intersection
		return T(1.0);

	T l = t0;
	while( true )
	{
		T dl = -log(T(rng.randomFloat()))/kt_max;
		l += dl;
		if( l > t1 )
			return T(1.0);
		cumath::Vec3<T> vsP = cumath::transform( ray.getPosition(l), kt->m_worldToVoxel );
		T k_t = kt->eval(vsP);
		if( T(rng.randomFloat()) < k_t/kt_max )
			return T(0.0);
	}

	return T(1.0);
}

*/

struct Framed
{

};

struct Intersection
{
	DEVICE void reset()
	{
		t = CUDART_INF_F;
	}

	cumath::V3f p;
	float t;
};

struct Vertex
{
	enum EVertexType
	{
		EInvalid = 0,
		ESurface = 1,
		EVolume = 2
	};

	DEVICE void setPosition( const cumath::V3f& position, cumath::V3f sigma_t, cumath::V3f albedo )
	{
		m_p = position;
		m_type = EVolume;
		m_sigma_t = sigma_t;
		m_albedo = albedo;
		m_sigma_s = sigma_t*albedo;
	}
	DEVICE void setPosition( const cumath::V3f& position, const Framed& frame)
	{
		m_p = position;
		m_type = ESurface;
	}

	DEVICE const cumath::V3f& getPosition()const
	{
		return m_p;
	}

	DEVICE EVertexType getType()const
	{
		return m_type;
	}

//		private:
	cumath::V3f m_p; // position in worldspace
	EVertexType m_type; // volume or surface scattering point
	cumath::V3f m_sigma_t;
	cumath::V3f m_sigma_s;
	cumath::V3f m_albedo;
};


// returns sigma_t at sampled position (is invalid when we exceeded maxt)
DEVICE float delta_tracking( const CudaVolume* volume, const cumath::Ray3f& ray, double maxt, int component, cumath::RNG& rng, cumath::V3f& sigma_t )
{
	float sigma_t_max = volume->getMaxExtinction()[component];

	float t = 0.0;
	while(true)
	{
		float step = -log( 1.0-rng.randomFloat() )/sigma_t_max;
		t += step;

		if(t>= maxt)
			break;

		sigma_t = volume->evalExtinction(ray.getPosition(t));

		// russian roulette
		if(rng.randomFloat()<sigma_t[component]/sigma_t_max)
			break;
	}

	return t;
}


struct TraceInfo
{
	bool debug;
	int depth;
	int maxDepth;
	Vertex current_vertex;
	cumath::V3f current_direction;
	Intersection its; // next surface intersection
	double current_distance;
	cumath::V3f throughput_over_pdf;
	CudaVolume* volume;
	CudaLight* light;


	DEVICE TraceInfo():
		throughput_over_pdf(1.0f, 1.0f, 1.0f),
		debug(false)
	{
	}


	//NB: the intersection routine assumes that the ray is within the bounding box...
	DEVICE bool intersect()
	{
		its.reset();

		float tmin, tmax;
		cumath::Ray3f ray(current_vertex.getPosition(), current_direction, 0.0, 100.0);
		bool result = volume->intersectBound(ray, tmin, tmax, debug);
		if(result)
		{
			its.t = tmax;
			its.p = ray.getPosition(its.t);
		}
		if(debug)
		{
			printf("TraceInfo::intersect: result=%i tmin=%f, tmax=%f\n", int(result), tmin, tmax);
			if(result==true)
				printf("result true!\n");
			else
				printf("result false!\n");
		}
		return result;
	}

	DEVICE bool propagate( cumath::RNG& rng )
	{
		cumath::Ray3f ray( current_vertex.getPosition(), current_direction, 0.0, 100.0 );
		cumath::V3f sigma_t;
		double distance = delta_tracking( volume, ray, its.t, 0, rng, sigma_t);
		if( distance < its.t )
		{
			cumath::V3f pWS = ray.getPosition(distance);
			// volume scattering event ---
			current_vertex.setPosition( pWS, sigma_t, volume->evalAlbedo(pWS));

			// apply pdf and volume attenuation
			throughput_over_pdf[0] *= 1.0/sigma_t.x;
			throughput_over_pdf[1] *= 1.0/sigma_t.y;
			throughput_over_pdf[2] *= 1.0/sigma_t.z;

			if( debug )
			{
				printf("got scattering event\n");
			}
		}else
		{
			//printf("got surface scattering event\n");
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



	DEVICE cumath::V3f integrator_sample_transmittance( const cumath::Ray3f& ray, float maxt, cumath::RNG& rng )
	{
		cumath::V3f sigma_t;
		float distance = delta_tracking( volume, ray, maxt, 0, rng, sigma_t );

		if( distance < maxt )
			// since we return 0, we dont need to produce the pdf in the denominator
			return cumath::V3f(0.0, 0.0, 0.0);

		//NB: transmittance sampling pdf cancels out with the transmittance term

		return cumath::V3f(1.0, 1.0, 1.0);
	}

	DEVICE bool intersect( const cumath::Ray3f& ray, Intersection& inters )
	{
		inters.reset();

		float mint, maxt;
		// TODO: take into account whether vertex is surface or volume (Epsilon)...
		bool result = volume->intersectBound(ray, mint, maxt);
		if( result )
		{
			//TODO: here we assume that mint is behind are exactly at the ray origin
			// if the ray is outside the box, then mint would be the correct choice
			inters.t = maxt;
			inters.p = ray.getPosition(inters.t);
			return true;
		}
		return false;
	}

	DEVICE cumath::V3f scene_sample_attenuated_directlight( LightSample& ls, cumath::RNG& rng )
	{
		// sample lightsource (including geometry term)
		cumath::V3f light_over_pdf = light->sample(ls, rng);

		float maxt = ls.distance;

		cumath::Ray3f ray_light( ls.refP, -ls.d, 0.0, 100.0 );

		// find next intersection point ---
		{
			Intersection its_light;
			if( intersect( ray_light, its_light ) )
			{
				// there is something between current position and lightsource
				//TODO: we currently assume its just the volume boundary and we are save to ignore it
				maxt = min( its_light.t, ls.distance );
			}else
			{
				// currently we expect to be within a volume
				return cumath::V3f(0.0f, 0.0f, 0.0f);
			}
		}

		//NB: the geometry term between light source position and current vertex is already been dealt with during light sampling

		// sample transmittance ---
		cumath::V3f transmittance_over_pdf = integrator_sample_transmittance( ray_light, maxt, rng );

		return light_over_pdf*transmittance_over_pdf;
	}

	DEVICE cumath::V3f nee(cumath::RNG& rng)
	{
		// sample light source and transmittance ---
		LightSample ls;
		ls.refP = current_vertex.getPosition();
		cumath::V3f attenuated_light_over_pdf = scene_sample_attenuated_directlight( ls, rng );

		// apply scattering ---
		cumath::V3f phase = volume->evalPhaseFunction(-ls.d, -current_direction )*current_vertex.m_sigma_s;

		return attenuated_light_over_pdf*phase*throughput_over_pdf;
	}

	DEVICE bool scatter(cumath::RNG& rng)
	{
		if( current_vertex.getType() == Vertex::EVolume )
		{
			cumath::V3f wi = current_direction;
			float pdf;
			cumath::V3f phase_over_pdf = volume->samplePhaseFunction(wi, current_direction, pdf, rng);
			throughput_over_pdf = throughput_over_pdf*phase_over_pdf*current_vertex.m_sigma_s;
			return true;
		}else
		{
			//throw std::runtime_error("surface sampling not implemented");
		}
		return false;
	}


};

__device__ cumath::V3f trace( TraceInfo& ti, cumath::RNG &rng )
{
	cumath::V3f L(0.0f, 0.0f, 0.0f);

	//if(ti.debug)
	//	std::cout <<"VolumePathTracer::trace\n";

	while( ti.depth < ti.maxDepth )
	{
		// get next intersection ---
		if( !ti.intersect() )
			// intersection test failed
			return cumath::V3f(0.0, 1.0, 0.0);

		// get next interaction vertex ---
		if( !ti.propagate(rng) )
			// distance sampling failed
			return cumath::V3f(0.0, 0.0, 0.0);


		if(ti.current_vertex.getType() == Vertex::ESurface)
			// we intersected the volume boundary...lets quit
			return L;


		// perform next event estimation ---
		L += ti.nee(rng);

		// handle scattering and find outgoing direction at next vertex ---
		if( !ti.scatter(rng) )
			// directional sampling failed
			return cumath::V3f(0.0, 0.0, 0.0);

		++ti.depth;
	}

	return L;
}


/*
// MonteCarlo single+multiple scattering for given ray and voxelgrid ====
template<typename T>
__device__ T pt_ms( const cumath::RNG &rng, const CudaField<T> *kt, T kt_max, T albedo, cumath::Vec3<T> lightDir, T lightIntensity, cumath::Ray3<T> startRay, int maxDepth, T pathpdf = T(1.0) )
{
	int rrDepth = 3; // indicates when we will start with russian roulette exit tests
	cumath::BoundingBox3<T> bound = kt->m_bound;
	T k_extinct_max = kt_max;

	cumath::Ray3<T> ray = startRay;
	T throughput = 1.0f;
	//T pathpdf = 1.0f;
	int depth = 0;
	T L = 0.0f;

	while( depth < maxDepth )
	{
		// sample ray (propagation event) ---
		T t0, t1;
		if( !cumath::intersectionRayBox<T>( ray, bound, t0, t1 ) )
			// no intersection
			break;

		T l = t0;
		T k_t = 0.0f;
		bool accepted = false;
		cumath::Vec3<T> p;

		while( true )
		{
			T dl = -log(T(rng.randomFloat()))/k_extinct_max;
			l += dl;
			if( l > t1 )
				break;
			p = ray.getPosition(l);

			cumath::Vec3<T> vsP = transform(p, kt->m_worldToVoxel );
			k_t = kt->eval(vsP);
			if( T(rng.randomFloat()) < k_t/k_extinct_max )
			{
				accepted = true;
				break;
			}
		}


		// if scatter event within slab ---
		if(accepted)
		{
			pathpdf *= k_t; // *T

			// direct lighting contribution at point p
			cumath::Ray3<T> lightRay;
			lightRay.o = p;
			lightRay.d = lightDir;
			if( (pathpdf>(T)0.0) && (pt_T_sample<T>( rng, kt, kt_max, lightRay ) > (T)0.0) )
			{
				// phase function
				T phase = ((T)1.0/((T)4.0*MATH_PI));
				T L_direct = lightIntensity*phase*albedo*k_t*throughput;
				L += L_direct / pathpdf;
			}


			// scatter event (isotropic phase function)
			T phase = ((T)1.0/((T)4.0*MATH_PI));
			T wnPDF = T(0.0);
			cumath::Vec3<T> wn = cumath::sample_sphere( rng, &wnPDF );
			throughput *= phase*albedo*k_t;
			pathpdf *= wnPDF;

			// proceed
			ray.o = p;
			ray.d = wn;
			++depth;
		}else
			break;

		if( depth > rrDepth )
		{
			T q = T(0.95);
			if( rng.randomFloat() > q )
				break;
			pathpdf *= q;
		}
	};

	return L;
}


*/
