#pragma once
#include <math/color.h>


#include <volume.h>
#include <light.h>
#include <camera.h>
//#include <util/bitmap.h>

#include <math/common.h>
#include <math/rng.h>

//#include <cache.h>



struct Scene;



struct RadianceQuery
{
	RadianceQuery():
		debug(false),
	    pixel(-1, -1),
	    volume(0)
	{
	}

	Ray3d ray;
	V2i pixel;
	bool debug;
	const Volume* volume;


	//V3d transmittance;
};


struct Integrator
{
	typedef std::shared_ptr<Integrator> Ptr;

	virtual V3d Li( const Scene* scene, RadianceQuery& rq, RNGd& rng )const=0;
	virtual V3d sample_transmittance( const Scene* scene, const Ray3d& ray, double maxt, RNGd& rng )const=0;
	//virtual std::string getId()const=0;
	virtual std::string toString()const=0;
};


/*
	struct Dummy :public Integrator
	{
		virtual Color3f Li( const RenderScene* scene, RadianceQuery& rq, RNGd& rng )const override
		{
			return Color3f(1.0);
		}

		virtual Color3f sample_transmittance( const RenderScene* scene, const Ray3d& ray, double maxt, RNGd& rng )const override
		{
			return Color3f(1.0);
		}

		virtual std::string getId()const override
		{
			return "dummy";
		}
	};

	Integrator::Ptr dummy()
	{
		return std::make_shared<Dummy>();
	}


  */
namespace integrators
{
/*
	//Integrator::Ptr dummy();
	Integrator::Ptr raymarcher(double stepsize = 0.1);
	Integrator::Ptr volume_path_tracer_cached(Cache::Ptr cache );
	Integrator::Ptr volume_path_tracer(int maxDepth = std::numeric_limits<int>::max(), bool doSingleScattering = true);
	Integrator::Ptr adrrs_volume_path_tracer( Bitmap::Ptr pixel_estimates, Field3d::Ptr fluence );
*/
}

