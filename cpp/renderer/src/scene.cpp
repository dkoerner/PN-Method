#include <scene.h>
#include <volume.h>
#include <integrator.h>

bool Scene::intersect( const Ray3d& ray, Intersection& its, bool debug )const
{
	its.reset();

	if(debug)
	{
		std::cout << "Scene::intersect ray=" << ray.toString() << std::endl;
	}

	double mint, maxt;
	// TODO: take into account whether vertex is surface or volume (Epsilon)...
	bool result = volume->intersectBound(ray, mint, maxt, debug);
	if( result )
	{
		if(debug)
		{
			std::cout << "Scene::intersect: got intersection at t=" << maxt << std::endl;
		}
		//TODO: here we assume that mint is behind are exactly at the ray origin
		// if the ray is outside the box, then mint would be the correct choice
		its.t = maxt;
		its.p = ray(its.t);
		return true;
	}
	return false;
}


// includes emission, geometry term and transmittance
V3d Scene::sample_attenuated_directlight(LightSample& ls, RNGd& rng , bool debug)const
{
	// sample lightsource (including geometry term)
	V3d light_over_pdf = light->sample(ls);

	if(debug)
	{
		std::cout << "\t\tlight_over_pdf=" <<  light_over_pdf.toString() << std::endl;
	}

	double maxt = ls.distance;

	// TEMP
	//V3d transmittance_over_pdf( std::exp(-10.0*ls.distance) );
	//V3d light_over_pdf( 1.0 );


	///*
	// find next intersection point ---
	{
		Intersection its;
		if( intersect( Ray3d( ls.refP, -ls.d ), its ) )
		{
			// there is something between current position and lightsource
			//TODO: we currently assume its just the volume boundary and we are save to ignore it
			maxt = std::min( its.t, ls.distance );
		}else
		{
			if(debug)
			{
				std::cout << "\t\t no intersection" << std::endl;
			}
			// currently we expect to be within a volume
			return V3d(0.0f, 0.0f, 0.0f);
		}
	}

	//NB: the geometry term between light source position and current vertex is already been dealt with during light sampling

	// sample transmittance ---
	V3d transmittance_over_pdf = integrator->sample_transmittance( this, Ray3d( ls.refP, -ls.d ), maxt, rng );
	if(debug)
	{
		std::cout << "\t\ttransmittance_over_pdf=" <<  transmittance_over_pdf.toString() << std::endl;
	}
	//*/
	//return light_over_pdf;
	return light_over_pdf.cwiseProduct(transmittance_over_pdf);
}
