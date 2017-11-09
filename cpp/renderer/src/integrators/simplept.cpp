#include <integrators/simplept.h>






// returns sigma_t at sampled position (is invalid when we exceeded maxt)
double delta_tracking( const Scene* scene, const Ray3d& ray, double maxt, int component, RNGd& rng, V3d& sigma_t, bool debug )
{
	double sigma_t_max = scene->volume->getMaxExtinction()[component];

	double t = 0.0;
	while(true)
	{
		double step = -log( 1.0-rng.next1D() )/sigma_t_max;
		t += step;

		if(t>= maxt)
			break;

		sigma_t = scene->volume->evalExtinction(ray(t), debug);

		// russian roulette
		if(rng.next1D()<sigma_t[component]/sigma_t_max)
			break;
	}

	return t;
}



V3d SimplePT::Li( const Scene* scene, RadianceQuery& rq, RNGd& rng )const
{
/*
	return V3d( double(rq.pixel[0])/double(scene->camera->getResolutionX()),
				double(rq.pixel[1])/double(scene->camera->getResolutionY()),
				0.0);
*/


	//rq.transmittance = V3d(1.0, 1.0, 1.0);

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
		P3d pWS = rq.ray(mint+Epsilon);
		if(rq.volume)
			ti.current_vertex.setPosition(pWS, rq.volume->evalExtinction(pWS), rq.volume->evalAlbedo(pWS));
		else
			ti.current_vertex.setPosition(pWS, V3d(0.0, 0.0, 0.0), V3d(0.0, 0.0, 0.0));
		ti.current_direction = rq.ray.d;
		ti.scene = scene;
		ti.debug = rq.debug;
		return trace( ti, rng );
	}else
	{
		// no intersection with the medium boundary
	}

	return V3d(0.0, 0.0, 0.0);
}
