#include <integrators/pnispt.h>




V3d PNISPT::Li( const Scene* scene, RadianceQuery& rq, RNGd& rng )const
{
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
		//TraceInfo ti;
		TraceInfoPNIS ti;
		ti.depth = 0;
		ti.current_vertex = Vertex();
		P3d pWS = rq.ray(mint+Epsilon);
		if(rq.volume)
			ti.current_vertex.setPosition(pWS, rq.volume->evalExtinction(pWS, rq.debug), rq.volume->evalAlbedo(pWS));
		else
			ti.current_vertex.setPosition(pWS, V3d(0.0, 0.0, 0.0), V3d(0.0, 0.0, 0.0));
		ti.current_direction = rq.ray.d;
		ti.scene = scene;
		ti.debug = rq.debug;
		ti.debug = rq.debug;
		ti.debug_info = &debug_info;
		ti.pns = m_pns.get();
		return trace( ti, rng );
	}else
	{
		// no intersection with the medium boundary
	}

	return V3d(0.0, 0.0, 0.0);
}




/*
V3d PNISPT::Li( const Scene* scene, RadianceQuery& rq, RNGd& rng )const
{

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
		//TraceInfoPNIS ti;
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
		//ti.debug = rq.debug;
		//ti.debug_info = &debug_info;
		//ti.pns = m_pns.get();
		return trace( ti, rng );
	}else
	{
		// no intersection with the medium boundary
	}

	return V3d(0.0, 0.0, 0.0);
}
*/

