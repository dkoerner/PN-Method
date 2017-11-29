#include <integrators/directpn.h>




#include <houio/HouGeoIO.h>
#include <fields/VoxelGridField.h>

houio::ScalarField::Ptr g_temp;


DirectPN::DirectPN( PNSolution* pns,
					bool doSingleScattering,
					bool doMultipleScattering):
					Integrator(),
					m_pns(pns),
					m_doSingleScattering(doSingleScattering),
					m_doMultipleScattering(doMultipleScattering)
{
	// debug
	/*
	{
		std::string filename = "C:/projects/epfl/epfl17/python/pnsolver/results/nebulae/nebulae_p1_2.pns.bgeo";
		g_temp = houio::HouGeoIO::importVolume(filename);

		P3d pWS( -0.034216, 0.321882, -0.204602 );


		houio::math::V3f pWS2( float(pWS.x()), float(pWS.y()), float(pWS.z()) );
		float e = g_temp->evaluate( g_temp->worldToVoxel(pWS2) );
		std::cout << "houio e=" << e << std::endl;


		double pns_e = m_pns->evalCoefficient(pWS, 0);
		std::cout << "pns e=" << pns_e << std::endl;
	}
	*/

	/*
	int coeff_index = 0;
	int numCoeffs = pns->getNumCoeffs();
	V3i res = pns->getResolution();
	std::vector<double> data(res[0]*res[1]*res[2]);

	for( int i=0;i<res[0];++i )
		for( int j=0;j<res[1];++j )
			for( int k=0;k<res[2];++k )
			{
				int index_dst = i*res[2]*res[1] + j*res[2] + k;
				int index_src = pns->getIndex(V3i(i, j, k))+coeff_index;
				data[index_dst] = pns->data()[index_src];
			}

	VoxelGridField3d t( res, data.data() );
	*/

}



V3d DirectPN::trace( TraceInfo& ti, RNGd& rng )const
{
	V3d L(0.0f, 0.0f, 0.0f);

	if(ti.debug)
		std::cout <<"VolumePathTracer::trace\n";

	// we only do a single bounce
	{
		if(ti.debug)
			std::cout <<"VolumePathTracer::trace depth=" << ti.depth << std::endl;

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
			L += nee(ti, rng);

		// instead of continuing the pathtracing loop in a random direction,
		// we stop and pick up an estimate of indirect light from PN-solution
		// TODO: sample according to SH distribution
		if(m_doMultipleScattering)
		{
			if( ti.current_vertex.getType() == Vertex::EVolume )
			{
				/*
				// sample sphere
				// TODO: sample according to SH distribution
				V3d dir = sampleSphere(rng);
				ti.throughput_over_pdf /= sampleSpherePDF();

				// evaluate indirect light from PN solution
				double e = m_pns->eval(ti.current_vertex.getPosition(), dir);

				V3d phase_times_sigma_s = ti.scene->volume->evalPhase( ti.current_vertex.getPosition(), dir, -ti.current_direction )*
										  ti.current_vertex.m_sigma_s;
				L+= V3d(e).cwiseProduct(phase_times_sigma_s).cwiseProduct(ti.throughput_over_pdf);

				ti.current_direction = dir;
				*/

				/*
				// variation of the previous version. here we integrate over the sphere explicitly
				{
					double zero_coefficient = 0.0;
					for( int i=0;i<1;++i )
					{
						// compute zero coefficient
						V3d dir = sampleSphere(rng);
						double dir_pdf = sampleSpherePDF();
						double e = m_pns->eval(ti.current_vertex.getPosition(), dir)/dir_pdf;
						zero_coefficient += (e-zero_coefficient)/double(i+1);
					}

					V3d phase_times_sigma_s = INV_FOURPI*ti.current_vertex.m_sigma_s;
					double e =zero_coefficient;
					L += V3d(e).cwiseProduct(phase_times_sigma_s).cwiseProduct(ti.throughput_over_pdf);
				}
				*/

				///*
				// constant phase function allows us to evaluate the zero coefficient directly
				P3d pWS = ti.current_vertex.getPosition();
				V3d phase_times_sigma_s = INV_FOURPI*ti.current_vertex.m_sigma_s;
				double e = m_pns->evalCoefficient(pWS, 0);
				L += V3d(e).cwiseProduct(phase_times_sigma_s).cwiseProduct(ti.throughput_over_pdf);
				//*/
			}else
				throw std::runtime_error("surface sampling not implemented");
		}

		++ ti.depth;
	}

	return L;
}

V3d DirectPN::Li( const Scene* scene, RadianceQuery& rq, RNGd& rng )const
{
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

