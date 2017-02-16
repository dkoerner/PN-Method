#include <pncache.h>






namespace tensor
{
	int ipow(int base, int exp)
	{
		int result = 1;
		while (exp)
		{
			if (exp & 1)
				result *= base;
			exp >>= 1;
			base *= base;
		}

		return result;
	}

	int numComponents( int order, int dimension )
	{
		return ipow(dimension, order);
	}
}


void PNCache::generate( const std::string& filename, const Scene* scene, int numMoments, int numSamples, int res )
{
	std::cout << "PNCache::generate..."<< std::endl;

	bool override = false;
	if( file_exists(filename) && !override )
	{
		std::cout << "PNCache::generate file " << filename << "already exists!"<<std::endl;
		throw std::runtime_error("PNCache::generate file " + filename + " already exists!");
	}

	m_resolution = V3i(res, res, res);
	m_numVoxels = m_resolution.x()*m_resolution.y()*m_resolution.z();
	m_numMoments = numMoments;
	m_localToWorld = scene->volume->getLocalToWorld();


	// create data
	m_numComponentsPerVoxel = 0;
	for( int m=0;m<m_numMoments;++m )
		m_numComponentsPerVoxel += tensor::numComponents(m);
	m_data.resize(m_numVoxels*m_numComponentsPerVoxel);

	// sample and accumulate moments (parallel version) ---
	Terminator terminator(numSamples);
	auto f = std::bind(&PNCache::sampleMoments, this, std::placeholders::_1, scene );
	runGenericTasks<MonteCarloTaskInfo>( f, terminator, ThreadPool::getNumSystemCores() );

	// save data to disk ---
	save(filename);
}


void PNCache::sampleMoments( MonteCarloTaskInfo& ti, const Scene* scene )
{
	int numVoxels = m_resolution.x()*m_resolution.y()*m_resolution.z();
	for (int voxel = ti.taskid; voxel< numVoxels; voxel += ti.numTasks)
	{
		// work out voxel coordinate
		int i,j, k;
		getCoordFromIndex(voxel, i, j, k);

		// sample position within voxel
		P3d pVS(i+ti.rng.next1D(), j+ti.rng.next1D(), k+ti.rng.next1D());
		P3d pLS = voxelToLocal(pVS);
		P3d pWS = m_localToWorld*pLS;

		// sample and accumulate direct light ---
		double L_direct_sample = 0.0;
		V3d d_direct;
		{
			LightSample ls;
			ls.refP = pWS;
			Color3f attenuated_light_over_pdf = scene->sample_attenuated_directlight( ls, ti.rng );
			L_direct_sample = attenuated_light_over_pdf.getLuminance();
			d_direct = -ls.d;
		}

		// sample and accumulate indirect light ---
		double L_indirect_sample = 0.0;
		V3d d_indirect;
		{
			// sample direction
			d_indirect = sampleSphere<double>(ti.rng);
			double d_pdf = sampleSpherePDF();

			// integrate (produce radiance field sample)
			RadianceQuery rq;
			rq.ray = Ray3d(pWS, d_indirect);
			L_indirect_sample = scene->integrator->Li( scene, rq, ti.rng ).getLuminance()/d_pdf;
		}

		// accumulate ---
		int moment_start_component = 0;
		for( int m=0;m<m_numMoments;++m )
		{
			double* current_moment_data = m_data.data() + voxel*m_numComponentsPerVoxel + moment_start_component;

			tensor::Tensor<double> t(current_moment_data, m);
			for( auto it = t.begin(), end=t.end();it != end; ++it )
			{
				// compute weights of current moment component
				double weight_direct = it.weight(d_direct);
				double weight_indirect = it.weight(d_indirect);
				double sample = weight_direct*L_direct_sample + weight_indirect*L_indirect_sample;

				// accumulate weighted radiance field value
				double& component = *it.ptr();
				component += (sample - component)/float(ti.samples+1);
			}

			moment_start_component += tensor::numComponents(m);
		} // for each moment
	}

	++ti.samples;
}


Color3f PNCache::eval( const P3d& pWS )const
{
	P3d pLS = m_localToWorld.inverse()*pWS;
	P3d pVS = localToVoxel(pLS);

	// take sample location within voxel into account
	pVS -= V3d(0.5, 0.5, 0.5);

	double tx = pVS.x() - floor(pVS.x());
	double ty = pVS.y() - floor(pVS.y());
	double tz = pVS.z() - floor(pVS.z());

	// lower left corner
	V3i c1;
	c1[0] = (int)floor(pVS.x());
	c1[1] = (int)floor(pVS.y());
	c1[2] = (int)floor(pVS.z());

	// upper right corner
	V3i c2 = c1+V3i(1);
	V3i res = m_resolution;

	// clamp the indexing coordinates
	c1[0] = std::max(0, std::min(c1.x(), res.x()-1));
	c2[0] = std::max(0, std::min(c2.x(), res.x()-1));
	c1[1] = std::max(0, std::min(c1.y(), res.y()-1));
	c2[1] = std::max(0, std::min(c2.y(), res.y()-1));
	c1[2] = std::max(0, std::min(c1.z(), res.z()-1));
	c2[2] = std::max(0, std::min(c2.z(), res.z()-1));

	// interpolate all moment coefficients ---
	double fluence = lerp( lerp( lerp( sample( c1.x(), c1.y(), c1.z() ),
									   sample( c2.x(), c1.y(), c1.z() ), tx ),
								 lerp( sample( c1.x(), c2.y(), c1.z() ),
									   sample( c2.x(), c2.y(), c1.z() ), tx ), ty ),
						   lerp( lerp( sample( c1.x(), c1.y(), c2.z() ),
									   sample( c2.x(), c1.y(), c2.z() ), tx ),
								 lerp( sample( c1.x(), c2.y(), c2.z() ),
									   sample( c2.x(), c2.y(), c2.z() ), tx ), ty ), tz );

	// evaluate and accumulate moments ---
	return Color3f(fluence, fluence, fluence);
}

Color3f PNCache::eval(const P3d& pWS, const V3d &d)const
{
	P3d pLS = m_localToWorld.inverse()*pWS;
	P3d pVS = localToVoxel(pLS);

	// take sample location within voxel into account
	pVS -= V3d(0.5, 0.5, 0.5);

	double tx = pVS.x() - floor(pVS.x());
	double ty = pVS.y() - floor(pVS.y());
	double tz = pVS.z() - floor(pVS.z());

	// lower left corner
	V3i c1;
	c1[0] = (int)floor(pVS.x());
	c1[1] = (int)floor(pVS.y());
	c1[2] = (int)floor(pVS.z());

	// upper right corner
	V3i c2 = c1+V3i(1);
	V3i res = m_resolution;

	// clamp the indexing coordinates
	c1[0] = std::max(0, std::min(c1.x(), res.x()-1));
	c2[0] = std::max(0, std::min(c2.x(), res.x()-1));
	c1[1] = std::max(0, std::min(c1.y(), res.y()-1));
	c2[1] = std::max(0, std::min(c2.y(), res.y()-1));
	c1[2] = std::max(0, std::min(c1.z(), res.z()-1));
	c2[2] = std::max(0, std::min(c2.z(), res.z()-1));

	// interpolate all moment coefficients ---
	double fluence = lerp( lerp( lerp( sample( c1.x(), c1.y(), c1.z() ),
									   sample( c2.x(), c1.y(), c1.z() ), tx ),
								 lerp( sample( c1.x(), c2.y(), c1.z() ),
									   sample( c2.x(), c2.y(), c1.z() ), tx ), ty ),
						   lerp( lerp( sample( c1.x(), c1.y(), c2.z() ),
									   sample( c2.x(), c1.y(), c2.z() ), tx ),
								 lerp( sample( c1.x(), c2.y(), c2.z() ),
									   sample( c2.x(), c2.y(), c2.z() ), tx ), ty ), tz );

	// evaluate and accumulate moments ---
	return Color3f(fluence, fluence, fluence)*INV_FOURPI;
	//return Color3f(0.0)*INV_FOURPI;
}
