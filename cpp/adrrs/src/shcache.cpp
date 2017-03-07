#include <shcache.h>
#include <util/moexp.h>



void SHCache::generate( const std::string& filename, const Scene* scene, int order, int numSamples, int res )
{
	std::cout << "SHCache::generate..."<< std::endl;

	bool override = m_override;
	if( file_exists(filename) && !override )
	{
		std::cout << "SHCache::generate file " << filename << "already exists!"<<std::endl;
		throw std::runtime_error("SHCache::generate file " + filename + " already exists!");
	}

	m_resolution = V3i(res, res, res);
	m_numVoxels = m_resolution.x()*m_resolution.y()*m_resolution.z();
	m_order = order;
	m_evaluationOrder = m_order;

	m_localToWorld = scene->volume->getLocalToWorld();

	// create data
	m_numSHCoefficientsPerVoxel = moexp::numSHCoefficients(m_order);
	m_data.resize(m_numVoxels*m_numSHCoefficientsPerVoxel);

	// sample and accumulate moments (parallel version) ---
	Terminator terminator(numSamples);
	auto f = std::bind(&SHCache::sample_SHCoefficients, this, std::placeholders::_1, scene );
	runGenericTasks<MonteCarloTaskInfo>( f, terminator, ThreadPool::getNumSystemCores() );

	// save data to disk ---
	//save(filename);
}



void SHCache::generate_images( const std::string& filename, const Scene* scene, int order, int numSamples, int res )
{
	std::cout << "SHCache::generate images..."<< std::endl;

	/*
	bool override = m_override;
	if( file_exists(filename) && !override )
	{
		std::cout << "SHCache::generate file " << filename << "already exists!"<<std::endl;
		throw std::runtime_error("SHCache::generate file " + filename + " already exists!");
	}
	*/

	m_resolution = V3i(res, res, res);
	m_numVoxels = m_resolution.x()*m_resolution.y()*m_resolution.z();
	m_order = order;
	m_evaluationOrder = m_order;

	m_localToWorld = scene->volume->getLocalToWorld();

	// create data
	m_numSHCoefficientsPerVoxel = moexp::numSHCoefficients(m_order);
	m_data.resize(m_numVoxels*m_numSHCoefficientsPerVoxel);

	/*
	// sample and accumulate moments (parallel version) ---
	Terminator terminator(numSamples);
	auto f = std::bind(&SHCache::sample_SHCoefficients, this, std::placeholders::_1, scene );
	runGenericTasks<MonteCarloTaskInfo>( f, terminator, ThreadPool::getNumSystemCores() );
*/

	// sample and accumulate moments (parallel version) ---
	Terminator terminator(1);
	auto f = std::bind(&SHCache::sample_images, this, std::placeholders::_1, scene );
	runGenericTasks<MonteCarloTaskInfo>( f, terminator, ThreadPool::getNumSystemCores() );


	// save data to disk ---
	//save(filename);
}

void SHCache::sample_SHCoefficients( MonteCarloTaskInfo& ti, const Scene* scene )
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
		Color3f L_direct_sample(0.0);
		/*
		V3d d_direct;
		if(m_doZeroScattering)
		{
			LightSample ls;
			ls.refP = pWS;
			Color3f attenuated_light_over_pdf = scene->sample_attenuated_directlight( ls, ti.rng );
			L_direct_sample = attenuated_light_over_pdf.getLuminance();
			d_direct = -ls.d;
		}
		*/

		// sample and accumulate indirect light ---
		Color3f L_indirect_sample(0.0);
		double theta, phi;
		V3d d_indirect;
		{
			// sample direction
			d_indirect = sampleSphere<double>(ti.rng);
			P2d theta_phi = sphericalCoordinates(d_indirect);
			theta = theta_phi.x();
			phi = theta_phi.y();
			double d_pdf = sampleSpherePDF();

			// integrate (produce radiance field sample)
			RadianceQuery rq;
			rq.ray = Ray3d(pWS, d_indirect);
			L_indirect_sample = scene->integrator->Li( scene, rq, ti.rng )/d_pdf;
		}

		// integrate and accumulate for each sh basis function ---
		Color3f* ptr = get_voxel_data(i, j, k);
		for( int l=0;l<=m_order;++l )
		{
			for( int m=-l;m<=l;++m, ++ptr )
			{
				double sh = moexp::Y_real(l, m, theta, phi);
				Color3f sample = L_indirect_sample*sh;
				Color3f& sh_coeff = *ptr;
				sh_coeff += (sample - sh_coeff)/float(i+1);
			}
		}
	} // for each voxel...

	++ti.samples;
}


void SHCache::sample_images( MonteCarloTaskInfo& ti, const Scene* scene )
{
	int numVoxels = m_resolution.x()*m_resolution.y()*m_resolution.z();
	for (int voxel = ti.taskid; voxel< numVoxels; voxel += ti.numTasks)
	{
		if(ti.taskid == 0)
			std::cout << "at voxel " << voxel << "/" << numVoxels << std::endl;

		// work out voxel coordinate
		int v_i,v_j,v_k;
		getCoordFromIndex(voxel, v_i, v_j, v_k);

		// sample position within voxel
		//P3d pVS(v_i+ti.rng.next1D(), v_j+ti.rng.next1D(), v_k+ti.rng.next1D());
		P3d pVS(v_i+0.5, v_j+0.5, v_k+0.5);
		P3d pLS = voxelToLocal(pVS);
		P3d pWS = m_localToWorld*pLS;


		// sample and accumulate direct light ---
		Color3f L_direct_sample(0.0);
		/*
		V3d d_direct;
		if(m_doZeroScattering)
		{
			LightSample ls;
			ls.refP = pWS;
			Color3f attenuated_light_over_pdf = scene->sample_attenuated_directlight( ls, ti.rng );
			L_direct_sample = attenuated_light_over_pdf.getLuminance();
			d_direct = -ls.d;
		}
		*/

		EnvMap map_indirect(V2i(128, 64));
		int resx = map_indirect.bitmap().cols();
		int resy = map_indirect.bitmap().rows();
		int numSamples = 100;


		for(int j=0;j<resy;++j)
			for(int i=0;i<resx;++i)
			{
				Color3f& c = map_indirect.bitmap().coeffRef(j, i);
				c = Color3f(1.0f, 0.0f, 0.0f);
				for( int s=0;s<numSamples;++s )
				{
					V3d d_indirect = map_indirect.xyToDirection( P2d(i+ti.rng.next1D(), j+ti.rng.next1D()) );

					// integrate (produce radiance field sample)
					RadianceQuery rq;
					rq.ray = Ray3d(pWS, d_indirect);

					Color3f L_indirect_sample = scene->integrator->Li( scene, rq, ti.rng );

					//if(i>200)
					//std::cout << L_indirect_sample.toString() << " " << d_indirect.toString() << std::endl;

					//c = L_indirect_sample;
					c += (L_indirect_sample-c)/double(s+1);
				}
			}

		//map_indirect = EnvMap("shcache_images/image_0_0_0.exr");

		double filter_width = 3.0;
		map_indirect.bitmap() = filter(map_indirect.bitmap(), gaussFilter(filter_width));

		std::string filename("shcache_images/image_$0_$1_$2_filtered.exr");
		filename = replace(filename, "$0", toString(v_i));
		filename = replace(filename, "$1", toString(v_j));
		filename = replace(filename, "$2", toString(v_k));
		map_indirect.bitmap().saveEXR(filename);

	} // for each voxel...

	++ti.samples;
}


/*
double contract_moment( tensor::Tensor<const double>& a, const V3d& d  )
{
	double result = 0.0;

	// iterate over all components of tensor a
	for( auto it = a.begin(), end = a.end(); it!=end;++it  )
		result += it.weight(d)*it.value();

	return result;
}
*/

const Color3f* SHCache::get_voxel_data( int i, int j, int k )const
{
	int voxelIndex = k*m_resolution.x()*m_resolution.y() + j*m_resolution.x() + i;
	return &m_data[voxelIndex * m_numSHCoefficientsPerVoxel];
}

Color3f* SHCache::get_voxel_data( int i, int j, int k )
{
	int voxelIndex = k*m_resolution.x()*m_resolution.y() + j*m_resolution.x() + i;
	return &m_data[voxelIndex * m_numSHCoefficientsPerVoxel];
}

/*
template<typename T>
struct SomeHelper
{
	SomeHelper( T* ptr, int size ):
		m_ptr(ptr),
		m_size(size)
	{
	}

	void multiply_add( const T* other_ptr, double coefficient )
	{
		T* ptr = m_ptr;
		for( int i=0;i<m_size;++i, ++ptr, ++other_ptr )
			*ptr += *other_ptr*coefficient;
	}

	void assign( const T* other_ptr )
	{
		T* ptr = m_ptr;
		for( int i=0;i<m_size;++i, ++ptr )
			*ptr = *other_ptr;
	}

private:
	T* m_ptr;
	int m_size;
};
*/

Color3f SHCache::eval(const P3d& pWS, const V3d &d, bool debug)const
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

	P2d theta_phi = sphericalCoordinates(d);
	const Color3f* ptr = get_voxel_data(c1[0], c1[1], c1[2]);
	Color3f gg = moexp::Y_real_sum<Color3f>( m_evaluationOrder, ptr, theta_phi.x(), theta_phi.y() );
	return gg;

	// interpolate all moment coefficients ---
//	double fluence = lerp( lerp( lerp( sample( c1.x(), c1.y(), c1.z() ),
//									   sample( c2.x(), c1.y(), c1.z() ), tx ),
//								 lerp( sample( c1.x(), c2.y(), c1.z() ),
//									   sample( c2.x(), c2.y(), c1.z() ), tx ), ty ),
//						   lerp( lerp( sample( c1.x(), c1.y(), c2.z() ),
//									   sample( c2.x(), c1.y(), c2.z() ), tx ),
//								 lerp( sample( c1.x(), c2.y(), c2.z() ),
//									   sample( c2.x(), c2.y(), c2.z() ), tx ), ty ), tz );
//	return Color3f(fluence)*INV_FOURPI;


//	double result = 0.0;
//	//
//	result += sample(c1.x(), c1.y(), c1.z())*(1.0-tx)*(1.0-ty)*(1.0-tz);
//	result += sample(c2.x(), c1.y(), c1.z())*tx*(1.0-ty)*(1.0-tz);
//	result += sample(c1.x(), c2.y(), c1.z())*(1.0-tx)*ty*(1.0-tz);
//	result += sample(c2.x(), c2.y(), c1.z())*tx*ty*(1.0-tz);
//	result += sample(c1.x(), c1.y(), c2.z())*(1.0-tx)*(1.0-ty)*tz;
//	result += sample(c2.x(), c1.y(), c2.z())*tx*(1.0-ty)*tz;
//	result += sample(c1.x(), c2.y(), c2.z())*(1.0-tx)*ty*tz;
//	result += sample(c2.x(), c2.y(), c2.z())*tx*ty*tz;
//	return Color3f(result)*INV_FOURPI;

	/*
	// evaluate and accumulate moments ---
	std::vector<double> interpolated_moments_data(m_numComponentsPerVoxel, 0.0);

	{
		SomeHelper<double> helper( interpolated_moments_data.data(), m_numComponentsPerVoxel );
		helper.multiply_add( get_voxel_data( c1.x(), c1.y(), c1.z() ), (1.0-tx)*(1.0-ty)*(1.0-tz) );
		helper.multiply_add( get_voxel_data( c2.x(), c1.y(), c1.z() ), tx*(1.0-ty)*(1.0-tz) );
		helper.multiply_add( get_voxel_data( c1.x(), c2.y(), c1.z() ), (1.0-tx)*ty*(1.0-tz) );
		helper.multiply_add( get_voxel_data( c2.x(), c2.y(), c1.z() ), tx*ty*(1.0-tz) );
		helper.multiply_add( get_voxel_data( c1.x(), c1.y(), c2.z() ), (1.0-tx)*(1.0-ty)*tz );
		helper.multiply_add( get_voxel_data( c2.x(), c1.y(), c2.z() ), tx*(1.0-ty)*tz );
		helper.multiply_add( get_voxel_data( c1.x(), c2.y(), c2.z() ), (1.0-tx)*ty*tz );
		helper.multiply_add( get_voxel_data( c2.x(), c2.y(), c2.z() ), tx*ty*tz );
	}

	double* ptr = interpolated_moments_data.data();
	double result = 0.0;

	if(debug)
	{
		std::cout << "numComponentsPerVoxel=" << m_numComponentsPerVoxel << std::endl;
	}

	int numMoments = 1;
	for( int i=0;i<numMoments;++i )
	{
		tensor::Tensor<const double> moment_tensor(ptr, 0);
		double coeff = 1.0;
		if(i==0)
			coeff = INV_FOURPI;
		else
		if(i==1)
			coeff = 3.0*INV_FOURPI;
		else
		if(i==2)
			coeff = 15.0/(8.0*M_PI);

		double t = contract_moment( moment_tensor, d );

		if(debug)
		{
			std::cout << "moment" << i << "=" << t << std::endl;
			std::cout << "numComponents=" << tensor::numComponents(i) << std::endl;
		}

		//if(i==2)
			result += coeff*t;
		ptr+=tensor::numComponents(i);
	}
	// first moment ---
	{
		//tensor::Tensor<const double> moment_tensor(ptr, 0);
		//result += 3.0*INV_FOURPI*contract_moment( moment_tensor, d );
		//ptr+=tensor::numComponents(1);
	}

	// second moment ---

	// third moment ---

	// fourth moment ---


	return Color3f(result);
	*/
	//
}

/*
tensor::Tensor<double> PNCache::getMoment( int i, int j, int k, int moment )
{
	double* ptr = get_voxel_data( i, j, k );
	for( int i=0;i<moment;++i )
		ptr+=tensor::numComponents(i);
	return tensor::Tensor<double>(ptr, moment);
}





*/

