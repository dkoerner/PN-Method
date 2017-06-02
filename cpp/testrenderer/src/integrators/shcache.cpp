#include <util/data.h>
#include <util/moexp.h>
#include <math/transform.h>
#include <math/color.h>

#include <volume.h>
#include <cache.h>


struct SHCache : public Cache
{
	typedef std::shared_ptr<SHCache> Ptr;

	SHCache( int order, V3i resolution, Transformd localToWorld ):
		Cache(),
		m_order(order),
		m_evaluationOrder(order),
		m_resolution(resolution),
		m_localToWorld(localToWorld),
		m_numSHCoefficientsPerVoxel(moexp::numSHCoefficients(order)),
		m_numVoxels(resolution.x()*resolution.y()*resolution.z())
	{
		m_data.resize(m_numVoxels*m_numSHCoefficientsPerVoxel);
		std::cout << "SHCache with " << m_data.size()*sizeof(Color3f) << " bytes\n";
	}

	SHCache( const std::string& filename ):
		Cache()
	{
		load(filename);
	}

	virtual std::string getId()const override
	{
		return toString<int>(m_evaluationOrder);
	}

	virtual Color3f eval( const P3d& pWS, const V3d& d, bool debug = false )const override;

	Color3f* getCoefficients( int voxel_index )
	{
		return &m_data[voxel_index*m_numSHCoefficientsPerVoxel];
	}





	void save(const std::string& filename)
	{
		std::cout << "SHCache::save: writing to file " << filename << std::endl;

		// dump all data to disk
		std::ofstream out( filename.c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::trunc );

		// write resolution
		out.write( (const char *)&m_resolution, sizeof(V3i) );

		// write order of the sh coefficient expansion
		out.write( (const char *)&m_order, sizeof(int) );

		// write localToWorld
		out.write( (const char*)m_localToWorld.getMatrix().data(), sizeof(double)*16 );

		// write data
		out.write( (const char*)m_data.data(), sizeof(Color3f)*m_data.size() );
	}

//private:
	void load(const std::string& filename)
	{
		std::cout << "SHCache::load: loading from file " << filename << std::endl;

		std::ifstream in( filename.c_str(), std::ios_base::in | std::ios_base::binary );
		if(!in.good())
		{
			throw std::runtime_error("SHCache::load unable to load...");
		}

		// read order of the sh coefficient expansion
		in.read( (char *)&m_resolution, sizeof(V3i) );
		m_numVoxels = m_resolution.x()*m_resolution.y()*m_resolution.z();

		// read order of the sh coefficient expansion
		in.read( (char *)&m_order, sizeof(int) );
		m_numSHCoefficientsPerVoxel = moexp::numSHCoefficients(m_order);
		m_evaluationOrder = m_order;

		// read localToWorld
		M44d m;
		in.read( (char *)m.data(), sizeof(double)*16 );
		m_localToWorld = Transformd(m);

		// read data
		m_data.resize(m_numVoxels*m_numSHCoefficientsPerVoxel);
		in.read( (char*)m_data.data(), sizeof(Color3f)*m_data.size() );

		std::cout << "resolution=" << m_resolution.toString() << std::endl;
		std::cout << "m_order=" << m_order<< std::endl;

	}

	void getCoordFromIndex( int index, int& x, int& y, int& z ) const
	{
		div_t divresult;

		divresult = div( index, m_resolution.x()*m_resolution.y() );

		z = divresult.quot;
		divresult = div( divresult.rem, m_resolution.x() );
		y = divresult.quot;
		x = divresult.rem;
	}

	P3d voxelToLocal(const P3d &p)const
	{
		return P3d( p.x()/m_resolution.x(),
					p.y()/m_resolution.y(),
					p.z()/m_resolution.z());
	}

	P3d localToVoxel( const P3d& p )const
	{
		return P3d( p.x()*m_resolution.x(),
					p.y()*m_resolution.y(),
					p.z()*m_resolution.z());
	}

	const Color3f* get_voxel_data( int i, int j, int k )const;
	Color3f* get_voxel_data( int i, int j, int k );




	int m_order; // order
	int m_evaluationOrder;

	V3i m_resolution;
	Transformd m_localToWorld;

	int m_numVoxels;
	int m_numSHCoefficientsPerVoxel;

	std::vector<Color3f> m_data; // sh coefficients (real sh)
};



template<typename T>
struct SHField : public Field<T>
{
	typedef std::shared_ptr<SHField<T>> Ptr;

	SHField( SHCache::Ptr cache ):
		Field<T>(),
		m_cache(cache)
	{

	}

	static Ptr create( SHCache::Ptr cache )
	{
		return std::make_shared<SHField>( cache );
	}

	virtual T eval( const P3d& p, bool debug = false )const override
	{
		// we evaluate the first sh band only, so direction doesnt matter
		V3d d(0.0, 0.0, 1.0);
		return m_cache->eval(p, d, debug).r();
	}

	virtual void getValueRange( T& min, T& max )const override
	{
	}

	SHCache::Ptr m_cache;
};



struct FluenceCache : public Cache
{
	typedef std::shared_ptr<FluenceCache> Ptr;


	FluenceCache( const std::string& filename ) : Cache()
	{
		m_fluence_field = field::read<double>(filename);
	}


	virtual Color3f eval( const P3d& pWS, const V3d& d, bool debug = false )const override
	{
		return m_fluence_field->eval(pWS)*INV_FOURPI;
		//return Color3f(0.0f);
	}

	virtual std::string getId()const override
	{
		return "fluencecache";
	}

private:
	Fieldd::Ptr m_fluence_field;
};










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

	// pick the first voxel
	//const Color3f* ptr = get_voxel_data(c1[0], c1[1], c1[2]);

	///*
	// linearly interpolate sh coefficients
	std::vector<Color3f> interpolated_shcoeffs(m_numSHCoefficientsPerVoxel, 0.0);
	const Color3f* ptr = interpolated_shcoeffs.data();
	{
		// evaluate and accumulate moments ---
		SomeHelper<Color3f> helper( interpolated_shcoeffs.data(), m_numSHCoefficientsPerVoxel );
		helper.multiply_add( get_voxel_data( c1.x(), c1.y(), c1.z() ), (1.0-tx)*(1.0-ty)*(1.0-tz) );
		helper.multiply_add( get_voxel_data( c2.x(), c1.y(), c1.z() ), tx*(1.0-ty)*(1.0-tz) );
		helper.multiply_add( get_voxel_data( c1.x(), c2.y(), c1.z() ), (1.0-tx)*ty*(1.0-tz) );
		helper.multiply_add( get_voxel_data( c2.x(), c2.y(), c1.z() ), tx*ty*(1.0-tz) );
		helper.multiply_add( get_voxel_data( c1.x(), c1.y(), c2.z() ), (1.0-tx)*(1.0-ty)*tz );
		helper.multiply_add( get_voxel_data( c2.x(), c1.y(), c2.z() ), tx*(1.0-ty)*tz );
		helper.multiply_add( get_voxel_data( c1.x(), c2.y(), c2.z() ), (1.0-tx)*ty*tz );
		helper.multiply_add( get_voxel_data( c2.x(), c2.y(), c2.z() ), tx*ty*tz );
	}
//*/



	Color3f gg = moexp::Y_real_sum<Color3f>( m_evaluationOrder, ptr, theta_phi.x(), theta_phi.y(), debug );
	return gg;

}










// load sh cache -----------------------
//SHCache::Ptr shcache = std::make_shared<SHCache>("nebulae200.shcache");

// verification
/*
for( int l=0;l<=shcache.m_order;++l )
{
	std::string filename = "sh_reconstruction_alt_$0.bgeo";
	filename = replace(filename, "$0", toString(l));

	Color3f* coeffs = shcache.getCoefficients(0);
	rasterizeSphericalFunctionSphere(filename, [&](double theta, double phi)->Color3f
	{
		return moexp::Y_real_sum<Color3f>(l, coeffs, theta, phi);
	}, 8.0);
}
*/


// pn analysis and debugging ---
/*
{
	Integrator::Ptr integrator = integrators::dummy();
	scene.integrator = integrator.get();

	PNCache cache;
	cache.m_override = true;
	cache.m_doZeroScattering = false;
	int numMoments = 4;
	int numSamples = 1000000;
	int res = 1;
	cache.generate( outpath + "/analysis", &scene, numMoments, numSamples, res );


	Wedge wedge;
	std::vector<int> moment_values = {0, 1, 2, 3, 4};
	wedge.addParm("moment", moment_values);
	wedge.build();

	std::cout << "4pi=" << 4.0*M_PI << std::endl;
	std::cout << "4pi/3=" << 4.0*M_PI/3.0 << std::endl;

	std::vector<Wedge::Iteration> iterations = wedge.iterations();
	for( auto it : iterations )
	{
		int moment_index = it.getInt("moment");
		std::cout << "moment=" << moment_index << std::endl;

		tensor::Tensor<double> moment_tensor = cache.getMoment( 0, 0, 0, moment_index );

		for(auto component = moment_tensor.begin(), end=moment_tensor.end(); component!=end;++component)
		{
			std::cout << "\t" << component.index_str() << "=" << component.value() << std::endl;

//				houio::Geometry::Ptr geo = houio::Geometry::createSphere(50, 50, 1.0);
//				houio::Attribute::Ptr pattr = geo->getAttr("P");
//				int numPoints = pattr->numElements();
//				for(int i=0;i<numPoints;++i)
//				{
//					houio::math::V3f d = pattr->get<houio::math::V3f>(i).normalized();
//					d *= std::abs(component.weight(V3d(d.x, d.y, d.z)));
//					pattr->set<houio::math::V3f>(i, d);
//				}
//				houio::HouGeoIO::xport(it.expand_value( outpath+"/analysis_$0_" + component.index_str() + ".bgeo"), geo);
		}
	}
}
*/
