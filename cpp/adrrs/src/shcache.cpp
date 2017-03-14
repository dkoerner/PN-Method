#include <shcache.h>



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

