#pragma once
#include <util/data.h>
#include <util/moexp.h>
#include <math/transform.h>
#include <math/color.h>

#include <volume.h>
#include <cache.h>



struct SHCache : public Cache
{
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
