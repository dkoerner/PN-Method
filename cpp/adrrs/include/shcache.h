#pragma once
#include <util/data.h>
#include <util/threadpool.h>
#include <util/moexp.h>
#include <math/transform.h>
#include <math/color.h>

#include <scene.h>
#include <volume.h>
#include <integrator.h>
#include <cache.h>



struct SHCache : public Cache
{
	SHCache():
		Cache(),
		m_override(false),
		m_doZeroScattering(true)
	{

	}

	SHCache( const std::string& filename ):
		Cache(),
		m_override(false),
		m_doZeroScattering(true)
	{
		//load(filename);
	}

	virtual Color3f eval( const P3d& pWS, const V3d& d, bool debug = false )const override;

	/*
	const double& sample( int i, int j, int k )const
	{
		int voxelIndex = k*m_resolution.x()*m_resolution.y() + j*m_resolution.x() + i;
		return m_data[voxelIndex * m_numComponentsPerVoxel];
	}
	*/

	const Color3f* get_voxel_data( int i, int j, int k )const;
	Color3f* get_voxel_data( int i, int j, int k );

	void generate( const std::string& filename, const Scene* scene, int order, int numSamples, int res );
	void generate_images( const std::string& filename, const Scene* scene, int order, int numSamples, int res );
	void sample_images( MonteCarloTaskInfo& ti, const Scene* scene );

	EnvMap m_current_map;
	int m_current_voxel;
	void sample_image( MonteCarloTaskInfo& ti, const Scene* scene );


//private:
	/*
	void save(const std::string& filename)
	{
		// dump all data to disk
		{
			std::ofstream out( filename.c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::trunc );

			// write dimension
			out.write( (const char *)&m_dimension, sizeof(int) );

			// write res
			out.write( (const char *)&m_resolution, sizeof(V3i) );

			// write numMoments
			out.write( (const char *)&m_numMoments, sizeof(int) );

			// write localToWorld
			out.write( (const char*)m_localToWorld.getMatrix().data(), sizeof(double)*16 );

			// write data
			out.write( (const char*)m_data.data(), sizeof(double)*m_data.size() );
		}


		// save bgeo voxelgrid of the red color channel of the zero moment (fluence)
		{
			field::VoxelGridField<double>::Ptr fluence_grid = std::make_shared<field::VoxelGridField<double>>(m_resolution);

			// set the first color channel
			double* ptr = m_data.data();
			for( int i=0;i<m_numVoxels;++i, ptr+=m_numComponentsPerVoxel )
				fluence_grid->grid.m_data[i] = *ptr;

			field::write(filename + ".fluence.bgeo", field::xform<double>(fluence_grid, m_localToWorld), m_resolution, m_localToWorld );
		}
	}

	void load(const std::string& filename)
	{
		std::cout << "PNCache::load: loading from file " << filename << std::endl;
		std::ifstream in( filename.c_str(), std::ios_base::in | std::ios_base::binary );
		if(!in.good())
		{
			throw std::runtime_error("PNCache::load unable to load...");
		}

		// read dimension
		in.read( (char *)&m_dimension, sizeof(int) );

		// read res
		in.read( (char *)&m_resolution, sizeof(V3i) );
		m_numVoxels = m_resolution.x()*m_resolution.y()*m_resolution.z();

		// read numMoments
		in.read( (char *)&m_numMoments, sizeof(int) );

		// read localToWorld
		M44d m;
		in.read( (char *)m.data(), sizeof(double)*16 );
		m_localToWorld = Transformd(m);

		// read data
		m_numComponentsPerVoxel = 0;
		for( int m=0;m<m_numMoments;++m )
			m_numComponentsPerVoxel += tensor::numComponents(m);
		m_data.resize(m_numVoxels*m_numComponentsPerVoxel);
		in.read( (char*)m_data.data(), sizeof(double)*m_data.size() );
	}
	*/

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


	void sample_SHCoefficients( MonteCarloTaskInfo& ti, const Scene* scene );


	V3i m_resolution;
	int m_order; // order
	int m_evaluationOrder;
	int m_numVoxels;
	int m_numSHCoefficientsPerVoxel;

	Transformd m_localToWorld;

	std::vector<Color3f> m_data; // sh coefficients (real sh)

	bool m_override;
	bool m_doZeroScattering;
};
