#pragma once
#include <util/data.h>
#include <util/threadpool.h>
#include <math/transform.h>
#include <math/color.h>

#include <scene.h>
#include <volume.h>
#include <integrator.h>
#include <cache.h>









namespace tensor
{
	int ipow(int base, int exp);

//// compute data index
//int index = 0;
//for( int i=0;i<m_order;++i )
//	index += indices[m_order-i-1]*ipow(m_dimension, i);
//return m_data[index];

	int numComponents( int order, int dimension=3 );

	template<typename T>
	struct Tensor
	{
		Tensor(T* data=0, int order=0, int dimension=3):
			m_data(data),
			m_order(order),
			m_dimension(dimension)
		{

		}

		struct Iterator
		{
			Iterator(T* ptr, int order, int dimension, int default_index = 0):
				m_ptr(ptr),
				m_order(order),
				m_dimension(dimension),
				m_indices(order, default_index)
			{
			}

			// returns value of the index with the specified index
			int index( int index_index )
			{
				return m_indices[index_index];
			}

			// returns pointer to the current component
			T* ptr()
			{
				return m_ptr;
			}

			// returns value of the current component
			T& value()
			{
				return *m_ptr;
			}

			double weight( const V3d& d )
			{
				double weight = 1.0;
				for( int j=0;j<m_order;++j )
					weight*=d[m_indices[j]];
				return weight;
			}

			void advance()
			{
				++m_ptr;

				int increment = 1;
				for( int i=0;i<m_order;++i )
				{
					m_indices[i] += increment;

					// check if index is smaller than dimension
					if(m_indices[i] < m_dimension)
						return;

					// otherwise we need to handle carry
					increment = m_indices[i]/m_dimension;
					m_indices[i] = m_indices[i] % m_dimension;
				}

				return;
			}

			Iterator& operator ++ ()
			{
				advance();
				return *this;
			}
			Iterator operator ++ (int)
			{
				Iterator it = *this;
				it.advance();
				return it;
			}

			bool operator==(const Iterator& other)
			{
				return m_ptr == other.m_ptr;
			}

			bool operator!=(const Iterator& other)
			{
				return m_ptr != other.m_ptr;
			}

			void print()
			{
				for( int i=0;i<m_order;++i )
				{
					std::cout << m_indices[i];
				}
				std::cout << std::endl;
			}

			std::string toString()const
			{
				std::string result = "";

				result += "index=";
				for( int i=0;i<m_order;++i )
					result += ::toString(m_indices[i]);
				return result;
			}

			std::string index_str()const
			{
				std::string result = "";
				//if(m_order==0)
				//	return "0";
				for( int i=0;i<m_order;++i )
					result += ::toString(m_indices[i]);
				return result;
			}

		private:
			T* m_ptr;
			std::vector<int> m_indices;
			int m_order;
			int m_dimension;
		};

		Iterator begin()
		{
			return Iterator(m_data, m_order, m_dimension);
		}

		Iterator end()
		{
			return Iterator(m_data+numComponents(m_order, m_dimension), m_order, m_dimension, m_dimension-1);
		}


	private:
		T* m_data;
		int m_order;
		int m_dimension;
	};

	/*
	template<typename T>
	void add( const SymmetricTensor<T>& a, const SymmetricTensor<T>& b )
	{

	}

	template<typename T>
	void contract( const SymmetricTensor<T>& a, int idx_a, const SymmetricTensor<T>& b, int idx_b )
	{

	}
	*/
} // namespace tensor



struct PNCache : public Cache
{
	PNCache():
		Cache(),
		m_dimension(3),
		m_override(false),
		m_doZeroScattering(true)
	{

	}

	PNCache( const std::string& filename ):
		Cache(),
		m_dimension(3),
		m_override(false),
		m_doZeroScattering(true)
	{
		load(filename);
	}

	const double& sample( int i, int j, int k )const
	{
		int voxelIndex = k*m_resolution.x()*m_resolution.y() + j*m_resolution.x() + i;
		return m_data[voxelIndex * m_numComponentsPerVoxel];
	}

	const double* get_voxel_data( int i, int j, int k )const;
	double* get_voxel_data( int i, int j, int k );

	virtual Color3f eval(const P3d& pWS, const V3d& d, bool debug = false)const override;

	bool m_override;
	bool m_doZeroScattering;

	void generate( const std::string& filename, const Scene* scene, int numMoments, int numSamples, int res );


	tensor::Tensor<double> getMoment( int i, int j, int k, int moment );
private:
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


	void sampleMoments( MonteCarloTaskInfo& ti, const Scene* scene );

	int m_dimension;
	V3i m_resolution;
	int m_numMoments; // number of moments this cache holds per voxel
	Transformd m_localToWorld;
	std::vector<double> m_data; // moment tensor components for all voxels, all moments, rgb color

	int m_numVoxels;
	int m_numComponentsPerVoxel;
};
