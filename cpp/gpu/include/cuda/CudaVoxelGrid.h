#pragma once
#include <cuda.h>
#include <cuda_runtime.h>
#include "cumath/cumath.h"


template<typename T>
struct CudaVoxelGrid
{
	HOST CudaVoxelGrid();
	HOST ~CudaVoxelGrid();



	HOST                        void initialize( const cumath::V3i &resolution, T *host_data = 0 );

	DEVICE cumath::V3f          localToVoxel(const cumath::V3f &pLS);
	DEVICE cumath::V3f          voxelToLocal(const cumath::V3f &pVS);

	DEVICE const T&             sample( int i, int j, int k )const; // returns discrete voxel value
	DEVICE T&                   lvalue( int i, int j, int k ); // returns discrete voxel value reference
	DEVICE T                    eval( const cumath::V3f &vsP )const; // returns linearly interpolated value at voxepspace position


	HOST void download( unsigned char *host_ptr );

	/*




	HOST void                   resize( const cumath::V3i &resolution );
	HOST void                   setBound( const cumath::Box3f &bound );
	HOST void                   fill( const T &value );

	HOST_DEVICE cumath::V3i     getResolution()const;
	HOST_DEVICE cumath::V3f     getVoxelSize()const; // returns voxelsize in worldspace




	cumath::M44f                m_localToVoxel; // transforms from worldspace to voxelspace
	cumath::M44f                m_voxelToLocal;



	cumath::Box3f               m_bound;
	*/
	T*                          m_data; // voxeldata (deviceptr)
	cumath::V3i                 m_resolution;
};



template<typename T>
HOST CudaVoxelGrid<T>::CudaVoxelGrid() :
	m_data(0)
{
}
template<typename T>
HOST CudaVoxelGrid<T>::~CudaVoxelGrid()
{
	if(m_data)
	{
		CudaSafeCall( cudaDeviceSynchronize() );
		CudaSafeCall( cudaFree( m_data ) );
		m_data = 0;
	}
}



template<typename T>
HOST void CudaVoxelGrid<T>::initialize( const cumath::V3i &resolution, T *host_data )
{
	m_resolution = resolution;

	CudaSafeCall( cudaDeviceSynchronize() );
	CudaSafeCall( cudaFree( m_data ) );
	int size = m_resolution.x*m_resolution.y*m_resolution.z*sizeof(T);

	CudaSafeCall( cudaMalloc( &m_data, size ) );

	if(host_data)
		CudaSafeCall( cudaMemcpy( m_data, host_data, size, cudaMemcpyHostToDevice ));
	else
		CudaSafeCall( cudaMemset( m_data, 0, size ));
}

template<typename T>
DEVICE cumath::V3f CudaVoxelGrid<T>::localToVoxel(const cumath::V3f &pLS)
{
	return cumath::V3f( pLS.x*m_resolution.x,
						pLS.y*m_resolution.y,
						pLS.z*m_resolution.z);
}

template<typename T>
DEVICE cumath::V3f CudaVoxelGrid<T>::voxelToLocal(const cumath::V3f &pVS)
{
	return cumath::V3f( pVS.x/m_resolution.x,
						pVS.y/m_resolution.y,
						pVS.z/m_resolution.z);
}


template<typename T>
DEVICE const T& CudaVoxelGrid<T>::sample( int i, int j, int k )const
{
	return m_data[k*m_resolution.x*m_resolution.y + j*m_resolution.x + i];
}

template<typename T>
DEVICE T& CudaVoxelGrid<T>::lvalue( int i, int j, int k )
{
	return m_data[k*m_resolution.x*m_resolution.y + j*m_resolution.x + i];
}


template<typename T>
DEVICE T CudaVoxelGrid<T>::eval( const cumath::V3f &vsP )const
{
	// voxelgrid ---
	cumath::Vec3<float> vs = vsP;

	// voxels defined at cell centers
	vs.x -= 0.5f;
	vs.y -= 0.5f;
	vs.z -= 0.5f;

	float tx = vs.x - floor(vs.x);
	float ty = vs.y - floor(vs.y);
	float tz = vs.z - floor(vs.z);

	// lower left corner
	cumath::V3i c1;
	c1.x = (int)floor(vs.x);
	c1.y = (int)floor(vs.y);
	c1.z = (int)floor(vs.z);

	// upper right corner
	cumath::V3i c2 = c1+cumath::V3i(1,1,1);

	// clamp the indexing coordinates
	c1.x = cumath::max(0, cumath::min(c1.x, m_resolution.x-1));
	c2.x = cumath::max(0, cumath::min(c2.x, m_resolution.x-1));
	c1.y = cumath::max(0, cumath::min(c1.y, m_resolution.y-1));
	c2.y = cumath::max(0, cumath::min(c2.y, m_resolution.y-1));
	c1.z = cumath::max(0, cumath::min(c1.z, m_resolution.z-1));
	c2.z = cumath::max(0, cumath::min(c2.z, m_resolution.z-1));

	//lerp...
	return cumath::lerp( cumath::lerp( cumath::lerp( sample( c1.x, c1.y, c1.z ),
													 sample( c2.x, c1.y, c1.z ), tx ),
									   cumath::lerp( sample( c1.x, c2.y, c1.z ),
													 sample( c2.x, c2.y, c1.z ), tx ), ty ),
						 cumath::lerp( cumath::lerp( sample( c1.x, c1.y, c2.z ),
													 sample( c2.x, c1.y, c2.z ), tx ),
									   cumath::lerp( sample( c1.x, c2.y, c2.z ),
													 sample( c2.x, c2.y, c2.z ), tx ), ty ), tz );
}

template<typename T>
HOST void CudaVoxelGrid<T>::download( unsigned char *host_ptr )
{
	int size =  m_resolution.x*m_resolution.y*m_resolution.z*sizeof(T);
	CudaSafeCall( cudaMemcpy( host_ptr, m_data, size, cudaMemcpyDeviceToHost ));
}


/*
template<typename T>
HOST void CudaVoxelgrid<T>::resize( const cumath::V3i &resolution )
{
	CudaSafeCall( cudaDeviceSynchronize() );
	CudaSafeCall( cudaFree( m_data ) );
	m_resolution = resolution;
	int size = sizeof( T )*resolution.x*resolution.y*resolution.z;
	CudaSafeCall( cudaMalloc( &m_data, size ) );
	CudaSafeCall( cudaMemset( m_data, 0, size ));

	m_worldToVoxel = m_worldToLocal*cumath::M44f().scale( cumath::V3f(m_resolution) );
	m_voxelToWorld = m_worldToVoxel.inverse();
}

template<typename T>
HOST void CudaVoxelgrid<T>::setBound( const cumath::Box3f &bound )
{
	m_bound = bound;
	cumath::V3f dim = m_bound.size();

	m_localToWorld = cumath::M44f().scale( cumath::V3f(dim.x, dim.y, dim.z) )*cumath::M44f().translate(m_bound.minPoint);
	m_worldToLocal = m_localToWorld.inverse();
	m_worldToVoxel = m_worldToLocal*cumath::M44f().scale( cumath::V3f((float)m_resolution.x, (float)m_resolution.y, (float)m_resolution.z) );
	m_voxelToWorld = m_worldToVoxel.inverse();
}

template<typename T>
HOST void CudaVoxelgrid<T>::fill( const T &value )
{
	// TODO
}

template<typename T>
cumath::V3i CudaVoxelgrid<T>::getResolution()const
{
	return m_resolution;
}

// returns voxelsize in worldspace
template<typename T>
HOST_DEVICE cumath::V3f CudaVoxelgrid<T>::getVoxelSize()const
{
	cumath::V3f v0 = cumath::transform( cumath::V3f(0.0f,0.0f, 0.0f), m_voxelToWorld );
	cumath::V3f v1 = cumath::transform( cumath::V3f(1.0f,1.0f, 1.0f), m_voxelToWorld );
	return cumath::V3f( v1.x - v0.x, v1.y - v0.y, v1.z - v0.z );
}


*/
