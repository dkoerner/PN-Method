#pragma once
#include <cuda.h>
#include "cumath/cumath.h"


struct CudaPixelBuffer// : public CudaData<CudaPixelBuffer>
{
	HOST CudaPixelBuffer() : /*CudaData<CudaPixelBuffer>(),*/ m_width(0), m_height(0), m_data(0)
	{
	}
	HOST virtual ~CudaPixelBuffer()
	{
		if(m_data)
		{
			CudaSafeCall( cudaDeviceSynchronize() );
			CudaSafeCall( cudaFree( m_data ) );
		}
	}

	// initializes fields and uploads host data if needed
	HOST void initialize( int width, int height, unsigned char *host_data = 0 )
	{
		m_width = width;
		m_height = height;

		CudaSafeCall( cudaDeviceSynchronize() );
		CudaSafeCall( cudaFree( m_data ) );
		int size =  m_width*m_height*sizeof(cumath::V3f);

		CudaSafeCall( cudaMalloc( &m_data, size ) );

		if(host_data)
			CudaSafeCall( cudaMemcpy( m_data, host_data, size, cudaMemcpyHostToDevice ));
		else
			CudaSafeCall( cudaMemset( m_data, 0, size ));
		CudaSafeCall( cudaDeviceSynchronize() );
	}

	HOST void download( unsigned char *host_ptr )
	{
		int size =  m_width*m_height*sizeof(cumath::V3f);
		CudaSafeCall( cudaMemcpy( host_ptr, m_data, size, cudaMemcpyDeviceToHost ));
	}

	DEVICE cumath::V3f& lvalue( int index )
	{
		return m_data[index];
	}
	DEVICE cumath::V3f& lvalue( int i, int j )
	{
		return m_data[j*m_width + i];
	}


	DEVICE const cumath::V3f& value( int i, int j )const
	{
		return m_data[j*m_width + i];
	}


	DEVICE cumath::V3f eval( const cumath::V2f &vsP )const
	{
		// voxelgrid ---
		cumath::Vec2<float> vs = vsP;

		// voxels defined at cell centers
		vs.x -= 0.5f;
		vs.y -= 0.5f;

		float tx = vs.x - floor(vs.x);
		float ty = vs.y - floor(vs.y);

		// lower left corner
		cumath::V2i c1;
		c1.x = (int)floor(vs.x);
		c1.y = (int)floor(vs.y);

		// upper right corner
		cumath::V2i c2 = c1+cumath::V2i(1,1);

		// clamp the indexing coordinates
		c1.x = cumath::max(0, cumath::min(c1.x, m_width-1));
		c2.x = cumath::max(0, cumath::min(c2.x, m_width-1));
		c1.y = cumath::max(0, cumath::min(c1.y, m_height-1));
		c2.y = cumath::max(0, cumath::min(c2.y, m_height-1));

		//lerp...
		return cumath::lerp( cumath::lerp( value( c1.x, c1.y ),
										   value( c2.x, c1.y ), tx ),
							 cumath::lerp( value( c1.x, c2.y ),
										   value( c2.x, c2.y ), tx ), ty );
	}

	void clear()
	{
		CudaSafeCall( cudaDeviceSynchronize() );
		CudaSafeCall( cudaMemset( m_data, 0, m_width*m_height*sizeof(cumath::V3f) ));
	}

	HOST_DEVICE int getWidth()const
	{
		return m_width;
	}

	HOST_DEVICE int getHeight()const
	{
		return m_height;
	}

	int                         m_width;
	int                         m_height;
	cumath::V3f*                m_data; // device ptr
	int                         m_pixelSize;
};

