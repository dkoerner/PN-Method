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

