#pragma once
#include <cuda.h>
#include "cumath/common.h"


template<class T>
struct CudaData
{
	__host__ CudaData<T>() : m_this_device(0)
	{
	}

	__host__ virtual ~CudaData<T>()
	{
		if(m_this_device)
			CudaSafeCall( cudaFree(m_this_device) );
	}

	T *synchronize()
	{
		cudaDeviceSynchronize();
		if(!m_this_device)
			CudaSafeCall( cudaMalloc( &m_this_device, sizeof( T ) ) );
		CudaSafeCall( cudaMemcpy( m_this_device, this, sizeof(T), cudaMemcpyHostToDevice ));
		return m_this_device;
	}

	const T *synchronize()const
	{
		cudaDeviceSynchronize();
		if(!m_this_device)
			CudaSafeCall( cudaMalloc( &m_this_device, sizeof( T ) ) );
		CudaSafeCall( cudaMemcpy( m_this_device, this, sizeof(T), cudaMemcpyHostToDevice ));
		return m_this_device;
	}
/*
	void bindToGlobal( const char *name )
	{
		CudaSafeCall( cudaDeviceSynchronize() );
		T *deviceptr = synchronize();

		// update global pointer
		CudaSafeCall( cudaMemcpyToSymbol(name, &deviceptr, sizeof(T*), 0, cudaMemcpyHostToDevice ) );

		CudaSafeCall( cudaDeviceSynchronize() );
	}
*/
	mutable T*                   m_this_device;
};



template<class T>
struct CudaPtr
{
	__host__ CudaPtr<T>( T* host_ptr ) :
		m_host_ptr(host_ptr),
		m_device_ptr(0)
	{
	}

	__host__ virtual ~CudaPtr<T>()
	{
		if(m_device_ptr)
		{
			CudaSafeCall( cudaFree(m_device_ptr) );
			m_device_ptr = 0;
		}
	}

	T *synchronize_device()
	{
		cudaDeviceSynchronize();
		if(!m_device_ptr)
			CudaSafeCall( cudaMalloc( &m_device_ptr, sizeof( T ) ) );
		CudaSafeCall( cudaMemcpy( m_device_ptr, m_host_ptr, sizeof(T), cudaMemcpyHostToDevice ));
		return m_device_ptr;
	}

	T* device()
	{
		if(!m_device_ptr)
			synchronize_device();
		return m_device_ptr;
	}
	T* host()
	{
		return m_host_ptr;
	}
	/*
	T* synchronize_host...
	*/

	const T *synchronize_device()const
	{
		cudaDeviceSynchronize();
		if(!m_device_ptr)
			CudaSafeCall( cudaMalloc( &m_device_ptr, sizeof( T ) ) );
		CudaSafeCall( cudaMemcpy( m_device_ptr, m_host_ptr, sizeof(T), cudaMemcpyHostToDevice ));
		return m_device_ptr;
	}
/*
	void bindToGlobal( const char *name )
	{
		CudaSafeCall( cudaDeviceSynchronize() );
		T *deviceptr = synchronize();

		// update global pointer
		CudaSafeCall( cudaMemcpyToSymbol(name, &deviceptr, sizeof(T*), 0, cudaMemcpyHostToDevice ) );

		CudaSafeCall( cudaDeviceSynchronize() );
	}
*/
	mutable T*                   m_device_ptr;
	mutable T*                   m_host_ptr;
};
