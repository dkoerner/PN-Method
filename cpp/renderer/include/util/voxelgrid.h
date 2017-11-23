#pragma once

#include <math/common.h>
#include <math/vector.h>

#include <fstream>
#include <memory>



template<typename T>
struct VoxelGrid
{
	typedef std::shared_ptr< VoxelGrid<T> > Ptr;
	typedef std::shared_ptr< const VoxelGrid<T> > CPtr;

	static Ptr        create( Vector3i res = Vector3i(10) );

	//static Ptr        load( const std::string &filename );
	void              load( const std::string &filename );
	void              save( const std::string &filename );

	VoxelGrid();
	VoxelGrid( const std::string& filename );

	T                 evaluate(const P3d &vsP )const;
	T                 sample( int i, int j, int k )const;
	T&                lvalue( int i, int j, int k );

	Vector3i          getResolution()const;
	void              getCoordFromIndex(int index, int& x, int& y, int& z)const;
	P3d               localToVoxel(const P3d &p)const;
	P3d               voxelToLocal(const P3d &p)const;
	void              resize( int x, int y, int z );
	void              resize( Vector3i resolution );
	void              fill(const T &value );

	T*                getRawPointer();
	const T*          getRawPointer()const;

	Vector3f          m_sampleLocation; // sample location within one voxel (allows staggered grids etc.)
	Vector3i          m_resolution;
	std::vector<T>    m_data;
	static const int  m_dataType; // e.g. float, double, v3f etc.
};


template<typename T>
typename VoxelGrid<T>::Ptr VoxelGrid<T>::create( Vector3i res )
{
	VoxelGrid<T>::Ptr field = std::make_shared<VoxelGrid<T>>();
	field->resize(res);
	return field;
}

template<typename T>
VoxelGrid<T>::VoxelGrid() : m_sampleLocation(0.5f)
{
	resize( Vector3i(1) );
}

template<typename T>
VoxelGrid<T>::VoxelGrid( const std::string& filename ) : m_sampleLocation(0.5f)
{
	load(filename);
}

/*
template<typename T>
typename VoxelGrid<T>::Ptr VoxelGrid<T>::load( const std::string &filename )
{
	VoxelGrid<T>::Ptr field = VoxelGrid<T>::create();

	// load resolution, data from file
	std::ifstream in( filename.c_str(), std::ios_base::in | std::ios_base::binary );

	if( !in.good() )
		return field;

	in.read( (char *)&field->m_resolution, sizeof(int)*3 );
	int dataType = 0;
	in.read( (char *)&dataType, sizeof(int) );

	if( dataType != m_dataType )
	{
		std::cerr << "VoxelGrid<T>::load: error: datatype in " << filename << " doesnt match.\n";
		return VoxelGrid<T>::Ptr();
	}

	int size = field->m_resolution.x()*field->m_resolution.y()*field->m_resolution.z();
	field->m_data.resize( size );
	in.read( (char *)&field->m_data[0], size*sizeof(T) );

	return field;
}
*/

template<typename T>
void VoxelGrid<T>::load( const std::string &filename )
{
	// load resolution, data from file
	std::ifstream in( filename.c_str(), std::ios_base::in | std::ios_base::binary );

	if( !in.good() )
		throw std::runtime_error("VoxelGrid<T>::load failed to open stream");

	in.read( (char *)&m_resolution, sizeof(int)*3 );
	int dataType = 0;
	in.read( (char *)&dataType, sizeof(int) );

	if( dataType != m_dataType )
		throw std::runtime_error( "VoxelGrid<T>::load: error: datatype in doesnt match." );

	resize(m_resolution);
	in.read( (char *)&m_data[0], m_resolution.x()*m_resolution.y()*m_resolution.z()*sizeof(T) );
}

template<typename T>
void VoxelGrid<T>::save( const std::string &filename )
{
	std::ofstream out( filename.c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::trunc );

	// save resolution, data to file
	out.write( (const char *)&m_resolution, sizeof(int)*3 );
	out.write( (const char *)&m_dataType, sizeof(int) );
	out.write( (const char *)&m_data[0], sizeof(T)*m_resolution.x()*m_resolution.y()*m_resolution.z() );
}

template<typename T>
void VoxelGrid<T>::resize( int x, int y, int z )
{
	resize(Vector3i(x, y, z));
}


template<typename T>
void VoxelGrid<T>::resize( Vector3i resolution )
{
	m_resolution = resolution;
	m_data.resize(m_resolution.x()*m_resolution.y()*m_resolution.z());
	memset( &m_data[0], 0, m_resolution.x()*m_resolution.y()*m_resolution.z()*sizeof(T));
}

template<typename T>
void VoxelGrid<T>::fill( const T& value )
{
	std::fill( m_data.begin(), m_data.end(), value );
}


template<typename T>
T VoxelGrid<T>::sample( int i, int j, int k )const
{
	//int index = k*m_resolution.x()*m_resolution.y() + j*m_resolution.x() + i;

	// numpy indexing (k ist fastest)
	int index =  i*m_resolution[1]*m_resolution[2] + j*m_resolution[2] + k;
	return m_data[index];
}

template<typename T>
T &VoxelGrid<T>::lvalue( int i, int j, int k )
{
	//int index = k*m_resolution.x()*m_resolution.y() + j*m_resolution.x() + i;

	// numpy indexing (k ist fastest)
	int index =  i*m_resolution[1]*m_resolution[2] + j*m_resolution[2] + k;
	return m_data[index];
}

template<typename T>
T VoxelGrid<T>::evaluate( const P3d &vsP )const
{
	typedef double real_t;
	typedef Vector3d Vector;

	Vector vs = vsP;

	// take sample location within voxel into account
	vs -= m_sampleLocation.cast<real_t>();

	real_t tx = vs.x() - floor(vs.x());
	real_t ty = vs.y() - floor(vs.y());
	real_t tz = vs.z() - floor(vs.z());

	// lower left corner
	Vector3i c1;
	c1[0] = (int)floor(vs.x());
	c1[1] = (int)floor(vs.y());
	c1[2] = (int)floor(vs.z());

	// upper right corner
	Vector3i c2 = c1+Vector3i(1);
	Vector3i res = getResolution();

	// clamp the indexing coordinates
	c1[0] = std::max(0, std::min(c1.x(), res.x()-1));
	c2[0] = std::max(0, std::min(c2.x(), res.x()-1));
	c1[1] = std::max(0, std::min(c1.y(), res.y()-1));
	c2[1] = std::max(0, std::min(c2.y(), res.y()-1));
	c1[2] = std::max(0, std::min(c1.z(), res.z()-1));
	c2[2] = std::max(0, std::min(c2.z(), res.z()-1));

	//lerp...
	return lerp( lerp( lerp( sample( c1.x(), c1.y(), c1.z() ),
							 sample( c2.x(), c1.y(), c1.z() ), (real_t)tx ),
					   lerp( sample( c1.x(), c2.y(), c1.z() ),
							 sample( c2.x(), c2.y(), c1.z() ), (real_t)tx ), (real_t)ty ),
				 lerp( lerp( sample( c1.x(), c1.y(), c2.z() ),
							 sample( c2.x(), c1.y(), c2.z() ), (real_t)tx ),
					   lerp( sample( c1.x(), c2.y(), c2.z() ),
							 sample( c2.x(), c2.y(), c2.z() ), (real_t)tx ), (real_t)ty ), (real_t)tz );
}

template<typename T>
Vector3i VoxelGrid<T>::getResolution()const
{
	return m_resolution;
}

template<typename T>
void VoxelGrid<T>::getCoordFromIndex( int index, int& x, int& y, int& z ) const
{
	div_t divresult;

	divresult = div( index, m_resolution.x()*m_resolution.y() );

	z = divresult.quot;
	divresult = div( divresult.rem, m_resolution.x() );
	y = divresult.quot;
	x = divresult.rem;
}

template<typename T>
P3d VoxelGrid<T>::localToVoxel( const P3d& p )const
{
	return P3d( p.x()*m_resolution.x(),
				p.y()*m_resolution.y(),
				p.z()*m_resolution.z());
}

template<typename T>
P3d VoxelGrid<T>::voxelToLocal(const P3d &p)const
{
	return P3d( p.x()/m_resolution.x(),
				p.y()/m_resolution.y(),
				p.z()/m_resolution.z());
}



template<typename T>
T *VoxelGrid<T>::getRawPointer()
{
	return &m_data[0];
}

template<typename T>
const T *VoxelGrid<T>::getRawPointer()const
{
	return &m_data[0];
}

// get maximum value
template<typename T>
T field_maximum( const VoxelGrid<T> &field )
{
	return *std::max_element( field.m_data.begin(), field.m_data.end() );
}



typedef VoxelGrid<float> VoxelGridf;
typedef VoxelGrid<float> ScalarVoxelGrid;
typedef VoxelGrid<Vector3f> VoxelGrid3f;
typedef VoxelGrid<Vector4f> VoxelGrid4f;
typedef VoxelGrid<double> VoxelGridd;
typedef VoxelGrid<Vector3d> VoxelGrid3d;





