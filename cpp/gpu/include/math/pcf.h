#pragma once



#include <vector>
#include <memory>
#include <algorithm>
#include <fstream>

#include <svpt/math/common.h>


// piecewise constant function
template<typename T, typename R=float>
struct PCF
{
	typedef std::shared_ptr<PCF<T>> Ptr;

	PCF()
		: m_numSamples(0)
	{
	}


	PCF( const PCF& other ):
		m_domain(other.m_domain),
		m_values(other.m_values),
		m_numSamples(other.m_numSamples)
	{
	}

	~PCF()
	{
	}

	void getDomainRange( T& minx, T& maxx )
	{
		if(m_numSamples==0)
		{
			minx = std::numeric_limits<T>::max();
			maxx = -std::numeric_limits<T>::max();
		}else
		{
			minx = m_domain[0];
			maxx = m_domain[m_numSamples-1];
		}
	}


	/*
	void getValueRange( R& miny, R& maxy )
	{
		if(m_numSamples==0)
		{
			miny = std::numeric_limits<T>::max();
			maxy = -std::numeric_limits<T>::max();
		}else
		{
			std::pair<typename std::vector<T>::iterator, typename std::vector<T>::iterator> minmax = std::minmax_element( m_values.begin(), m_values.end() );
			miny = *minmax.first;
			maxy = *minmax.second;
		}
	}
	*/

	void clear()
	{
		m_domain.clear();
		m_values.clear();
		m_numSamples = 0;
	}

	int getNumPoints()const
	{
		return m_numSamples;
	}
	int size()const
	{
		return m_numSamples;
	}
	void addSample( T x, R y )
	{
		m_domain.push_back(x);
		m_values.push_back(y);
		m_numSamples = (int)m_domain.size();
	}

	void setPosition( int pointIndex, T x )
	{
		m_domain[pointIndex] = x;
	}

	void setValue( int pointIndex, R y )
	{
		m_values[pointIndex] = y;
	}


	R evaluate( T x )const
	{
		if( m_numSamples == 0 )
			return R(0.0f);

		// out of bound cases
		if( x < m_domain[0] )
			return m_values[0];
		if( x >= m_domain[m_numSamples-1] )
			return m_values[m_numSamples-1];

		// find interval using binary search http://en.wikipedia.org/wiki/Binary_search_algorithm#Deferred_detection_of_equality
		int imin = 0;
		int imax = m_numSamples - 1;

		while( imin < imax )
		{
			int imid = (imax + imin)/2 + 1;
			if( m_domain[imid] > x )
				imax = imid - 1;
			else
				imin = imid;
		};

		return m_values[imax];
	}

	void save( const std::string& filename )const
	{
		std::ofstream out( filename.c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::trunc );

		// number of samples
		out.write( (const char *)&m_numSamples, sizeof(int) );
		// domain
		out.write( (const char *)&m_domain[0], sizeof(T)*m_numSamples );
		// values
		out.write( (const char *)&m_values[0], sizeof(R)*m_numSamples );
	}
	static Ptr load( const std::string &filename )
	{

		Ptr plf = std::make_shared<PCF<T, R>>();

		// load resolution, data from file
		std::ifstream in( filename.c_str(), std::ios_base::in | std::ios_base::binary );

		if( !in.good() )
			return plf;

		in.read( (char *)&plf->m_numSamples, sizeof(int));
		plf->m_domain.resize( plf->m_numSamples );
		in.read( (char *)&plf->m_domain[0], plf->m_numSamples*sizeof(T) );
		plf->m_values.resize( plf->m_numSamples );
		in.read( (char *)&plf->m_values[0], plf->m_numSamples*sizeof(R) );
		return plf;
	}

	std::vector<T>                      m_domain;
	std::vector<R>                      m_values;
	int                                 m_numSamples;
};


typedef PCF<float, float> PCFf;
typedef PCF<float, V3f> PCF3f;
typedef PCF<double, V3d> PCF3d;
typedef PCF<double, double> PCFd;
