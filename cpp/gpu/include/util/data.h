#pragma once
#include <sys/stat.h>
#include <fstream>
#include <memory>
#include <math/vector.h>
#include <math/plf.h>

#include <util/bitmap.h>
#include <camera.h>




inline bool file_exists(const std::string& name)
{
	struct stat buffer;
	return (stat (name.c_str(), &buffer) == 0);
}


// utility stuff ===================================
struct Data
{
	enum Type
	{
		EUnknown,
		ERadial,
		EImage
	};

	void addSensor( const P3d& sensor )
	{
		sensors.push_back(sensor);
		estimates.push_back(0.0);
	}

	void addRay( const Ray3d& ray, double estimate = 0.0 )
	{
		sensors.push_back(ray.o);
		directions.push_back(ray.d);
		estimates.push_back(estimate);
		rays.push_back(ray);
	}

	void addRay( const Ray3d& ray, const Ray3d& ray2, double estimate = 0.0 )
	{
		sensors.push_back(ray.o);
		directions.push_back(ray.d);
		estimates.push_back(estimate);
		rays.push_back(ray2);
	}

	int numElements()const
	{
		return int(sensors.size());
	}

	void clear()
	{
		m_type = Type::EUnknown;
		sensors.clear();
		directions.clear();
		estimates.clear();
	}

	std::vector<P3d>    sensors;
	std::vector<V3d>    directions;
	std::vector<Ray3d>  rays;
	std::vector<double> estimates;

	// information related to the datatype
	Type m_type;
	int m_xres; // image resolution in x (in case type is image)
	int m_yres; // image resolution in y (in case type is image)
	P3d m_o; // origin of radial data

	// type dependent save function
	void save( const std::string& filename )const;
};



// utility stuff ===================================
template<typename T = double>
struct Data2
{
	enum Type
	{
		EUnknown,
		ERadial,
		EImage
	};

	void addSensor( const P3d& sensor )
	{
		sensors.push_back(sensor);
		estimates.push_back(0.0);
	}

	void addRay( const Ray3d& ray, double estimate = 0.0 )
	{
		sensors.push_back(ray.o);
		directions.push_back(ray.d);
		estimates.push_back(estimate);
		rays.push_back(ray);
	}

	void addRay( const Ray3d& ray, const Ray3d& ray2, double estimate = 0.0 )
	{
		sensors.push_back(ray.o);
		directions.push_back(ray.d);
		estimates.push_back(estimate);
		rays.push_back(ray2);
	}

	int numElements()const
	{
		return int(sensors.size());
	}

	void clear()
	{
		m_type = Type::EUnknown;
		sensors.clear();
		directions.clear();
		estimates.clear();
	}

	std::vector<P3d>    sensors;
	std::vector<V3d>    directions;
	std::vector<Ray3d>  rays;
	std::vector<T> estimates;

	// information related to the datatype
	Type m_type;
	int m_xres; // image resolution in x (in case type is image)
	int m_yres; // image resolution in y (in case type is image)
	P3d m_o; // origin of radial data
};


struct Slice
{
	Slice() : o(0.0, 0.0, 0.0), right(1.0, 0.0, 0.0), up(0.0, 1.0, 0.0)
	{}
	Slice( P3d o, double width, double height, V3d right=V3d(1.0, 0.0, 0.0), V3d up=V3d(0.0, 1.0, 0.0) ) :
		o(o),
		right(right),
		up(up),
		width(width),
		height(height)
	{}
	P3d o;
	V3d right;
	V3d up;
	double width;
	double height;
};


void initializeData(Camera::Ptr camera, Data& data );
void initializeData( const Slice& slice, int resX, int resY, Data& data );

template<class Data>
void initializeData(const std::vector<double>& radii, Data& data, const V3d& d=V3d(1.0, 0.0, 0.0), const P3d& o = P3d(0.0, 0.0, 0.0))
{
	data.clear();
	data.m_type = Data::Type::ERadial;
	data.m_o = o;
	for( auto r : radii )
		data.addSensor(o+r*d);
}

template<class Data>
void initializeData(double minRadius, double maxRadius, double dr, Data& data)
{
	data.clear();
	int numBins = int((maxRadius-minRadius)/dr)+1;

	data.addSensor(V3d(0.01,0.0,0.0));
	data.addSensor(V3d(0.02,0.0,0.0));
	data.addSensor(V3d(0.03,0.0,0.0));
	data.addSensor(V3d(0.04,0.0,0.0));
	data.addSensor(V3d(0.05,0.0,0.0));
	data.addSensor(V3d(0.06,0.0,0.0));
	for( int i=0;i<numBins;++i )
	{
		double r = minRadius+dr*0.5 + i*dr;
		data.addSensor(V3d(r,0.0,0.0));
	}
}

template<class Data>
void initializeData(double minRadius, double maxRadius, int num, Data& data)
{
	data.clear();
	int numBins = num;
	double dr = (maxRadius-minRadius)/double(numBins);
	for( int i=0;i<numBins;++i )
	{
		double r = minRadius+dr*0.5 + i*dr;
		data.addSensor(V3d(r,0.0,0.0));
	}
}

template<typename T>
void writeSamples( const std::string& filename, const std::vector<T>& x, const std::vector<T>& y, const std::string& sep = " " )
{
	std::cout << "writing samples to: " << filename << std::endl;
	std::ofstream out_values( filename.c_str(), std::ios_base::out | std::ios_base::trunc );
	int numSamples = int(x.size());
	for( int i=0;i<numSamples;++i )
		out_values << x[i] << sep << y[i] << std::endl;
}


template<typename T>
void writeSamplesMATLAB(const std::string &filename, const std::vector<T>& samples, const std::string& id = "samples" )
{
	std::ofstream f(filename);

	f << id << " = [ ";
	int numSamples = int(samples.size());
	for (int j=0; j<numSamples; ++j)
	{
		f << samples[j];
		if (j+1 < numSamples)
			f << ", ";
	}
	f << " ];" << std::endl;
	f << id << " = transpose(" << id << ");" << std::endl;

/*
		<< "colormap(jet);" << std::endl
		<< "clf; subplot(2,1,1);" << std::endl
		<< "imagesc(obsFrequencies_svpt);" << std::endl
		<< "title('Observed frequencies');" << std::endl
		<< "axis equal;" << std::endl
		<< "subplot(2,1,2);" << std::endl
		<< "imagesc(expFrequencies_svpt);" << std::endl
		<< "axis equal;" << std::endl
		<< "title('Expected frequencies');" << std::endl;
*/
	//f.close();
}

template<typename T>
void readSamples2( const std::string& filename, std::vector<T>& values, int& numRows, int& numCols, char delim = ' ')
{
	std::cout << "reading samples from: " << filename << std::endl;
	std::ifstream in_values( filename.c_str(), std::ios_base::in );

	values.clear();

	numRows = 0;
	numCols = 0;

	if(!in_values.good())
		throw std::runtime_error("readSamples2: unable to read file " + filename );

	std::string line = "";
	while (getline(in_values, line))
	{
		int numCols_cur = 0;
		std::stringstream strstr(line);
		std::string word = "";
		while (getline(strstr,word, delim))
		{
			values.push_back(fromString<T>(word));
			++numCols_cur;
		}

		if(numRows == 0)
			numCols = numCols_cur;
		else
		if(numCols_cur != numCols)
		{
			std::cout << "readSamples2: warning: varying number of columns per row.\n";
		}

		++numRows;
	}

	/*
	T x, y;
	while (in_values >> x)
	{
		// process pair (a,b)
		x_list.push_back(x);
		y_list.push_back(y);
	}
	*/
}


template<typename T>
void readSamples( const std::string& filename, std::vector<T>& x_list, std::vector<T>& y_list )
{
	std::cout << "reading samples from: " << filename << std::endl;
	std::ifstream in_values( filename.c_str(), std::ios_base::in );

	T x, y;
	while (in_values >> x >> y)
	{
		// process pair (a,b)
		x_list.push_back(x);
		y_list.push_back(y);
	}
}

template<typename T>
void readSamples( const std::string& filename, PLF<T, T>& plf )
{
	std::vector<T> x_list;
	std::vector<T> y_list;
	readSamples( filename, x_list, y_list );

	int numSamples = int(x_list.size());
	plf.clear();
	for( int i = 0;i<numSamples;++i )
		plf.addSample( x_list[i], y_list[i] );
}

template<typename T>
void writeSamples( const std::string& filename, const std::vector<T>& y, const std::string& sep = " " )
{
	std::vector<T> x;
	for( int i=0;i<int(y.size());++i )
		x.push_back(i);
	writeSamples( filename, x, y, sep );
}

std::shared_ptr<Bitmap> toImage(const std::vector<double> &samples, int resX, int resY);
void writeImage( const std::vector<double>& samples, int resX, int resY, const std::string &filename );

template<typename T>
void saveRadialData( const std::vector<P3d>& positions, const std::vector<T>& values, const std::string& filename, const P3d& o = P3d(0.0, 0.0, 0.0) )
{
	std::vector<T> radii;
	for( auto& p:positions )
	{
		T r = (p-o).norm();
		radii.push_back(r);
	}

	writeSamples( filename, radii, values);
}

template<typename T>
void linearSamples( T min, T max, int num, std::vector<T>& samples )
{
	samples.clear();
	if( num == 1 )
	{
		samples.push_back(min);
	}else
	{
		T dr = (max-min)/T(num-1);
		for( int i=0;i<num;++i )
		{
			//T r = minRadius+dr*T(0.5) + i*dr;
			T r = min + i*dr;
			samples.push_back(r);
		}
	}
}

template<typename T>
std::vector<T> linearSamples( T min, T max, int num )
{
	std::vector<T> samples;
	if( num == 1 )
	{
		samples.push_back(min);
	}else
	{
		T dr = (max-min)/T(num-1);
		for( int i=0;i<num;++i )
		{
			//T r = minRadius+dr*T(0.5) + i*dr;
			T r = min + i*dr;
			samples.push_back(r);
		}
	}
	return samples;
}


template<typename T>
void logSamples( T minRadius, T maxRadius, int num, std::vector<T>& samples, T base = T(10.0) )
{
	linearSamples<T>( minRadius, maxRadius, num, samples );
	for( auto& s:samples )
		s = std::pow(base, s);
}

template<typename T>
void writeHistogram( const std::string& filename, const std::vector<T>& samples, int numBins )
{
	auto minmax = std::minmax_element( samples.begin(), samples.end() );
	T min = *minmax.first;
	T max = *minmax.second;
	T dbin = (max-min)/T(numBins);

	//int numSamples = 0;
	std::vector<double> bins(numBins, 0);
	for( auto sample:samples )
	{
		int bin = int((sample-min)/dbin);
		if((bin >=0)&&(bin<numBins))
		{
			bins[bin] += 1.0;
			//++numSamples;
		}
	}

	// normalize
	for( int i=0;i<numBins;++i )
		bins[i] /= T(samples.size());

	std::vector<double> x_list;
	std::vector<double> y_list;

	for( int i=0;i<numBins;++i )
	{
		x_list.push_back( min + i*dbin );
		y_list.push_back( bins[i] );
	}
	writeSamples( filename, x_list, y_list );


}

/*
void logSamples(double minRadius, double maxRadius, int num, std::vector<double>& samples , double base = 10.0);
void natlogSamples( double minRadius, double maxRadius, int num, std::vector<double>& samples );

void saveRadialData(const Data& data, const std::string& filename , const P3d &refP = P3d(0.0, 0.0, 0.0));

*/



//http://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
template <typename T, typename Compare>
std::vector<std::size_t> sort_permutation(
	const std::vector<T>& vec,
	Compare& compare)
{
	std::vector<std::size_t> p(vec.size());
	std::iota(p.begin(), p.end(), 0);
	std::sort(p.begin(), p.end(),
		[&](std::size_t i, std::size_t j){ return compare(vec[i], vec[j]); });
	return p;
}

template <typename T>
std::vector<T> apply_permutation(
	const std::vector<T>& vec,
	const std::vector<std::size_t>& p)
{
	std::vector<T> sorted_vec(p.size());
	std::transform(p.begin(), p.end(), sorted_vec.begin(),
		[&](std::size_t i){ return vec[i]; });
	return sorted_vec;
}

// ==============================================================================




void xport_observations(const std::string &filename, std::vector<std::pair<P3d, double>>& obervations);
void xport_line(const std::string &filename, P3d& a, P3d& b);











