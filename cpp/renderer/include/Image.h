#pragma once
#include <string>
#include <memory>
#include <math/vector.h>
#include <math/rng.h>
#include <math/dpdf.h>








struct Image// : public Eigen::Array<Eigen::Vector3d, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
{
	typedef std::shared_ptr<Image> Ptr;
	typedef Eigen::Array<Eigen::Vector3d, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Base;

	// construction
	Image(const Vector2i &size = Vector2i(0, 0));
	Image(const std::string& filename);
	Image(const Image& image);

	// evaluation
	Eigen::Vector3d& pixel( int x, int y );
	const Eigen::Vector3d& pixel( int x, int y )const;
	V3d eval( const P2d& pRaster )const;

	// io
	const Eigen::Vector3d* data()const;
	Eigen::Vector3d* data();
	Base& getArray();
	void save(const std::string& filename); // stores content as exr file
	std::string toString()const;

	// info
	V2i getResolution()const;
	static double luminance( const V3d& c );
	P2d uvToRaster(const P2d& uv)const;
	P2d rasterToUV(const P2d& pRaster)const;

	// manipulation

private:
	Base m_data;
};

Image::Ptr blur_image( Image::Ptr image, double stddev );

struct ImageSampler
{
	typedef std::shared_ptr<ImageSampler> Ptr;

	ImageSampler( Image::Ptr image )
	{
		buildPDF(image);
	}

	P2d sample( double& pdf, RNGd& rng ); // samples raster position
	P2d sample(double &pdf, P2d sample);

	/*
	P2d directionToUV( const V3d& d )const;
	V3d uvToDirection( const P2d& uv )const;
	*/
	P2d rasterToUV(const P2d& xy)const;
	/*
	P2d uvToRaster(const P2d& uv)const;

	V3d rasterToDirection( const P2d& xy )const;
	P2d directionToRaster( const V3d& d )const;
	*/

	void buildPDF(Image::Ptr image);

private:
	V2i m_resolution;
	DiscretePDFd m_colDpdf;
	std::vector<DiscretePDFd> m_rowDpdf;

};
