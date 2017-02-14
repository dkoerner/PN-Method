#pragma once

#include <string>
#include <memory>
#include <math/color.h>
#include <math/vector.h>



/**
 * \brief Stores a RGB high dynamic-range bitmap
 *
 * The bitmap class provides I/O support using the OpenEXR file format
 */
class Bitmap : public Eigen::Array<Color3f, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> {
public:
	typedef std::shared_ptr<Bitmap> Ptr;
	typedef Eigen::Array<Color3f, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Base;

	//
	// \brief Allocate a new bitmap of the specified size
	//
	// The contents will initially be undefined, so make sure
	// to call \ref clear() if necessary
	//
	Bitmap(const Vector2i &size = Vector2i(0, 0))
		: Base(size.y(), size.x())
	{
	}

	Bitmap(const Bitmap& bitmap)
		: Base(bitmap.rows(), bitmap.cols())
	{
		memcpy( data(), bitmap.data(), sizeof(Color3f)*rows()*cols() );
	}


	/// Load an OpenEXR file with the specified filename
	Bitmap(const std::string &filename);

	Color3f eval( const P2d& uv ) const;
	Color3f eval( const P2d& uv, bool debug ) const;

	void getCoordFromIndex( int index, int &x, int &y )const
	{
		getCoordFromIndex(index, cols(), x, y);
	}

	static void getCoordFromIndex( int index, int resX, int &x, int &y )
	{
		div_t divresult  = div( index, resX );
		y = divresult.quot;
		x = divresult.rem;
	}

	P2i uvToRaster( const P2d& uv )const
	{
		P2d xy(
			uv.x()*(cols()-1)-0.5,
			uv.y()*(rows()-1)-0.5
			);

		// lower left corner
		return P2i( int(floor(xy.x())),
					int(floor(xy.y())));
	}

	void resizeInterpolated( int xres, int yres );
	void resizeInterpolated( double scale );

	void saveEXR(const std::string &filename); // Save the bitmap as an EXR file with the specified filename
	void makeAbs();
	double computeMSE( const Bitmap& groundtruth )const; // computes mean square error

	void saveTXT(const std::string &filename) const;
};

void compositeImages( std::vector<Bitmap::Ptr>& images, const std::string& filename, int maxCols = 6, int spacing = 5 );
Bitmap::Ptr compositeImages2( std::vector<Bitmap::Ptr>& images, int maxCols = 6, int spacing = 5 );

void rasterizeDot( Bitmap& map, const V2d& center, double radius, Color3f value );
void putBitmap( Bitmap& dest, Bitmap& src, const V2d& dest_uv_min, const V2d& dest_uv_max );
void expose( Bitmap& map, float exposure );
void flip( Bitmap& map );
Bitmap::Ptr makeLeftRightImage(Bitmap &left, Bitmap &right );
Color3f variance( Bitmap& map );
Color3f mean( Bitmap& map );
Color3f variance(Bitmap &map, Color3f mean);



struct Filter
{
	int size;
	std::function<double(Bitmap&, int, int, int, int)> weight;
	std::function<double(int, int, int, int)> weight2;
};


Bitmap filter( Bitmap& bm, Filter filter );


Filter gaussFilter(double stddev);
// alpha and beta are smoothing parameters per feature (alpha=distance, beta=color similarity)
Filter bilateralFilter(double alpha, double beta);
