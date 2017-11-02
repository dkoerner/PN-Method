#pragma once
#include <string>
#include <memory>
#include <math/vector.h>








struct Image// : public Eigen::Array<Eigen::Vector3d, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
{
public:
	typedef std::shared_ptr<Image> Ptr;
	typedef Eigen::Array<Eigen::Vector3d, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Base;

	// construction
	Image(const Vector2i &size = Vector2i(0, 0));
	Image(const Image& image);

	// io
	const Eigen::Vector3d* data()const;
	Base& getArray();
	Eigen::Vector3d& pixel( int x, int y );
	void save(const std::string& filename); // stores content as exr file
	std::string toString()const;

	// info
	V2i getResolution()const;

private:
	Base m_data;
};
