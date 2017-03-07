/*
  This header file is concerned about everything related to moment expansion of spherical functions.
  Since the moment expansion is tightly linked to the spherical harmonics expansion, you will also find
  some spherical harmonics code here.






  Implemented functions and equations most closely follow: http://authors.library.caltech.edu/11159/







*/

#pragma once


#include <vector>
#include <string>
#include <iostream>
#include <memory>
#include <math/vector.h>
#include <math/rng.h>
#include <math/color.h>
#include <math/transform.h>
#include <util/bitmap.h>




namespace sh {

// A spherical function, the first argument is phi, the second is theta.
// See EvalSH(int, int, double, double) for a description of these terms.
typedef std::function<double(double, double)> SphericalFunction;

const int kDefaultSampleCount = 10000;

// Get the total number of coefficients for a function represented by
// all spherical harmonic basis of degree <= @order (it is a point of
// confusion that the order of an SH refers to its degree and not the order).
inline int GetCoefficientCount(int order) {
  return (order + 1) * (order + 1);
}

// Get the one dimensional index associated with a particular degree @l
// and order @m. This is the index that can be used to access the Coeffs
// returned by SHSolver.
inline int GetIndex(int l, int m) {
  return l * (l + 1) + m;
}

// Convert from spherical coordinates to a direction vector. @phi represents
// the rotation about the Z axis and is from [0, 2pi]. @theta represents the
// angle down from the Z axis, from [0, pi].
Eigen::Vector3d ToVector(double phi, double theta);

// Convert from a direction vector to its spherical coordinates. The
// coordinates are written out to @phi and @theta. This is the inverse of
// ToVector.
// Check will fail if @dir is not unit.
void ToSphericalCoords(const Eigen::Vector3d& dir, double* phi, double* theta);

// Convert the (x, y) pixel coordinates into spherical coordinates (phi, theta)
// suitable for use with spherical harmonic evaluation. The x image axis maps
// to phi (0 to 2pi) and the y image axis maps to theta (0 to pi). A pixel index
// maps to the center of the pixel, so phi = 2pi (x + 0.5) / width and
// theta = pi (y + 0.5) / height. This is consistent with ProjectEnvironmentMap.
//
// x and y are not bounds checked against the image, but given the repeated
// nature of trigonometry functions, out-of-bounds x/y values produce reasonable
// phi and theta values (e.g. extending linearly beyond 0, pi, or 2pi).
// Results are undefined if the image dimensions are less than or equal to 0.
//
// The x and y functions are separated because they can be computed
// independently, unlike ToImageCoords.
double ImageXToPhi(int x, int width);
double ImageYToTheta(int y, int height);

// Convert the (phi, theta) spherical coordinates (using the convention
// defined spherical_harmonics.h) to pixel coordinates (x, y). The pixel
// coordinates are floating point to allow for later subsampling within the
// image. This is the inverse of ImageCoordsToSphericalCoords. It properly
// supports angles outside of the standard (0, 2pi) or (0, pi) range by mapping
// them back into it.
Eigen::Vector2d ToImageCoords(double phi, double theta, int width, int height);

// Evaluate the spherical harmonic basis function of degree @l and order @m
// for the given spherical coordinates, @phi and @theta.
// For low values of @l this will use a hard-coded function, otherwise it
// will fallback to EvalSHSlow that uses a recurrence relation to support all l.
double EvalSH(int l, int m, double phi, double theta);

// Evaluate the spherical harmonic basis function of degree @l and order @m
// for the given direction vector, @dir.
// Check will fail if @dir is not unit.
// For low values of @l this will use a hard-coded function, otherwise it
// will fallback to EvalSHSlow that uses a recurrence relation to support all l.
double EvalSH(int l, int m, const Eigen::Vector3d& dir);

// As EvalSH, but always uses the recurrence relationship. This is exposed
// primarily for testing purposes to ensure the hard-coded functions equal the
// recurrence relation version.
double EvalSHSlow(int l, int m, double phi, double theta);

// As EvalSH, but always uses the recurrence relationship. This is exposed
// primarily for testing purposes to ensure the hard-coded functions equal the
// recurrence relation version.
// Check will fail if @dir is not unit.
double EvalSHSlow(int l, int m, const Eigen::Vector3d& dir);

// Fit the given analytical spherical function to the SH basis functions
// up to @order. This uses Monte Carlo sampling to estimate the underlying
// integral. @sample_count determines the number of function evaluations
// performed. @sample_count is rounded to the greatest perfect square that
// is less than or equal to it.
//
// The samples are distributed uniformly over the surface of a sphere. The
// number of samples required to get a reasonable sampling of @func depends on
// the frequencies within that function. Lower frequency will not require as
// many samples. The recommended default kDefaultSampleCount should be
// sufficiently high for most functions, but is also likely overly conservative
// for many applications.
std::unique_ptr<std::vector<double>> ProjectFunction( int order, const SphericalFunction& func, int sample_count);


// Evaluate the already computed coefficients for the SH basis functions up
// to @order, at the coordinates @phi and @theta. The length of the @coeffs
// vector must be equal to GetCoefficientCount(order).
// There are explicit instantiations for double, float, and Eigen::Array3f.
template <typename T>
T EvalSHSum(int order, const std::vector<T>& coeffs, double phi, double theta);

// As EvalSHSum, but inputting a direction vector instead of spherical coords.
// Check will fail if @dir is not unit.
template <typename T>
T EvalSHSum(int order, const std::vector<T>& coeffs, const Eigen::Vector3d& dir);


}  // namespace sh



// complex Y <- complex spherical harmonics
// P <- legrende polynomial
// std::vector<Tensor<complex>> Y_tensors(int l, int m);

namespace moexp
{
	using complex = std::complex<double>;

	template<typename T>
	using SphericalFunction = std::function<T (double, double)>;

	// we start by giving some basic utility functions which are used all over the place

	// pow of integer with non-negative exponent (http://stackoverflow.com/a/101613)
	int ipow(int base, int exp);

	// returns the number of components for a tensor of given order/rank
	int numTensorComponents( int order, int dimension=3 );

	// ...
	double C(int l, int m);
	double a( int l, int m, int j );
	double P(int l, int m, double x); // associated Legendre polynomial
	complex Y( int l, int m, double theta, double phi ); // complex spherical harmonics
	double Y_real(int l, int m, double theta, double phi); // real spherical harmonics


	// tensor is basically a multidimensional array with dynamic number of axes
	// what makes it useful, is the iterator which one can query for its tensor
	// component coordinates etc.
	template<typename T>
	struct Tensor
	{
		Tensor(int order=0):
			m_data(numTensorComponents(order, 3)),
			m_order(order)
		{
		}

		Tensor( const Tensor<T>& other ):
			m_data(other.m_data.begin(), other.m_data.end()),
			m_order(other.m_order)
		{
		}

		// returns rank of the tensor, for naming consistency we call it order
		int getOrder()const
		{
			return m_order;
		}

		T& componentRef( int index )
		{
			return m_data[index];
		}

		int numComponents()
		{
			return m_data.size();
		}

		void multiplyAdd( const T& coefficient, const Tensor<T>& tensor )
		{
			if( m_order != tensor.getOrder() )
			{
				std::cout << "Tensor::multiplyAdd: error: tensors of different rank\n";
				throw std::runtime_error("sfjfrez");
			}
			const T* other_ptr = tensor.getData();
			T* ptr = m_data.data();
			T* end = ptr+m_data.size();
			for( ;ptr!=end;++ptr, ++other_ptr )
			{
				*ptr += *other_ptr*coefficient;
			}
//			for( T* ptr = m_data.data(), end=ptr+m_data.size();ptr!=end;++ptr, ++other_ptr )
//				*ptr += *other_ptr*coefficient;
		}

		struct Iterator
		{
			Iterator(T* ptr, int order, int default_index = 0):
				m_ptr(ptr),
				m_order(order),
				m_indices(order, default_index)
			{
			}

			// returns value of the index with the specified index
			int index( int index_index )
			{
				return m_indices[index_index];
			}

			//
			V3i code()
			{
				V3i c(0,0,0);
				for( int j=0;j<m_order;++j )
					++c[m_indices[j]];
				return c;
			}


			// returns pointer to the current component
			T* ptr()
			{
				return m_ptr;
			}

			// returns value of the current component
			T& value()
			{
				return *m_ptr;
			}

			// this computes the value of the tensor d_i d_j d_k...d_order associated
			// with the current tensor component the iterator points at.
			// it is a shortcut used when contracting moment tensor with direction tensor...
			double weight( const V3d& d )
			{
				double weight = 1.0;
				for( int j=0;j<m_order;++j )
					weight*=d[m_indices[j]];
				return weight;
			}

			void advance()
			{
				const int dimension = 3;
				++m_ptr;

				int increment = 1;
				for( int i=0;i<m_order;++i )
				{
					m_indices[i] += increment;

					// check if index is smaller than dimension
					if(m_indices[i] < dimension)
						return;

					// otherwise we need to handle carry
					increment = m_indices[i]/dimension;
					m_indices[i] = m_indices[i] % dimension;
				}

				return;
			}

			Iterator& operator ++ ()
			{
				advance();
				return *this;
			}
			Iterator operator ++ (int)
			{
				Iterator it = *this;
				it.advance();
				return it;
			}

			bool operator==(const Iterator& other)
			{
				return m_ptr == other.m_ptr;
			}

			bool operator!=(const Iterator& other)
			{
				return m_ptr != other.m_ptr;
			}

			void print()
			{
				for( int i=0;i<m_order;++i )
				{
					std::cout << m_indices[i];
				}
				std::cout << std::endl;
			}

			std::string toString()const
			{
				std::string result = "";

				result += "index=";
				for( int i=0;i<m_order;++i )
					result += ::toString(m_indices[i]);
				return result;
			}

			std::string index_str()const
			{
				std::string result = "";
				for( int i=0;i<m_order;++i )
					result += ::toString(m_indices[i]);
				return result;
			}

		private:
			T* m_ptr;
			std::vector<int> m_indices;
			int m_order;
		};

		Iterator begin()
		{
			return Iterator(m_data.data(), m_order);
		}

		Iterator end()
		{
			return Iterator(m_data.data()+numComponents(), m_order, 2);
		}

		const T* getData()const
		{
			return m_data.data();
		}

	private:
		std::vector<T> m_data;
		int m_order;
	};


	template<typename T>
	T contract( Tensor<T>& a, const V3d& d  )
	{
		T result = 0.0;

		// iterate over all components of tensor a
		for( auto it = a.begin(), end = a.end(); it!=end;++it  )
		{
			result += it.weight(d)*it.value();
		}

		return result;
	}

	std::vector<Tensor<complex>> Y_tensors( int l, int m ); // Y in tensor form (required for moment expansion)

	// Get the total number of coefficients for a function represented by
	// all spherical harmonic basis of degree <= @order (it is a point of
	// confusion that the order of an SH refers to its degree and not the order).
	inline int numSHCoefficients(int order)
	{
		return (order + 1) * (order + 1);
	}

	// Get the one dimensional index associated with a particular degree @l
	// and order @m. This is the index that can be used to access the Coeffs
	// returned by SHSolver.
	inline int shIndex(int l, int m)
	{
		return l * (l + 1) + m;
	}

	std::unique_ptr<std::vector<complex>> project_Y(int order, const SphericalFunction<double>& func, int sample_count);

	template<typename T>
	std::unique_ptr<std::vector<T>> project_Y_real(int order, const SphericalFunction<T>& func, int sample_count)
	{
		std::unique_ptr<std::vector<T>> coeffs(new std::vector<T>());
		coeffs->assign(numSHCoefficients(order), 0.0);

		RNGd rng;
		// integrate func times sh basis function to get the coefficient associated with that basis
		for( int i=0;i<sample_count;++i )
		{
			V3d d = sampleSphere<double>(rng);
			P2d theta_phi = sphericalCoordinates(d);
			double theta = theta_phi.x();
			double phi = theta_phi.y();
			double d_pdf = sampleSpherePDF();

			T f = func( theta_phi.x(), theta_phi.y() );

			// integrate and accumulate for each sh basis function
			for( int l=0;l<=order;++l )
				for( int m=-l;m<=l;++m )
				{
					double sh = Y_real(l, m, theta, phi);
					T sample = f*sh/d_pdf;
					T& sh_coeff = (*coeffs)[shIndex(l, m)];
					sh_coeff += (sample - sh_coeff)/float(i+1);
				}
		}
		return coeffs;
	}

	complex Y_sum(int order, const std::vector<complex>& coeffs, double theta, double phi);


	template <typename T>
	T Y_real_sum(int order, const std::vector<T>& coeffs, double theta, double phi)
	{
		//CHECK(GetCoefficientCount(order) == coeffs.size(), "Incorrect number of coefficients provided.");
		T sum = T(0.0);//Zero<T>();
		for (int l = 0; l <= order; l++)
		{
			for (int m = -l; m <= l; m++)
			{
				sum += Y_real(l, m, theta, phi) * coeffs[shIndex(l, m)];
			}
		}
		return sum;
	}

	template <typename T>
	T Y_real_sum(int order, const Color3f* coeffs, double theta, double phi)
	{
		//CHECK(GetCoefficientCount(order) == coeffs.size(), "Incorrect number of coefficients provided.");
		T sum = T(0.0);//Zero<T>();
		for (int l = 0; l <= order; l++)
		{
			for (int m = -l; m <= l; m++)
			{
				sum += Y_real(l, m, theta, phi) * coeffs[shIndex(l, m)];
			}
		}
		return sum;
	}

}//namespace moexp


struct EnvMap
{
	typedef std::shared_ptr<EnvMap> Ptr;

	EnvMap( const std::string& filename )
	{
		m_bitmap = Bitmap(filename);
		m_transform = Transformd();
	}

	EnvMap( const V2i& res = V2i(512, 256) )
	{
		m_bitmap = Bitmap(res);
		m_transform = Transformd();
	}

	// evaluate environment map
	Color3f eval( double theta, double phi )const
	{
		V3d d = sphericalDirection<double>(theta, phi);
		P2d uv = directionToUV(d);
		return m_bitmap.eval(uv);
	}


	P2d directionToUV( const V3d& d )const
	{
		// using formulas given in http://gl.ict.usc.edu/Data/HighResProbes/
		// with the difference that u=[0,1] (instead of [0,2]) and we negate z
		P2d uv( (1+std::atan2(d.x(), d.z())/M_PI)/2,
				 safe_acos(d.y())/M_PI );
		return uv;
	}
	V3d uvToDirection( const P2d& uv )const
	{
		// using formulas given in http://gl.ict.usc.edu/Data/HighResProbes/
		// with the difference that u=[0,1] (instead of [0,2]) and we negate z
		// azimuthal angle
		double theta = M_PI*(uv.x()*2.0-1.0);
		// elevation angle
		double phi = M_PI*uv.y();
		return V3d( std::sin(phi)*std::sin(theta), std::cos(phi), std::sin(phi)*cos(theta) );
	}
	P2d uvToXY(const P2d& uv)const
	{
		P2d xy(
			(uv.x()*(m_bitmap.cols()-1)),
			(uv.y()*(m_bitmap.rows()-1))
			);
		return xy;
	}

	P2d xyToUV(const P2d& xy)const
	{
		return P2d(
			(xy.x())/double(m_bitmap.cols()-1),
			(xy.y())/double(m_bitmap.rows()-1)
			);
	}

	V3d xyToDirection( const P2d& xy )const
	{
		return uvToDirection( xyToUV(xy) );
	}
	P2d directionToXY( const V3d& d )const
	{
		return uvToXY(directionToUV(d));
	}

	Bitmap& bitmap()
	{
		return m_bitmap;
	}

	void saveGeo( const std::string& filename, double exposure = 0.0);
//private:
	Transformd m_transform;
	Bitmap m_bitmap;
};


void rasterizeSphericalFunctionMap(const std::string& filename, moexp::SphericalFunction<Color3f> func);
void rasterizeSphericalFunctionSphere(const std::string& filename, moexp::SphericalFunction<Color3f> func, double exposure =0.0);

