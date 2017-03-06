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
#include <math/RNG.h>


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

}//namespace moexp
