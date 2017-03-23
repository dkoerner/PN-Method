#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <memory>
#include <math/vector.h>
#include <math/rng.h>
#include <math/color.h>











namespace moexp
{
namespace _2d
{
	template<typename T>
	using Function = std::function<T (double)>;


	int ipow(int base, int exp);
	int numFourierCoefficients( int order );
	int numTensorComponents( int order );


	std::unique_ptr<std::vector<double>> projectFourier( int order, const Function<double>& func );
	double fourierSum( int order, double* coeffs, double t );
	double fourierBasis(int order, double* coeffs, double t);







	// tensor is basically a multidimensional array
	// what makes it useful, is the iterator which one can query for its tensor
	// component coordinates etc.
	template<typename T>
	struct Tensor
	{
		Tensor(int order=0):
			m_data(numTensorComponents(order)),
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
			V2i code()
			{
				V2i c(0,0);
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
			double weight( const V2d& d )
			{
				double weight = 1.0;
				for( int j=0;j<m_order;++j )
					weight*=d[m_indices[j]];
				return weight;
			}

			void advance()
			{
				const int dimension = 2;
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
			int dimension = 2;
			return Iterator(m_data.data()+numComponents(), m_order, dimension-1);
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
	T contract( Tensor<T>& a, const V2d& d  )
	{
		T result = 0.0;

		// iterate over all components of tensor a
		for( auto it = a.begin(), end = a.end(); it!=end;++it  )
		{
			result += it.weight(d)*it.value();
		}

		return result;
	}

	std::vector<Tensor<double>> convertFourierCoefficientsToMoments(double* coeffs);
}
}
