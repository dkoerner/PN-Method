#pragma once
#include <iostream>
#include <map>

#include <Eigen/Sparse>

#include <common/Domain.h>
#include <field/Field.h>
#include <field/Constant.h>


struct PNSystem
{
	typedef std::complex<double> Complex;
	typedef Eigen::SparseMatrix<Complex> ComplexMatrix; // declares a column-major sparse matrix type of double
	typedef Eigen::SparseMatrix<double> RealMatrix; // declares a column-major sparse matrix type of double
	typedef Eigen::SparseMatrix<std::complex<double>> ComplexVector;
	typedef Eigen::SparseMatrix<double> RealVector;
	typedef Eigen::Triplet<std::complex<double>> ComplexTriplet;


	PNSystem( const Domain& domain,
			  int order);

	struct MatrixAccessHelper
	{
		MatrixAccessHelper(std::vector<ComplexTriplet>* triplets):
			m_triplets(triplets)
		{

		}

		MatrixAccessHelper& operator+=(const Complex& rhs)
		{
			// add triplet
			//std::cout << "setting matrix at i=" << m_global_i << "j=" << m_global_j << " to " << rhs << std::endl;
			if(m_triplets)
				m_triplets->push_back(ComplexTriplet(m_global_i, m_global_j, rhs));
			return *this; // return the result by reference
		}

		int m_global_i;
		int m_global_j;
		std::vector<ComplexTriplet>* m_triplets;
	};


	struct SHCoefficientArray
	{
		typedef std::shared_ptr<SHCoefficientArray> Ptr;

		SHCoefficientArray( int order );

		void setField( int l, int m, Field::Ptr field );
		std::complex<double> eval( int l, int m, const V2d& pWS )const;

	private:
		std::vector<Field::Ptr> m_coeff_fields; // a field for each coefficient
	};

	struct Fields
	{
		Fields( int order ) :
			sigma_t(std::make_shared<Constant>(0.0)),
			sigma_a(std::make_shared<Constant>(0.0)),
			sigma_s(std::make_shared<Constant>(0.0)),
			f_p( std::make_shared<SHCoefficientArray>(order) ),
			q( std::make_shared<SHCoefficientArray>(order) )
		{
		}

		Field::Ptr sigma_t;
		Field::Ptr sigma_a;
		Field::Ptr sigma_s;
		SHCoefficientArray::Ptr f_p;
		SHCoefficientArray::Ptr q;
	};

	struct Voxel
	{
		Voxel(PNSystem* pns,
			  const V2i& voxel_i );

		MatrixAccessHelper A(   int coefficient_i,
								V2i voxel_j,
								int coefficient_j );
		MatrixAccessHelper b( int coefficient_i );

		V2d voxelToWorld(const V2d& pVS)const;
		const V2i& getVoxel()const;
	private:
		V2i m_voxel_i;
		PNSystem* m_pns;
	};


	void setDebugVoxel( const V2i& dv )
	{
		debugVoxel = dv;
	}

	Voxel getVoxel( const V2i& voxel );
	const Domain& getDomain()const;
	int getGlobalIndex( V2i voxel, int coeff )const;
	void getVoxelAndCoefficient( int global_index, V2i& voxel, int& coeff )const;
	int getIndex( int l, int m )const;
	int getNumCoefficients()const;
	int getNumVoxels()const;


	MatrixAccessHelper A(   V2i voxel_i,
							int coefficient_i,
							V2i voxel_j,
							int coefficient_j );
	MatrixAccessHelper b( V2i voxel_i,
							 int coefficient_i );

	void setField( const std::string& id, Field::Ptr field );
	void build();
	ComplexMatrix& get_A();
	ComplexVector& get_b();
	RealMatrix& get_A_real();
	RealVector& get_b_real();

	ComplexMatrix& get_S()
	{
		return m_debug_S;
	}

	ComplexVector& get_S_inv()
	{
		return m_debug_S_inv;
	}

	RealMatrix& get_M()
	{
		return m_debug_M;
	}

private:
	void build_S();

	const Domain m_domain;
	int m_order;
	int m_numCoeffs;
	Fields m_fields;
	ComplexMatrix m_A_complex;
	RealMatrix m_A_real;
	ComplexVector m_b_complex;
	RealVector m_b_real;
	ComplexMatrix m_S;
	ComplexMatrix m_S_inv;

	ComplexMatrix m_debug_S;
	ComplexMatrix m_debug_S_inv;
	RealMatrix m_debug_M;


	std::vector<ComplexTriplet> m_triplets_A;
	std::vector<ComplexTriplet> m_triplets_b;

	std::map<std::pair<int, int>, int> m_lm_to_index;

	V2i debugVoxel;
};


void set_system_row(PNSystem::Voxel& sys,
					PNSystem::Fields& fields);























