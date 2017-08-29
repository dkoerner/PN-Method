#pragma once
#include <iostream>
#include <map>

#include <Eigen/Sparse>

#include <common/Domain.h>
#include <field/Field.h>
#include <field/Constant.h>


struct PNSystem
{
	typedef std::complex<double> ValueType;
	typedef Eigen::SparseMatrix<ValueType> Matrix; // declares a column-major sparse matrix type of double
	typedef Eigen::Triplet<ValueType> Triplet;


	PNSystem( const Domain& domain,
			  int order);

	struct MatrixAccessHelper
	{
		MatrixAccessHelper(PNSystem* pns):
			m_pns(pns)
		{

		}

		MatrixAccessHelper& operator+=(const ValueType& rhs)
		{
			// add triplet
			//std::cout << "setting matrix at i=" << m_global_i << "j=" << m_global_j << " to " << rhs << std::endl;
			m_pns->m_triplets.push_back(Triplet(m_global_i, m_global_j, rhs));
			return *this; // return the result by reference
		}

		int m_global_i;
		int m_global_j;
		PNSystem* m_pns;
	};

	struct Fields
	{
		Fields()
		{
			sigma_t = std::make_shared<Constant>(0.0);
			sigma_a = std::make_shared<Constant>(0.0);
			sigma_s = std::make_shared<Constant>(0.0);
		}

		Field::Ptr sigma_t;
		Field::Ptr sigma_a;
		Field::Ptr sigma_s;
	};

	struct Voxel
	{
		Voxel(PNSystem* pns,
			  const V2i& voxel_i );

		MatrixAccessHelper A(   int coefficient_i,
								V2i voxel_j,
								int coefficient_j );
		//std::complex<double>& b( int coefficient_i );

		V2d voxelToWorld(const V2d& pVS)const;
		const V2i& getVoxel()const;
	private:
		V2i m_voxel_i;
		PNSystem* m_pns;
	};



	Voxel getVoxel( const V2i& voxel );
	const Domain& getDomain()const;
	int getGlobalIndex( V2i voxel, int coeff )const;
	int getNumCoefficients()const;
	int getNumVoxels()const;


	MatrixAccessHelper A(   V2i voxel_i,
							int coefficient_i,
							V2i voxel_j,
							int coefficient_j );
	//std::complex<double>& b( V2i voxel_i,
	//						 int coefficient_i );

	void setField( const std::string& id, Field::Ptr field );
	void build();
	Matrix& get_A();
private:
	const Domain m_domain;
	int m_order;
	int m_numCoeffs;
	Fields m_fields;
	Matrix m_matrix;

	std::vector<Triplet> m_triplets;
};


void set_system_row(PNSystem::Voxel& sys,
					PNSystem::Fields& fields);























