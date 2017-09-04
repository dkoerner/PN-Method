#pragma once
#include <iostream>
#include <map>

#include <Eigen/Sparse>

#include <common/Domain.h>
#include <field/Field.h>
#include <field/Constant.h>

// PNSystem is the core class of the pnsolver module. It is mostly concerned with building the system
// Ax=b for PN problems, hence the name.
// 1. The system is initialized with a domain (describing spatial descretization of a rectangular region
//    in worldspace) and an order (effecting the amount of unknowns per voxel).
// 2. Then the setField functions should be used to specify RTE parameters over the computational domain.
// 3. build() will generate A and b
// 4. solve() will return x
//
// NB: this system is for 2d problems.
// NB: the global function set_system_row is being used during the build step to specify the coefficients
//     in A and b for a given voxel (set of rows in the global system). The function is being generated
//     by an external python script directly from the expanded RTE terms.
//
struct PNSystem
{
	typedef std::complex<double> Complex;
	typedef Eigen::SparseMatrix<Complex> ComplexMatrix; // declares a column-major sparse matrix type of double
	typedef Eigen::SparseMatrix<double> RealMatrix; // declares a column-major sparse matrix type of double
	typedef Eigen::SparseMatrix<std::complex<double>> ComplexVector;
	typedef Eigen::SparseMatrix<double> RealVector;
	typedef Eigen::Triplet<std::complex<double>> ComplexTriplet;


	PNSystem( const Domain& domain, int order);

	// This class is used to allow setting coefficients in A and b as if it were dense
	// matrices. The access to the components are stored as triplets which will be
	// used to efficiently build the sparse matrix after all coefficients have been set.
	// The class stores a pointer to a triplet list. Every call to += will generate
	// a new triplet, which is being stored in that list.
	struct MatrixAccessHelper
	{
		MatrixAccessHelper(std::vector<ComplexTriplet>* triplets):
			m_triplets(triplets)
		{

		}

		MatrixAccessHelper& operator+=(const Complex& rhs)
		{
			// add triplet
			if(m_triplets)
				m_triplets->push_back(ComplexTriplet(m_global_i, m_global_j, rhs));
			return *this; // return the result by reference
		}

		int m_global_i;
		int m_global_j;
		std::vector<ComplexTriplet>* m_triplets;
	};




	// This structure holds all (potentially spatially varying) RTE parameters,
	// such as sigma_s, sigma_t etc.
	// For this it uses a generic Field class, which can be a constant, Voxelgrid, etc.
	// Some RTE parameters (such as the 1d phase function and the source term q) will
	// also depend on the spherical harmonics band (l and m) and therefore will have
	// a seperate field for each SH coefficient. This is represented by SHCoefficientFieldArray.
	struct Fields
	{
		Fields( int order ) :
			sigma_t(std::make_shared<Constant>(0.0)),
			sigma_a(std::make_shared<Constant>(0.0)),
			sigma_s(std::make_shared<Constant>(0.0)),
			f_p( std::make_shared<SHCoefficientFieldArray>(order) ),
			q( std::make_shared<SHCoefficientFieldArray>(order) )
		{
		}

		Field::Ptr sigma_t;
		Field::Ptr sigma_a;
		Field::Ptr sigma_s;
		SHCoefficientFieldArray::Ptr f_p;
		SHCoefficientFieldArray::Ptr q;
	};

	// This class serves as an interface to the set_system_row function, which is used to
	// build A and b. The idea is that this function sets the coefficients for all the rows
	// in the global system, which belong to the same voxel. This is because
	// there is a coupling between coefficients within a single voxel and this coupling
	// is expressed by the expanded RTE terms which are used to generate set_system_row.
	// So instead of generating a function, which sets a single row in the global system,
	// it is more intuitive to have a function which sets a complete block-row per voxel.
	struct VoxelSystem
	{
		VoxelSystem(PNSystem* pns,
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


	// This method is used to set RTE parameters. Allowed ids are:
	// sigma_t (Field)
	// sigma_a (Field)
	// sigma_s (Field)
	// f_p  (setting only l=0 and m=0 field of a SHCoefficientFieldArray)
	// q (setting only l=0 and m=0 field of a SHCoefficientFieldArray)
	void setField( const std::string& id, Field::Ptr field );

	// This method builds the matrices A and b and also converts them
	// from complex-valued to real-valued matrices
	void build();

	// this function uses Eigen to solve the system
	RealVector solve();


	const Domain& getDomain()const; // returns the spatial discretization used by the system
	int getGlobalIndex( V2i voxel, int coeff )const; // converts voxel and coefficient into a global index
	void getVoxelAndCoefficient( int global_index, V2i& voxel, int& coeff )const; // converts global index into voxel and coefficient
	int getIndex( int l, int m )const; // returns a linear index from given SH band l, m
	int getNumCoefficients()const; // returns the number of SH coefficients per voxel
	int getNumVoxels()const; // returns the number of voxels


	// the following methods are used to allow python to access the build result
	// (used for storing to disk etc.)
	ComplexMatrix& get_A_complex();
	RealMatrix& get_A_real();

	ComplexVector& get_b_complex();
	RealVector& get_b_real();



	void setDebugVoxel( const V2i& dv ); // this is used for debugging

private:
	// These are used by VoxelSystem during construction of Ax=b
	MatrixAccessHelper A(   V2i voxel_i,
							int coefficient_i,
							V2i voxel_j,
							int coefficient_j );
	MatrixAccessHelper b( V2i voxel_i,
							 int coefficient_i );

	// returns an interface for working with rows in Ax=b which belong to the given voxel
	VoxelSystem getVoxelSystem( const V2i& voxel );

	// builds the matrix S (and its inverse) which is used to convert from complex-valued
	// to real valued matrices (see starmap paper p.5 for details)
	void build_S();


	const Domain m_domain; // defines the spatial discretization
	int m_order;  // PN truncation order
	int m_numCoeffs; // number of coefficients per voxel


	Fields m_fields; // the RTE parameters. Those have to be set by the client code through ::setField

	// the following matrices define the system Ax=b
	ComplexMatrix m_A_complex;
	RealMatrix m_A_real;
	ComplexVector m_b_complex;
	RealVector m_b_real;

	ComplexMatrix m_S; // converts given complex-valued vector to real-valued vector
	ComplexMatrix m_S_inv; // the inverse...


	std::vector<ComplexTriplet> m_triplets_A; // triplets for construction of sparse matrix A (used by MatrixAccessHelper)
	std::vector<ComplexTriplet> m_triplets_b; // triplets for construction of sparse vector b (used by MatrixAccessHelper)
	std::map<std::pair<int, int>, int> m_lm_to_index; // used to convert from l,m to linear index in 2d

	V2i debugVoxel;
};

// This function implments the PN stencil. IT is being generated by a python script directly from the RTE terms.
void set_system_row(PNSystem::VoxelSystem& sys,
					PNSystem::Fields& fields);























