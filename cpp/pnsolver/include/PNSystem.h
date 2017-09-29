#pragma once
#include <iostream>
#include <map>

#include <Eigen/Sparse>

#include <common/Domain.h>
#include <field/Field.h>
#include <field/Constant.h>


#define REGISTER_STENCIL(name, order, width) \
	static struct name ##_{ \
		name ##_() { \
			PNSystem::registerStencil(#name, PNSystem::Stencil( order, width, name )); \
		} \
	}name ##__;


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
	typedef Eigen::Triplet<double> RealTriplet;
	struct VoxelSystem;
	struct Fields;

	struct Stencil
	{
		typedef std::function<void(PNSystem&, const V2i&)> StencilFunction;

		Stencil( int _order = -1, int _width = 0, StencilFunction fun = StencilFunction() ):
			order(_order),
			width(_width),
			apply(fun)
		{

		}

		StencilFunction apply;
		int order;

		// this tells us how many neighbouring voxels the stencil needs to see
		// this value is required for setting up boundary voxels
		// widht=0 means: no neighbours are touched
		int width;

		static Stencil noop()
		{
			return Stencil(0, 2, [](PNSystem& sys, const V2i& voxel){});
		}
	};


	PNSystem(Stencil stencil, const Domain& domain);



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

	template<typename T>
	struct MatrixBuilder
	{
		typedef Eigen::SparseMatrix<T> Matrix;
		typedef Eigen::Triplet<T> Triplet;


		// This class is used to allow setting coefficients in A and b as if it were dense
		// matrices. The access to the components are stored as triplets which will be
		// used to efficiently build the sparse matrix after all coefficients have been set.
		// The class stores a pointer to a triplet list. Every call to += will generate
		// a new triplet, which is being stored in that list.
		struct MatrixAccessHelper
		{
			MatrixAccessHelper(std::vector<Triplet>* triplets = 0):
				triplets(triplets)
			{

			}

			MatrixAccessHelper& operator+=(const T& rhs)
			{
				// add triplet
				if(triplets)
					triplets->push_back(Triplet(i, j, rhs));
				return *this; // return the result by reference
			}

			int i;
			int j;
			std::vector<Triplet>* triplets;
		};

		MatrixAccessHelper coeff(int i,int j)
		{
			MatrixAccessHelper mah(&triplets);
			mah.i = i;
			mah.j = j;
			return mah;
		}

		void reset()
		{
			triplets.clear();
			matrix = Matrix();
		}

		void build( int numRows, int numCols )
		{
			matrix = Matrix(numRows, numCols);

			// build sparse matrix from triplets
			matrix.setFromTriplets(triplets.begin(), triplets.end());
		}

		void add( const MatrixBuilder<T>& other )
		{
			for( auto t:other.triplets )
				triplets.push_back(t);
		}

		std::vector<Triplet> triplets;
		V2i                  dimensions;
		Matrix               matrix;
	};
	typedef MatrixBuilder<double> MatrixBuilderd;



	// This structure holds all (potentially spatially varying) RTE parameters,
	// such as sigma_s, sigma_t etc.
	// For this it uses a generic Field class, which can be a constant, Voxelgrid, etc.
	// Some RTE parameters (such as the 1d phase function and the source term q) will
	// also depend on the spherical harmonics band (l and m) and therefore will have
	// a seperate field for each SH coefficient. This is represented by SHCoefficientFieldArray.
	struct Fields
	{
		//Fields()
		//{
		//}

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


		MatrixBuilderd::MatrixAccessHelper coeff_A(int coefficient_i, V2i voxel_j, int coefficient_j);
		MatrixBuilderd::MatrixAccessHelper coeff_b(int coefficient_i);

		V2d voxelToWorld(const V2d& pVS)const;
		const V2i& getVoxel()const;
		P2d getVoxelSize()const;
	private:
		V2i m_voxel_i;
		PNSystem* m_pns;
	};

	struct Voxel
	{
		Voxel():type(-1),tmp0(-100),tmp1(-100)
		{

		}

		enum EVoxelType
		{
			EVT_Interior=0, // contains no boundary coefficients
			EVT_E=1, // below are all voxels containing boundary and interior coefficients
			EVT_W=2,
			EVT_N=3,
			EVT_S=4,
			EVT_NE=5,
			EVT_NW=6,
			EVT_SW=7,
			EVT_SE=8,
			EVT_Boundary=9 // contains only boundary coefficients
		};

		V2i coord;
		int globalOffset; // points at the global index of the first coefficient
		int type; // identifies the voxel type (interior, east-boundary, etc.)
		int tmp0;
		int tmp1;
	};
	std::vector<Voxel> m_voxels;

	struct VoxelManager
	{
		VoxelManager();
		void init(PNSystem* sys);
		Voxel createVoxel( const V2i& coord, int globalOffset );
		int getNumCoeffs(const Voxel& voxel);
		int getGlobalIndex( const Voxel& voxel, int coeff );
	private:
		PNSystem* sys;
		std::map<std::pair<int,int>, int> mixedTypes;

		// maps [type][coeff_index] to local index within voxel
		// this mapping is different per type
		std::vector<std::vector<int>> m_localCoefficientIndices;
	};

	MatrixBuilderd::MatrixAccessHelper coeff_A(V2i voxel_i, int coefficient_i, V2i voxel_j, int coefficient_j);
	MatrixBuilderd::MatrixAccessHelper coeff_b(V2i voxel_i, int coefficient_i);



	// computes the groundtruth by doing 2d path tracing
	Eigen::MatrixXd computeGroundtruth( int numSamples );


	Fields& getFields();

	// This method is used to set RTE parameters. Allowed ids are:
	// sigma_t (Field)
	// sigma_a (Field)
	// sigma_s (Field)
	// f_p  (setting only l=0 and m=0 field of a SHCoefficientFieldArray)
	// q (setting only l=0 and m=0 field of a SHCoefficientFieldArray)
	void setField( const std::string& id, Field::Ptr field );

	void setNeumannBoundaryConditions( bool flag );


	// This method builds the matrices A and b and also converts them
	// from complex-valued to real-valued matrices
	void build();

	RealMatrix getIndexMatrix();
	RealMatrix getVoxelInfo( const std::string& info );

	// this function uses Eigen to solve the system
	RealVector solve();
	RealVector solveWithGuess(RealVector &x0 );
	//ComplexVector solve2();
	RealVector solve_boundary();


	const Domain& getDomain()const; // returns the spatial discretization used by the system
	int getGlobalIndex( V2i voxel, int coeff )const; // converts voxel and coefficient into a global index
	void getVoxelAndCoefficient( int global_index, V2i& voxel, int& coeff )const; // converts global index into voxel and coefficient
	int getIndex( int l, int m )const; // returns a linear index from given SH band l, m
	void getLM( int sh_index, int& l, int& m )const; // returns the sh band l,m from linear index
	int getNumCoefficients()const; // returns the number of SH coefficients per voxel
	int getNumVoxels()const; // returns the number of voxels
	int getOrder()const;


	// the following methods are used to allow python to access the build result
	// (used for storing to disk etc.)
	RealMatrix& get_A_real();
	RealVector& get_b_real();
	//ComplexMatrix& get_A_complex();
	//ComplexVector& get_b_complex();




	void setDebugVoxel( const V2i& dv ); // this is used for debugging

	// temp boundary test
	RealMatrix& get_boundary_A_real();
	RealMatrix& get_boundary_b_real();

	static void registerStencil( const std::string& name, Stencil stencil )
	{
		std::cout << "PNSystem::registerStencil: registering stencil \"" << name << "\" with order=" << stencil.order << std::endl;
		g_stencils[name] = stencil;
	}
	static Stencil findStencil( const std::string& id )
	{
		// find stencil
		auto it = g_stencils.find(id);
		if( it == g_stencils.end() )
			throw std::runtime_error( "PNSystem::findStencil: unable to find stencil " + id );
		return it->second;
	}


public:
	PNSystem::Voxel &getVoxel2(const V2i &voxel);
	int getVoxel3(const V2i &voxel);
	PNSystem::Voxel createVoxel(const V2i &coord, int globalOffset);
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



	std::map<std::pair<int, int>, int> m_lm_to_index; // used to convert from l,m to linear index in 2d
	std::map<int, std::pair<int, int>> m_index_to_lm; // used to convert from linear index to l,m in 2d

	V2i debugVoxel;


	// TEMP: testing neumann BC hack
	// the following matrices define the system Ax=b
	int getGlobalBoundaryIndex( const V2i boundaryVoxel, int coeff );
	int m_numBoundaryVoxels;
	std::map<std::pair<int,int>, int> m_voxel_to_global_boundary_index;

	MatrixBuilderd m_builder_A_boundary;
	MatrixBuilderd m_builder_b_boundary;

	/*
	ComplexMatrix m_boundary_A_complex;
	RealMatrix m_boundary_A_real;
	ComplexVector m_boundary_b_complex;
	RealVector m_boundary_b_real;
	ComplexMatrix m_boundary_S; // converts given complex-valued vector to real-valued vector
	ComplexMatrix m_boundary_S_inv; // the inverse...
	*/


	// computational domain, problem and stencil
	int m_numCoeffs; // number of coefficients per voxel
	int m_numBoundaryLayers;
	const Domain m_domain; // defines the spatial discretization
	Fields m_fields; // the RTE parameters. Those have to be set by the client code through ::setField
	bool m_neumannBC;
	Stencil m_stencil;


	// the following matrices define the system Ax=b
	MatrixBuilderd m_builder_A;
	MatrixBuilderd m_builder_b;

	VoxelManager m_voxelManager;
	int m_numSystemEquations; // number of cols and rows of A


	static std::map<std::string, Stencil> g_stencils; // global register of stencils
};
























