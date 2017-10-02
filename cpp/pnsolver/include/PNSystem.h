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
			PNSystem::registerStencil(#name, PNSystem::Stencil( order, width, name, name##_get_offset )); \
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
	struct Voxel
	{
		Voxel():
			type(-1)
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

		static std::string typeStr( int type )
		{
			switch(type)
			{
			case PNSystem::Voxel::EVT_Interior:return "Interior";
			case PNSystem::Voxel::EVT_Boundary:return "Boundary";
			case PNSystem::Voxel::EVT_E:return "E";
			case PNSystem::Voxel::EVT_S:return "S";
			case PNSystem::Voxel::EVT_W:return "W";
			case PNSystem::Voxel::EVT_N:return "N";
			case PNSystem::Voxel::EVT_NE:return "NE";
			case PNSystem::Voxel::EVT_SE:return "SE";
			case PNSystem::Voxel::EVT_SW:return "SW";
			case PNSystem::Voxel::EVT_NW:return "NW";
			};
		}

		V2i coord;
		int globalOffset; // points at the global index of the first coefficient
		int type; // identifies the voxel type (interior, east-boundary, etc.)
	};

	struct VoxelManager
	{
		VoxelManager();
		void init(PNSystem* sys);
		int getNumCoeffs(const Voxel& voxel);
		int getGlobalIndex( const Voxel& voxel, int coeff );
		bool isBoundaryCoefficient(const Voxel &voxel, int coeff ); // returns true, if the given coefficient is a boundary coefficient

		Voxel &getVoxel(const V2i &voxel);
		bool voxelIsValid( const P2i& voxel )const
		{
			if( (voxel[0] < -m_numBoundaryLayers)||(voxel[0] >= m_resolution[0]+m_numBoundaryLayers)||
				(voxel[1] < -m_numBoundaryLayers)||(voxel[1] >= m_resolution[1]+m_numBoundaryLayers))
				return false;
			return true;
		}
		std::vector<Voxel>& getVoxels()
		{
			return m_voxels;
		}

		int getNumUnknowns()const
		{
			return m_numUnknowns;
		}

		V2i getVoxelMin()
		{
			return V2i(-m_numBoundaryLayers, -m_numBoundaryLayers);
		}

		V2i getVoxelMax()
		{
			return V2i(m_resolution[0]+m_numBoundaryLayers, m_resolution[1]+m_numBoundaryLayers);
		}
	private:

		// returns boundary layer indices of the current voxel in x and y
		// the indices are signed, indicating left or right side of the domain.
		// if no layer is touched: zero is returned
		V2i getBoundaryLayer( const Voxel& voxel );

		PNSystem* sys;
		std::map<std::pair<int,int>, int> mixedTypes;


		int m_numUnknowns; // number of cols and rows of A
		int m_numBoundaryLayers;
		V2i m_resolution;


		// maps [type][coeff_index] to local index within voxel
		// this mapping is different per type
		std::vector<std::vector<int>> m_localCoefficientIndices;
		std::vector<std::vector<bool>> m_isBoundaryCoefficient; // this table flags each boundary coefficient
		std::vector<int> m_numNonBoundaryCoefficients; // stores the number of non-boundary coefficients for every voxeltype
		std::vector<Voxel> m_voxels;
	};







	struct Stencil
	{
		struct Context
		{
			Context( PNSystem& sys, Voxel& voxel ):
				sys(sys),
				voxel(voxel)
			{
				if(!sys.m_voxelManager.voxelIsValid(voxel.coord))
					throw std::runtime_error("Stencil::Context::Context: voxel is expected to be valid.");
			}

			V2i getVoxel()
			{
				return voxel.coord;
			}

			const Domain& getDomain()
			{
				return sys.getDomain();
			}

			Fields& getFields()
			{
				return sys.getFields();
			}

			MatrixBuilderd::MatrixAccessHelper coeff_A( int coeff_i, const V2i& voxel_j, int coeff_j );
			MatrixBuilderd::MatrixAccessHelper coeff_b( int coeff_i );

		private:
			PNSystem& sys;
			Voxel& voxel;
		};

		typedef std::function<void(Context&)> StencilFunction;
		typedef std::function<V2i(int)> GetCoefficientOffsetFunction;

		Stencil( int _order = -1,
				 int _width = 0,
				 StencilFunction fun = StencilFunction(),
				 GetCoefficientOffsetFunction getOffset_fun = GetCoefficientOffsetFunction() ):
			order(_order),
			width(_width),
			apply(fun),
			getOffset(getOffset_fun)
		{

		}

		StencilFunction apply;
		GetCoefficientOffsetFunction getOffset;
		int order;

		// this tells us how many neighbouring voxels the stencil needs to see
		// this value is required for setting up boundary voxels
		// widht=0 means: no neighbours are touched
		int width;

		static Stencil noop()
		{
			return Stencil(0, 2,
						   [](Context& sys){},
						   [](int coeff){return V2i(1,1);}); // noop places all coefficients at the cell center
		}
	};

	PNSystem(Stencil stencil, const Domain& domain, bool neumannBC);

	MatrixBuilderd::MatrixAccessHelper coeff_A(V2i voxel_i, int coefficient_i, V2i voxel_j, int coefficient_j);
	MatrixBuilderd::MatrixAccessHelper coeff_b(V2i voxel_i, int coefficient_i);



	// computes the groundtruth by doing 2d path tracing
	Eigen::MatrixXd computeGroundtruth( int numSamples );


	Stencil& getStencil();

	Fields& getFields();
	// This method is used to set RTE parameters. Allowed ids are:
	// sigma_t (Field)
	// sigma_a (Field)
	// sigma_s (Field)
	// f_p  (setting only l=0 and m=0 field of a SHCoefficientFieldArray)
	// q (setting only l=0 and m=0 field of a SHCoefficientFieldArray)
	void setField( const std::string& id, Field::Ptr field );

	int getBoundaryConditions();


	// This method builds the matrices A and b and also converts them
	// from complex-valued to real-valued matrices
	void build();

	RealMatrix getVoxelInfo( const std::string& info );

	// this function uses Eigen to solve the system
	RealVector solve();
	RealVector solveWithGuess(RealVector &x0 );
	//ComplexVector solve2();


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


	VoxelManager& getVoxelManager();

private:

	std::map<std::pair<int, int>, int> m_lm_to_index; // used to convert from l,m to linear index in 2d
	std::map<int, std::pair<int, int>> m_index_to_lm; // used to convert from linear index to l,m in 2d

	V2i debugVoxel;


	// TEMP: testing neumann BC hack
	// the following matrices define the system Ax=b
	int getGlobalBoundaryIndex( const V2i boundaryVoxel, int coeff );
	int m_numBoundaryVoxels;
	std::map<std::pair<int,int>, int> m_voxel_to_global_boundary_index;



	// computational domain, problem and stencil
	int m_numCoeffs; // number of coefficients per voxel
	const Domain m_domain; // defines the spatial discretization
	Fields m_fields; // the RTE parameters. Those have to be set by the client code through ::setField
	bool m_neumannBC;
	Stencil m_stencil;
	bool m_debug;


	// the following matrices define the system Ax=b
	MatrixBuilderd m_builder_A;
	MatrixBuilderd m_builder_b;

	VoxelManager m_voxelManager;

	static std::map<std::string, Stencil> g_stencils; // global register of stencils
};
























