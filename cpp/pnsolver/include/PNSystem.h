#pragma once
#include <iostream>
#include <map>

#include <Eigen/Sparse>

#include <common/Domain.h>
#include <field/Field.h>
#include <field/Constant.h>


#define REGISTER_STENCIL(name, order, numCoeffs, width) \
	static struct name ##_{ \
		name ##_() { \
			PNSystem::registerStencil(#name, PNSystem::Stencil( order, numCoeffs, width, name, name##_get_offset )); \
		} \
	}name ##__;

// PNSystem is the core class of the pnsolver module. It is mostly concerned with building the system
// Ax=b for PN problems, hence the name.
// 1. The system is initialized with a domain (describing spatial descretization of a rectangular region
//    in worldspace) and a stencil (which effects the number of unknowns)
// 2. Then the setField functions should be used to specify RTE parameters over the computational domain.
// 3. build() will generate A and b
// 4. solve() will return x
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
			MatrixAccessHelper(std::vector<Triplet>* triplets = 0, bool debug = false):
				triplets(triplets),
				debug(debug)
			{

			}

			MatrixAccessHelper& operator+=(const T& rhs)
			{
				// add triplet
				if(triplets)
				{
					if(debug)
						std::cout << "MatrixAccessHelper: adding triplet i=" << i << " j=" << j << " value=" << rhs << std::endl;
					triplets->push_back(Triplet(i, j, rhs));
				}
				return *this; // return the result by reference
			}

			int i;
			int j;
			bool debug;
			std::vector<Triplet>* triplets;
		};

		MatrixAccessHelper coeff(int i,int j, bool debug = false)
		{
			MatrixAccessHelper mah(&triplets, debug);
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

	// This class serves as an interface to the stencil, which is used to
	// build A and b. The idea is that the stencils apply function sets the coefficients for all the rows
	// in the global system, which belong to the same voxel. This is because
	// there is a coupling between coefficients within a single voxel and this coupling
	// is expressed by the expanded RTE terms which are used to generate the stencil.
	// So instead of generating a function, which sets a single row in the global system,
	// it is more intuitive to have a function which sets a complete block-row per voxel.
	struct Voxel
	{
		Voxel():
			type(-1),
			debug(false)
		{
		}

		V3i coord;
		int globalOffset; // points at the global index of the first coefficient
		int type; // identifies the voxel type (interior, east-boundary, etc.)
		bool debug;
	};

	struct VoxelManager
	{
		struct VoxelType
		{
			VoxelType(int index = -1, int maxNumCoefficients = 0, const std::tuple<int, int, int> &boundaryLayer = std::make_tuple(0,0,0) );

			int getIndex()const;

			std::tuple<int, int, int>& getBoundaryLayer();

			// This method is called during voxelmanager initialization when the voxeltypes are created.
			// It changes the local coefficient indices and flags accordingly.
			void registerActiveCoefficient( int coeff_index );
			int getNumActiveCoefficients()const;
			// returns true if the given coefficient is active and therefore
			// has an index within the global system
			bool isActiveCoefficient(int coeff_index)const;

			// returns the coefficient index within a single voxel (local) for the specific voxel type
			// returns -1 if the coefficient is not active and therefore has no index within the global system
			int getLocalCoefficientIndex( int coeff_index )const;

		private:
			int m_index;

			// this identifies the location of the voxel in respect to the domain for which we want to solve
			// 0,0,0 means, that the voxel is an internal voxel and lies somewhere within
			// 1,0,0 means that the voxel sits on the first layer of voxels north to the domain
			// 2,0,0 means that the voxel sits on the second layer of voxels north to the domain etc.
			std::tuple<int, int, int> m_boundaryLayer;

			// this maps coefficient indices to a local coefficient index the reason is,
			// that for some voxeltypes (on the boundary), not all coefficients have a
			// row in the global system. Therefore voxels of these types skip certain
			// coefficients which requires a remapping
			std::vector<int> m_localCoefficientIndices;
			// this stores a flag for every coefficient, saying if the coefficient is a boundary coefficient
			// TODO: this is not necessary anymore. if the local coefficient is -1, we can conclude that we
			// have a boundary coefficient
			std::vector<bool> m_isBoundaryCoefficient;

			// the number of active coefficients (coefficients with local index >=0)
			int m_numNonBoundaryCoefficients;
		};


		VoxelManager();
		void init(PNSystem* sys);
		int getNumCoeffs(const Voxel& voxel);
		int getGlobalIndex( const Voxel& voxel, int coeff );
		bool isBoundaryCoefficient(const Voxel &voxel, int coeff ); // returns true, if the given coefficient is a boundary coefficient
		V3i getNumBoundaryLayers()const;

		Voxel &getVoxel(const V3i &coord);
		bool voxelIsValid( const P3i& voxel )const
		{
			if( (voxel[0] < -m_numBoundaryLayers[0])||(voxel[0] >= m_resolution[0]+m_numBoundaryLayers[0])||
				(voxel[1] < -m_numBoundaryLayers[1])||(voxel[1] >= m_resolution[1]+m_numBoundaryLayers[1])||
				(voxel[2] < -m_numBoundaryLayers[2])||(voxel[2] >= m_resolution[2]+m_numBoundaryLayers[2]))
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

		V3i getVoxelMin()
		{
			return V3i(-m_numBoundaryLayers[0], -m_numBoundaryLayers[1], -m_numBoundaryLayers[2]);
		}

		V3i getVoxelMax()
		{
			return V3i(m_resolution[0]+m_numBoundaryLayers[0], m_resolution[1]+m_numBoundaryLayers[1], m_resolution[2]+m_numBoundaryLayers[2]);
		}
	private:

		// returns boundary layer indices of the current voxel in x and y
		// the indices are signed, indicating left or right side of the domain.
		// if no layer is touched: zero is returned
		V3i getBoundaryLayer( const Voxel& voxel );
		VoxelType& getVoxelType(int boundaryLayer_x, int boundaryLayer_y, int boundaryLayer_z);
		VoxelType& getVoxelType(const std::tuple<int, int, int>& boundaryLayer);

		PNSystem* sys;


		int m_numUnknowns; // number of cols and rows of A
		V3i m_numBoundaryLayers;
		V3i m_resolution;


		std::vector<Voxel> m_voxels;
		std::vector<VoxelType> m_voxelTypes;
		std::map<std::tuple<int,int,int>, int> m_layerToVoxelTypeIndex;
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

			V3i getVoxelCoord()
			{
				return voxel.coord;
			}

			V3d getOffset( int coeff )
			{
				return sys.getStencil().getOffset(coeff).cast<double>()*0.5;
			}

			const Domain& getDomain()
			{
				return sys.getDomain();
			}

			Fields& getFields()
			{
				return sys.getFields();
			}

			// return the offset of the given staggered grid in voxelspace
			const V3d& getGridOffset(int grid_index)
			{
				return g_grid_offsets[grid_index].cast<double>()*0.5;
			}

			const V3d& getGridOffset2(int grid_index)
			{
				return g_grid_offsets2[grid_index];
			}

			MatrixBuilderd::MatrixAccessHelper coeff_A(int coeff_i, const V3i &voxel_j, int coeff_j );
			MatrixBuilderd::MatrixAccessHelper coeff_b( int coeff_i );

		private:
			PNSystem& sys;
			Voxel& voxel;
			static V3i g_grid_offsets[8];
			static V3d g_grid_offsets2[8];
		};


		typedef std::function<void(Context&)> StencilFunction;
		typedef std::function<V3i(int)> GetCoefficientOffsetFunction;

		Stencil( int _order = -1,
				 int _numCoeffs = -1, // is different for the same order in 2d or 3d
				 int _width = 0,
				 StencilFunction fun = StencilFunction(),
				 GetCoefficientOffsetFunction getOffset_fun = GetCoefficientOffsetFunction() ):
			order(_order),
			numCoeffs(_numCoeffs),
			width(_width),
			apply(fun),
			getOffset(getOffset_fun)
		{

		}

		StencilFunction apply;
		GetCoefficientOffsetFunction getOffset;

		// the truncation order of the spherical harmonics expansion
		int order;
		// The number of SH coefficients per voxel. This number if different for the same order
		// depending on 2d or 3d
		int numCoeffs;

		// this tells us how many neighbouring voxels the stencil needs to see
		// this value is required for setting up boundary voxels
		// widht=0 means: no neighbours are touched
		int width;

		// this vector contains for every coefficient the staggered grid location index and therefore defines
		// where each coefficient is located
		//std::vector<int> coeffGridIndex;

		static Stencil noop()
		{
			return Stencil(0, 1, 2,
						   [](Context& sys){},
						   [](int coeff){return V3i(1,1,1);}); // noop places all coefficients at the cell center
		}
	};

	PNSystem(Stencil stencil, const Domain& domain, bool neumannBC);

	MatrixBuilderd::MatrixAccessHelper coeff_A(V3i voxel_i, int coefficient_i, V3i voxel_j, int coefficient_j, bool debug = false);
	MatrixBuilderd::MatrixAccessHelper coeff_b(V3i voxel_i, int coefficient_i, bool debug = false);



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
	Eigen::VectorXd solve_cg(const Eigen::VectorXd &x0);
	Eigen::VectorXd stripBoundary(const Eigen::VectorXd& x);


	const Domain& getDomain()const; // returns the spatial discretization used by the system
	int getNumCoefficients()const; // returns the number of SH coefficients per voxel
	int getNumVoxels()const; // returns the number of voxels
	int getOrder()const;


	// the following methods are used to allow python to access the build result
	// (used for storing to disk etc.)
	RealMatrix& get_A_real();
	RealVector& get_b_real();

	RealMatrix get_A_real_test();
	Eigen::VectorXd get_b_real_test();
	Eigen::VectorXd get_solve_convergence();

	void setDebugVoxel( const V3i& dv ); // this is used for debugging


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
	V3i debugVoxel;


	// computational domain, problem and stencil
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
	std::vector<double> m_solve_convergence;
};


















