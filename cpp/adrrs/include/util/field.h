#pragma once

#include <util/voxelgrid.h>
#include <util/data.h>
#include <math/plf.h>
#include <math/bbox.h>
#include <math/transform.h>

#include <houio/HouGeo.h>
#include <houio/HouGeoIO.h>


template<typename T>
struct Field
{
	typedef std::shared_ptr<Field> Ptr;

	// p is given in local space
	virtual T eval( const P3d& p, bool debug = false )const=0;
	virtual void getValueRange( T& min, T& max )const = 0;

	T getMax()const
	{
		T min, max;
		getValueRange(min, max);
		return max;
	}
};

typedef Field<float> Fieldf;
typedef Field<double> Fieldd;
typedef Field<TVector<double, 3>> Field3d;

namespace field
{

template<typename T>
struct ConstantField : public Field<T>
{
	typedef std::shared_ptr<ConstantField<T>> Ptr;

	ConstantField( T value ):
		Field<T>(),
		value(value)
	{

	}

	static Ptr create( T value )
	{
		return std::make_shared<ConstantField>( value );
	}

	virtual T eval( const P3d& p, bool debug = false )const override
	{
		return value;
	}

	virtual void getValueRange( T& min, T& max )const override
	{
		min = max = value;
	}

	T value;
};

typedef ConstantField<float> ConstantFieldf;
typedef ConstantField<double> ConstantFieldd;
typedef ConstantField<TVector<double,3>> ConstantField3d;

template<typename T>
typename Field<T>::Ptr constant( T value )
{
	return std::make_shared<ConstantField<T>>(value);
}


// ===========================================

template<typename T, typename R>
struct ScalarToVectorField : public Field<TVector<T, 3>>
{
	typedef std::shared_ptr<ScalarToVectorField<T, R>> Ptr;

	ScalarToVectorField():
		Field<TVector<T, 3>>()
	{

	}

//	static Ptr create( T value )
//	{
//		return std::make_shared<ConstantField>( value );
//	}

	virtual TVector<T, 3> eval( const P3d& p, bool debug = false )const override
	{
		return TVector<T, 3>(T(m_scalarField->eval(p)));
	}

	virtual void getValueRange( TVector<T, 3>& min, TVector<T, 3>& max )const override
	{
		R scalar_min, scalar_max;
		m_scalarField->getValueRange(scalar_min, scalar_max);
		min = TVector<T, 3>(T(scalar_min));
		max = TVector<T, 3>(T(scalar_max));
	}

	typename Field<R>::Ptr m_scalarField;
};

template<typename T, typename R>
typename ScalarToVectorField<T, R>::Ptr scalar_to_vector( typename Field<R>::Ptr scalarField )
{
	typename ScalarToVectorField<T, R>::Ptr vectorField = std::make_shared<ScalarToVectorField<T, R>>();
	vectorField->m_scalarField = scalarField;
	return vectorField;
}

// ===========================================


template<typename T, typename R>
struct VectorComponentField : public Field<T>
{
	typedef std::shared_ptr<VectorComponentField<T, R>> Ptr;

	VectorComponentField():
		Field<T>(),
		m_component(0)
	{

	}

	virtual T eval( const P3d& p, bool debug = false )const override
	{
		return T(m_vectorField->eval(p)[m_component]);
	}

	virtual void getValueRange( T& min, T& max )const override
	{
		TVector<R, 3> vector_min, vector_max;
		m_vectorField->getValueRange(vector_min, vector_max);
		min = T(vector_min[m_component]);
		max = T(vector_max[m_component]);
	}

	typename Field<TVector<R, 3>>::Ptr m_vectorField;
	int m_component;
};

template<typename T, typename R>
typename VectorComponentField<T, R>::Ptr vector_component( typename Field<TVector<R, 3>>::Ptr vectorField, int component )
{
	typename VectorComponentField<T, R>::Ptr scalarField = std::make_shared<VectorComponentField<T, R>>();
	scalarField->m_vectorField = vectorField;
	scalarField->m_component = component;
	return scalarField;
}

// ===========================================

template<typename T>
struct SphereField : public Field<T>
{
	typedef std::shared_ptr<SphereField<T>> Ptr;

	SphereField( T value ):
		Field<T>(),
		value(value)
	{

	}

	virtual T eval( const P3d& p, bool debug = false )const override
	{
		if( (P3d(0.5, 0.5, 0.5)-p).norm() < 0.5 )
			return value;
		return 0.0;
	}

	virtual void getValueRange( T& min, T& max )const override
	{
		min = max = value;
	}

	T value;
};

typedef SphereField<double> SphereFieldd;

template<typename T, typename TGRID=T>
struct VoxelGridField : public Field<T>
{
	typedef std::shared_ptr<VoxelGridField> Ptr;

	VoxelGridField( V3i res, TGRID *data = 0 ):
		Field<T>()
	{
		grid.resize( res.x(), res.y(), res.z() );
		if(data)
			memcpy( &grid.m_data[0], data, res.x()*res.y()*res.z()*sizeof(TGRID) );
	}

	//p is given in localspace
	virtual T eval( const P3d& p, bool debug = false )const override
	{
		//return T(grid.evaluate(grid.localToVoxel(p)));

		// grid.evaluate will cast its result to T
		// which means it will return int after interpolation
		// to avoid this we do the interpolation here

		typedef double real_t;
		typedef Vector3d Vector;

		Vector vs = grid.localToVoxel(p);

		// take sample location within voxel into account
		vs -= grid.m_sampleLocation.template cast<real_t>();

		real_t tx = vs.x() - floor(vs.x());
		real_t ty = vs.y() - floor(vs.y());
		real_t tz = vs.z() - floor(vs.z());

		// lower left corner
		V3i c1;
		c1[0] = (int)floor(vs.x());
		c1[1] = (int)floor(vs.y());
		c1[2] = (int)floor(vs.z());

		// upper right corner
		V3i c2 = c1+V3i(1);
		V3i res = grid.getResolution();

		// clamp the indexing coordinates
		c1[0] = std::max(0, std::min(c1.x(), res.x()-1));
		c2[0] = std::max(0, std::min(c2.x(), res.x()-1));
		c1[1] = std::max(0, std::min(c1.y(), res.y()-1));
		c2[1] = std::max(0, std::min(c2.y(), res.y()-1));
		c1[2] = std::max(0, std::min(c1.z(), res.z()-1));
		c2[2] = std::max(0, std::min(c2.z(), res.z()-1));

		//lerp...
		return lerp<T>( lerp<T>( lerp<T>( grid.sample( c1.x(), c1.y(), c1.z() ),
										  grid.sample( c2.x(), c1.y(), c1.z() ), (real_t)tx ),
								 lerp<T>( grid.sample( c1.x(), c2.y(), c1.z() ),
										  grid.sample( c2.x(), c2.y(), c1.z() ), (real_t)tx ), (real_t)ty ),
						lerp<T>( lerp<T>( grid.sample( c1.x(), c1.y(), c2.z() ),
										  grid.sample( c2.x(), c1.y(), c2.z() ), (real_t)tx ),
								 lerp<T>( grid.sample( c1.x(), c2.y(), c2.z() ),
										  grid.sample( c2.x(), c2.y(), c2.z() ), (real_t)tx ), (real_t)ty ), (real_t)tz );
	}

	virtual void getValueRange( T& min, T& max )const override
	{
		auto result = std::minmax_element(grid.m_data.begin(), grid.m_data.end());
		min = T(*result.first);
		max = T(*result.second);
	}


	VoxelGrid<TGRID> grid;
};
typedef VoxelGridField<float> VoxelGridFieldf;

struct VoxelGridField4f : public Field<TVector<double, 4>>
{
	typedef std::shared_ptr<VoxelGridField4f> Ptr;

	VoxelGridField4f( V3i res, float *data = 0 ):
		Field()
	{
		grid.resize( res.x(), res.y(), res.z() );
		if(data)
			memcpy( &grid.m_data[0], data, res.x()*res.y()*res.z()*sizeof(float)*4 );
		else
			grid.fill(Vector4f(0.0,0.0,0.0,0.0));

	}

	VoxelGridField4f( const VoxelGrid4f& other ):
		Field()
	{
		V3i res = other.getResolution();
		grid.resize( res.x(), res.y(), res.z() );
		memcpy( &grid.m_data[0], other.getRawPointer(), res.x()*res.y()*res.z()*sizeof(float)*4 );
	}

	//p is given in localspace
	virtual TVector<double, 4> eval( const P3d& p, bool debug = false )const override
	{
		return grid.evaluate(grid.localToVoxel(p)).cast<double>();
	}

	virtual void getValueRange( TVector<double, 4>& min, TVector<double, 4>& max )const override
	{
//		auto result = std::minmax_element(grid.m_data.begin(), grid.m_data.end());
//		min = *result.first;
//		max = *result.second;
	}


	VoxelGrid4f grid;
};


template<typename T, typename R>
struct TransferFunctionField3 : public Field<TVector<T, 3>>
{
	typedef std::shared_ptr<TransferFunctionField3<T, R>> Ptr;

	TransferFunctionField3( typename Field<R>::Ptr scalarField = typename Field<R>::Ptr() ):
		Field<TVector<T, 3>>(),
		m_scalarfield(scalarField)
	{
	}

	virtual TVector<T, 3> eval( const P3d& p, bool debug = false )const override
	{
		double scalar = m_scalarfield->eval(p);
		if(debug)
			std::cout << "TransferFunctionField3 scalar=" << scalar << std::endl;
		return m_mapping.evaluate( scalar );
	}

	virtual void getValueRange( TVector<T, 3>& min, TVector<T, 3>& max )const override
	{
		// returns min and max for each channel
		min = TVector<T, 3>(std::numeric_limits<T>::max());
		max = TVector<T, 3>(-std::numeric_limits<T>::max());
		for( const auto& v:m_mapping.m_values )
		{
			for(int i=0;i<3;++i)
			{
				min[i] = std::min(min[i], v[i]);
				max[i] = std::max(max[i], v[i]);
			}
		}
	}


	typename Field<R>::Ptr m_scalarfield; // 1d values which will be mapped to 3d values
	PLF<R, TVector<T, 3>> m_mapping; // transfer function maps T to vector<T>
};

typedef TransferFunctionField3<double, float> TransferFunctionField3df;

template<typename T, typename R>
typename TransferFunctionField3<T, R>::Ptr transferFunction( typename Field<R>::Ptr scalarField, const std::string& filename, double scale = 1.0 )
{
	typename TransferFunctionField3<T, R>::Ptr result = std::make_shared<TransferFunctionField3<T, R>>(scalarField);

	// load samples
	std::vector<double> samples;
	int numRows, numCols;
	readSamples2<T>( filename, samples, numRows, numCols );
	//std::cout << "transferFunction: numRows=" << numRows << " numCols=" << numCols << std::endl;

	// add mapping
	result->m_mapping.clear();
	for( int j=0;j<numRows;++j )
	{
		double in = samples[j*4];
		TVector<T, 3> temp(samples[j*4+1]*scale,
						   samples[j*4+2]*scale,
						   samples[j*4+3]*scale);
		//std::cout << "transferFunction: in=" << in << "temp=" << temp.toString() << std::endl;
		result->m_mapping.addSample(in, temp);
	}

	return result;
}




/*
HeterogeneousMedium::Ptr HeterogeneousMedium::load(const std::string &path)
{
	std::cout << "HeterogeneousMedium::load: loading " << path;
	houio::ScalarField::Ptr volume = houio::HouGeoIO::importVolume( path );

	// extract transform
	Eigen::Affine3d transform = Eigen::Affine3d::Identity();
	Eigen::Matrix4d m;
	m << volume->m_localToWorld.ma[0], volume->m_localToWorld.ma[1], volume->m_localToWorld.ma[2], volume->m_localToWorld.ma[3],
		 volume->m_localToWorld.ma[4], volume->m_localToWorld.ma[5], volume->m_localToWorld.ma[6], volume->m_localToWorld.ma[7],
		 volume->m_localToWorld.ma[8], volume->m_localToWorld.ma[9], volume->m_localToWorld.ma[10], volume->m_localToWorld.ma[11],
		 volume->m_localToWorld.ma[12], volume->m_localToWorld.ma[13], volume->m_localToWorld.ma[14], volume->m_localToWorld.ma[15];
	transform = Eigen::Affine3d(m.transpose());

	// extract voxeldata
	VoxelGridField::Ptr sigmaT = std::make_shared<VoxelGridField>( V3i(volume->getResolution().x, volume->getResolution().y, volume->getResolution().z),
														volume->getRawPointer());

	std::cout << " resolution: " << volume->getResolution().x << "x" << volume->getResolution().y << "x" << volume->getResolution().z << std::endl;

	PhaseFunction::Ptr phase = std::make_shared<Isotropic>();
	Fieldd::Ptr albedo = std::make_shared<ConstantFieldd>(0.9);
	Transformd localToWorld(transform);
	return std::make_shared<HeterogeneousMedium>(  phase, sigmaT, albedo, localToWorld );
	return HeterogeneousMedium::Ptr();
}
*/


	//VoxelGridFieldf::Ptr bgeo(const std::string &path, Transformd *localToWorld=0, double *voxelSize=0);


	template<typename T>
	struct SlabField : public Field<T>
	{
		typedef std::shared_ptr<SlabField<T>> Ptr;

		SlabField( const V3d& plane_normal, double plane_distance_origin, double thickness ):
			Field<T>(),
			plane_normal(plane_normal),
			plane_distance_origin(plane_distance_origin),
			thickness(thickness)
		{

		}

		static Ptr create( const V3d& plane_normal, double plane_distance_origin, double thickness )
		{
			return std::make_shared<SlabField>( plane_normal, plane_distance_origin, thickness );
		}

		virtual T eval( const P3d& p, bool debug = false )const override
		{
			double plane_distance = p.dot(plane_normal) + plane_distance_origin;
			if( (plane_distance > 0.0) && (plane_distance < thickness) )
				return T(1.0);
			return T(0.0);
		}

		virtual void getValueRange( T& min, T& max )const override
		{
			min = T(0.0);
			max = T(1.0);
		}

		V3d plane_normal;
		double plane_distance_origin;
		double thickness;
	};

	template<typename T>
	typename SlabField<T>::Ptr slab( const V3d& plane_normal, double plane_distance_origin, double thickness )
	{
		return SlabField<T>::create(plane_normal, plane_distance_origin, thickness);
	}

	// ==================================================
	template<typename T, typename R>
	struct MultiplyConstantField : public Field<T>
	{
		typedef std::shared_ptr<MultiplyConstantField<T, R>> Ptr;

		MultiplyConstantField( typename Field<T>::Ptr input, R constant ):
			Field<T>(),
			m_input(input),
			m_constant(constant)
		{
		}

		static Ptr create( typename Field<T>::Ptr input, R constant )
		{
			return std::make_shared<MultiplyConstantField>(input, constant);
		}

		virtual T eval( const P3d& p, bool debug = false )const override
		{
			return m_input->eval(p, debug)*m_constant;
		}

		virtual void getValueRange( T& min, T& max )const override
		{
			m_input->getValueRange(min, max);
			min *= m_constant;
			max *= m_constant;
		}


		typename Field<T>::Ptr m_input;
		R m_constant;
	};

	template<typename T, typename R>
	typename MultiplyConstantField<T, R>::Ptr multiply_const( typename Field<T>::Ptr input, R constant )
	{
		return MultiplyConstantField<T, R>::create(input, constant);
	}


	// ========================================================
	template<typename T>
	struct MultiplyField : public Field<T>
	{
		typedef std::shared_ptr<MultiplyField<T>> Ptr;

		MultiplyField( typename Field<T>::Ptr a, typename Field<T>::Ptr b ):
			Field<T>(),
			m_a(a),
			m_b(b)
		{
		}

		static Ptr create( typename Field<T>::Ptr a, typename Field<T>::Ptr b )
		{
			return std::make_shared<MultiplyField>(a, b);
		}

		virtual T eval( const P3d& p, bool debug = false )const override
		{
			return m_a->eval(p, debug)*m_b->eval(p, debug);
		}

		virtual void getValueRange( T& min, T& max )const override
		{
			m_a->getValueRange(min, max);
			T b_min, b_max;
			m_b->getValueRange(b_min, b_max);
			min *= b_min;
			max *= b_max;
		}


		typename Field<T>::Ptr m_a;
		typename Field<T>::Ptr m_b;
	};

	template<typename T>
	typename MultiplyField<T>::Ptr multiply( typename Field<T>::Ptr a, typename Field<T>::Ptr b )
	{
		return MultiplyField<T>::create(a, b);
	}


	// bgeo import ======================================
	template<typename T>
	typename VoxelGridField<T, float>::Ptr bgeo(const std::string &path, Transformd* localToWorld=0, double* voxelSize=0 )
	{
		std::cout << "bgeo: loading " << path << std::endl;
		houio::ScalarField::Ptr grid = houio::HouGeoIO::importVolume( path );

		if(!grid)
			throw std::runtime_error("bgeo: unable to load " + path);

		// extract transform
		if(localToWorld)
		{
			Eigen::Matrix4d m;
			m << grid->m_localToWorld.ma[0], grid->m_localToWorld.ma[1], grid->m_localToWorld.ma[2], grid->m_localToWorld.ma[3],
				 grid->m_localToWorld.ma[4], grid->m_localToWorld.ma[5], grid->m_localToWorld.ma[6], grid->m_localToWorld.ma[7],
				 grid->m_localToWorld.ma[8], grid->m_localToWorld.ma[9], grid->m_localToWorld.ma[10], grid->m_localToWorld.ma[11],
				 grid->m_localToWorld.ma[12], grid->m_localToWorld.ma[13], grid->m_localToWorld.ma[14], grid->m_localToWorld.ma[15];
			*localToWorld = Transformd(Eigen::Affine3d(m.transpose()));
		}

		if(voxelSize)
		{
			*voxelSize = grid->getVoxelSize().getLength();
		}


		// extract voxeldata
		V3i res(grid->getResolution().x, grid->getResolution().y, grid->getResolution().z);
		float* data = grid->getRawPointer();
		typename VoxelGridField<T, float>::Ptr grid_new = std::make_shared<VoxelGridField<T, float>>(res, data);

		return grid_new;
	}

	template<typename T>
	typename Field<T>::Ptr read(const std::string &path )
	{
		std::cout << "bgeo: loading " << path << std::endl;
		houio::ScalarField::Ptr grid = houio::HouGeoIO::importVolume( path );

		if(!grid)
			throw std::runtime_error("bgeo: unable to load " + path);

		// extract transform
		Eigen::Matrix4d m;
		m << grid->m_localToWorld.ma[0], grid->m_localToWorld.ma[1], grid->m_localToWorld.ma[2], grid->m_localToWorld.ma[3],
			 grid->m_localToWorld.ma[4], grid->m_localToWorld.ma[5], grid->m_localToWorld.ma[6], grid->m_localToWorld.ma[7],
			 grid->m_localToWorld.ma[8], grid->m_localToWorld.ma[9], grid->m_localToWorld.ma[10], grid->m_localToWorld.ma[11],
			 grid->m_localToWorld.ma[12], grid->m_localToWorld.ma[13], grid->m_localToWorld.ma[14], grid->m_localToWorld.ma[15];
		Transformd localToWorld = Transformd(Eigen::Affine3d(m.transpose()));

		// extract voxeldata
		V3i res(grid->getResolution().x, grid->getResolution().y, grid->getResolution().z);
		float* data = grid->getRawPointer();
		typename VoxelGridField<T, float>::Ptr grid_new = std::make_shared<VoxelGridField<T, float>>(res, data);

		return xform<T>(grid_new, localToWorld);
	}

	// export
	//void write(const std::string& filename, Field<double>::Ptr field, const Box3d &bound, const V3i& res);
	//void write(const std::string& filename, const Transformd& localToWorld, const V3i& res, const Box3d &bound_ls=Box3d(P3d(0.0, 0.0, 0.0), P3d(1.0, 1.0, 1.0)) );
	void write(const std::string& filename, Field<double>::Ptr field, const V3i& res, const Transformd& localToWorld = Transformd(), const Box3d &bound_ls=Box3d(P3d(0.0, 0.0, 0.0), P3d(1.0, 1.0, 1.0)));


	// =============================================

	template<typename T>
	struct TransformField : public Field<T>
	{
		typedef std::shared_ptr<TransformField<T>> Ptr;

		TransformField( typename Field<T>::Ptr input, Transformd xform ):
			Field<T>(),
			m_input(input),
			m_xform(xform)
		{
		}

		static Ptr create( typename Field<T>::Ptr input, Transformd xform )
		{
			return std::make_shared<TransformField>(input, xform);
		}

		virtual T eval( const P3d& p, bool debug = false )const override
		{
			P3d p2 = m_xform.inverse()*p;
			if(debug)
				std::cout << "TransformField::eval p=" << p2.toString() << std::endl;
			return m_input->eval(p2, debug);
		}

		virtual void getValueRange( T& min, T& max )const override
		{
			m_input->getValueRange(min, max);
		}


		typename Field<T>::Ptr m_input;
		Transformd m_xform;
	};

	template<typename T>
	typename Field<T>::Ptr xform( typename Field<T>::Ptr input, Transformd localToWorld )
	{
		return std::make_shared<TransformField<T>>( input, localToWorld );
	}

	// transforms the input field such that the local space will be mapped to the given aa-bound
	template<typename T>
	typename Field<T>::Ptr xform( typename Field<T>::Ptr input, const Box3d& bound_dest, const Box3d& bound_src = Box3d(P3d(0.0), P3d(1.0)) )
	{
		V3d extent_src = bound_src.getExtents();
		V3d extent_dest = bound_dest.getExtents();

		Eigen::Transform<T,3,Eigen::Affine> test =
				Eigen::Translation<double, 3>(bound_dest.min-bound_src.min)*
				Eigen::DiagonalMatrix<double, 3>(extent_dest.cwiseQuotient(extent_src));
		Transformd xform( test );
		return std::make_shared<TransformField<T>>( input, xform );
	}


	// ================================================================================



	template<typename T>
	struct GaussianField : public Field<T>
	{
		typedef std::shared_ptr<GaussianField> Ptr;

		struct Gaussian
		{
			Gaussian():
				weight(1.0),
				center(0.0),
				stddev()
			{
			}
			Gaussian( const TPoint<T, 3>& center, T stddev, T weight ):
				weight(weight),
				center(center),
				stddev(stddev)
			{
			}

			T eval(const TPoint<T, 3> &p, bool debug=false) const
			{
				// transform into local space
				TVector<T, 3> d = p-center;
				T g = 0.5/(stddev*stddev);

				T t = d.dot(d); // squared distance
				return weight*std::exp(-g*t);
			}

			std::string toString()const
			{
				std::stringstream ss;
				//ss << "weight=" << weight << " center=" << center.toString() << " scale=" << scale.toString() << std::endl;
				return ss.str();
			}

			T weight;
			TPoint<T, 3> center;
			T stddev;
		};

		GaussianField():
			Field<T>(),
			m_gaussian()
		{
		}

		GaussianField( const Gaussian& gaussian ):
			Field<T>(),
			m_gaussian(gaussian)
		{
		}



		virtual T eval(const TPoint<T, 3> &p, bool debug = 0) const override
		{
			return m_gaussian.eval(p, debug);

		}



		virtual void getValueRange( T& min, T& max )const override
		{
			min = 0.0;
			max = 0.0;
		}

		Gaussian m_gaussian;
	};

	template<typename T>
	typename GaussianField<T>::Ptr gaussian( const TPoint<T, 3>& center, T stddev, T weight )
	{
		return std::make_shared<GaussianField<T>>( GaussianField<T>::Gaussian( center, stddev, weight ) );
	}

	// ======================================================================================================

	template<typename T>
	struct AnisotropicGaussianField : public Field<T>
	{
		typedef std::shared_ptr<AnisotropicGaussianField> Ptr;

		struct AnisotropicGaussian
		{
			AnisotropicGaussian():
				weight(1.0),
				center(0.0),
				cov(Eigen::Matrix<double, 3, 3>::Identity()),
				cov_inverse(Eigen::Matrix<double, 3, 3>::Identity())
			{
			}
			AnisotropicGaussian( const TPoint<T, 3>& center, Eigen::Matrix<T, 3, 3> cov, T weight ):
				weight(weight),
				center(center),
				cov(cov),
				cov_inverse(cov.inverse())
			{

			}

			T eval(const TPoint<T, 3> &p, bool debug=false) const
			{
				// transform into local space
				V3d d = p-center;

				// now evaluate rbf ---
				double t = d.dot(cov_inverse*d);
				return weight*std::exp(-0.5*t);
			}

			std::string toString()const
			{
				std::stringstream ss;
				//ss << "weight=" << weight << " center=" << center.toString() << " scale=" << scale.toString() << std::endl;
				return ss.str();
			}

			T weight;
			TPoint<T, 3> center;
			Eigen::Matrix<T, 3, 3> cov; // covariance matrix
			Eigen::Matrix<T, 3, 3> cov_inverse;
		};

		AnisotropicGaussianField():
			Field<T>(),
			m_anisotropic_gaussian()
		{
		}

		AnisotropicGaussianField( const AnisotropicGaussian& gaussian ):
			Field<T>(),
			m_anisotropic_gaussian(gaussian)
		{
		}



		virtual T eval(const TPoint<T, 3> &p, bool debug = 0) const override
		{
			return m_anisotropic_gaussian.eval(p, debug);

		}



		virtual void getValueRange( T& min, T& max )const override
		{
			min = 0.0;
			max = 0.0;
		}

		AnisotropicGaussian m_anisotropic_gaussian;
	};

	template<typename T>
	typename AnisotropicGaussianField<T>::Ptr anisotropic_gaussian( const TPoint<T, 3>& center, const Eigen::Matrix<double, 3, 3>& stddev, T weight )
	{
		return std::make_shared<AnisotropicGaussianField<T>>( AnisotropicGaussianField<T>::AnisotropicGaussian( center, stddev, weight ) );
	}

} // namespace field



struct RBFSField : public Field<double>
{
	typedef std::shared_ptr<RBFSField> Ptr;
	struct RBF
	{
		RBF():
			weight(1.0),
			scale(1.0),
			rotation(0.0, 0.0, 0.0),
			center(P3d(0.0, 0.0, 0.0))
		{
		}

		RBF( const P3d& center, const V3d& scale, V3d rotation, double weight=1.0 ):
			weight(weight),
			rotation(rotation),
			scale(scale),
			center(center)
		{
			/*
			Eigen::Matrix<double, 3, 3> r;
			r =  Eigen::AngleAxis<double>(rotation.x(), TVector<double, 3>::UnitX())
				*Eigen::AngleAxis<double>(rotation.y(), TVector<double, 3>::UnitY())
				*Eigen::AngleAxis<double>(rotation.z(), TVector<double, 3>::UnitZ());

			Eigen::Matrix<double, 3, 3> s = V3d(1.0/(scale.x()*scale.x()), 1.0/(scale.y()*scale.y()), 1.0/(scale.z()*scale.z())).asDiagonal();

			toLocal = r*s*r.transpose();
			*/
		}

		template<typename T>
		static void sincosx(T theta, T *_sin, T *_cos)
		{
			*_sin = std::sin(theta);
			*_cos = std::cos(theta);
		}

		template<typename T>
		static TVector<T, 3> sphericalDirectionx(T theta, T phi)
		{
			T sinTheta, cosTheta, sinPhi, cosPhi;

			sincosx(theta, &sinTheta, &cosTheta);
			sincosx(phi, &sinPhi, &cosPhi);

			return TVector<T, 3>(
			   sinTheta * cosPhi,
			   sinTheta * sinPhi,
			   cosTheta
			);
		}

		template<typename T>
		static T evalx( const T* const p, const T* const center, const T* const scale, const T* const rotation, const T& weight )
		{
			// build ellipsoidal matrix ---

			// rotation
			Eigen::Matrix<T, 3, 3> r;

			// angle axis parameterization
			r =  Eigen::AngleAxis<T>(rotation[2], sphericalDirectionx<T>(rotation[0], rotation[1]));

			// euler angles
			//r =  Eigen::AngleAxis<T>(rotation[0], TVector<T, 3>::UnitX())
			//	*Eigen::AngleAxis<T>(rotation[1], TVector<T, 3>::UnitY())
			//	*Eigen::AngleAxis<T>(rotation[2], TVector<T, 3>::UnitZ());

			Eigen::Matrix<T, 3, 3> s = TVector<T, 3>(1.0/(scale[0]*scale[0]), 1.0/(scale[1]*scale[1]), 1.0/(scale[2]*scale[2])).asDiagonal();
			Eigen::Matrix<T, 3, 3> toLocal = r*s*r.transpose();

			// now evaluate rbf ---
			TVector<T, 3> d( p[0]-center[0], p[1]-center[1], p[2]-center[2] );
			//return weight*ceres::exp(-T(0.5)*d.dot(toLocal*d));
			return weight*std::exp(-T(0.5)*d.dot(toLocal*d));
		}


		double eval( const P3d& p )const
		{
			return evalx<double>( &p[0], &center[0], &scale[0], &rotation[0], weight );
		}

		std::string toString()const
		{
			std::stringstream ss;
			ss << "weight=" << weight << " center=" << center.toString() << " scale=" << scale.toString() << " rotation=" << rotation.toString() << std::endl;
			return ss.str();
		}

		void test()
		{
			// angle axis parameterization
			rotationMatrix =  Eigen::AngleAxis<double>(rotation[2], sphericalDirectionx<double>(rotation[0], rotation[1]));
		}

		// parameters
		double weight;
		V3d rotation; // rotation in angle axis format (theta, phi, angle in rad)
		V3d scale;
		P3d center;

		//Eigen::Matrix3d toLocal;
		//
		Eigen::Matrix<double, 3, 3> rotationMatrix;
	};


	RBFSField():
		Field<double>()
	{
	}

	virtual double eval(const P3d &p, bool debug = 0) const override
	{
		double result = 0.0;
		for( auto& rbf:m_rbfs )
			result+=rbf.eval(p);
		return result;
	}

	virtual void getValueRange( double& min, double& max )const override
	{
		min = 0.0;
		max = 0.0;
		for( auto& rbf:m_rbfs )
			max = std::max(max, rbf.eval(rbf.center));
	}



	std::vector<RBF> m_rbfs;
};


struct RBFField : public Field<double>
{
	typedef std::shared_ptr<RBFField> Ptr;

	// rbf is based on 3d anisotropic gauss
	struct RBF
	{
		RBF():
			enabled(false),
			weight(1.0),
			center(0.0),
			cov_inverse(Eigen::Matrix<double, 3, 3>::Identity())
		{
		}

		double eval(const P3d &p, bool debug=false) const
		{
			// transform into local space
			V3d d = p-center;

			// now evaluate rbf ---
			double t = d.dot(cov_inverse*d);
			return weight*std::exp(-0.5*t);
		}

		std::string toString()const
		{
			std::stringstream ss;
			//ss << "weight=" << weight << " center=" << center.toString() << " scale=" << scale.toString() << std::endl;
			return ss.str();
		}




		bool enabled;
		double weight;
		Point3d center;
		Eigen::Matrix<double, 3, 3> cov_inverse; // inverse of the covariance matrix
	};

	RBFField():
		Field<double>(),
		m_rbf()
	{
	}

	RBFField( const RBF& rbf ):
		Field<double>(),
		m_rbf(rbf)
	{
	}



	virtual double eval(const P3d &p, bool debug = 0) const override
	{
		return m_rbf.eval(p, debug);

	}



	virtual void getValueRange( double& min, double& max )const override
	{
		min = 0.0;
		max = 0.0;
	}

	RBF m_rbf;
};
