#include <volume.h>









namespace volumes
{
	struct Isotropic : public PhaseFunction
	{
		typedef std::shared_ptr<Isotropic> Ptr;

		virtual Color3f sample( const V3d& wi, V3d& wo, double& pdf, RNGd& rng ) const override
		{
			wo = sampleSphere<double>(rng);
			pdf = INV_FOURPI;
			return Color3f(1.0); // eval/pdf = 1.0
		}

		virtual Color3f eval( const V3d& wi, const V3d& wo )const override
		{
			return Color3f(INV_FOURPI, INV_FOURPI, INV_FOURPI);
		}
	};

	struct FieldVolume : public Volume
	{
		typedef std::shared_ptr<FieldVolume> Ptr;

		FieldVolume():
			Volume(),
			extinction_min(std::numeric_limits<double>::max()),
			extinction_max(-std::numeric_limits<double>::max()),
			stepSize(0.1)
		{
			bboxLS.reset();
			bboxLS.min = P3d(0.0f,0.0f,0.0f);
			bboxLS.max = P3d(1.0f,1.0f,1.0f);

			albedo_field = albedo_transformed = std::make_shared<field::TransformField<V3d>>( Field<V3d>::Ptr(), Transformd() );
			extinction_field = extinction_transformed = std::make_shared<field::TransformField<V3d>>( Field<V3d>::Ptr(), Transformd() );
			phase_function = std::make_shared<Isotropic>();
		}

		// Volume overrides ----
		virtual Color3f evalExtinction( const P3d& pWS, bool debug = false )const override
		{
			V3d value = extinction_field->eval(pWS);
			return Color3f( value.x(), value.y(), value.z() );
		}

		virtual Color3f evalAlbedo( const P3d& pWS, bool debug = false )const override
		{
			V3d value = albedo_field->eval(pWS);
			return Color3f( value.x(), value.y(), value.z() );
		}

		virtual const PhaseFunction* getPhaseFunction()const
		{
			return phase_function.get();
		}

		virtual Box3d getBound()const override
		{
			return bboxWS;
		}

		virtual bool intersectBound(const Ray3d& rayWS, double& mint, double& maxt, bool debug = false)const override
		{
			// intersect ray in local space to allow rotations etc.
			Ray3d rayLS = localToWorld.inverse() * rayWS;

			if(debug)
			{
				std::cout << "FieldVolume::intersectBound: rayLS=" << rayLS.toString() << std::endl;
			}

			double tnearLS, tfarLS;
			if( bboxLS.rayIntersect(rayLS, tnearLS, tfarLS) )
			{
				mint = (localToWorld*rayLS(tnearLS) - rayWS.o).norm()*sign(tnearLS);
				maxt = (localToWorld*rayLS(tfarLS) - rayWS.o).norm()*sign(tfarLS);


				if(debug)
				{
					std::cout << "FieldVolume::intersectBound: mint=" << mint << std::endl;
					std::cout << "FieldVolume::intersectBound: maxt=" << maxt << std::endl;
					std::cout << "FieldVolume::intersectBound: bboxLS=" << bboxLS.toString() << std::endl;
				}


				mint = std::max( mint, rayWS.mint );
				maxt = std::min( maxt, rayWS.maxt );
				if(maxt-mint>0)
					return true;
			}
			return false;
		}

		virtual V3d getMaxExtinction()const override
		{
			return this->extinction_max;
		}


		// FieldVolume methods ----
		void setExtinction( Field3d::Ptr _extinction )
		{
			_extinction->getValueRange(extinction_min, extinction_max);

			extinction_transformed->m_input = _extinction;
		}

		void setAlbedo( Field3d::Ptr _albedo )
		{
			//_albedo->getValueRange(extinction_min, extinction_max);
			albedo_transformed->m_input = _albedo;
		}


		void setLocalToWorld( const Transformd& _localToWorld )
		{
			localToWorld = _localToWorld;

			albedo_transformed->m_xform = _localToWorld;
			extinction_transformed->m_xform = _localToWorld;

			// compute worldspace boundingbox
			bboxWS.reset();
			for( int i=0;i<8;++i )
				bboxWS.expandBy( localToWorld*bboxLS.getCorner(i) );
		}

		virtual Transformd getLocalToWorld()const override
		{
			return localToWorld;
		}

		void setStepsize( double _stepSize )
		{
			stepSize = _stepSize;
		}


		// rgb dependent albedo and extinction values
		Field3d::Ptr extinction_field;
		Field3d::Ptr albedo_field;
		PhaseFunction::Ptr phase_function;

		field::TransformField<V3d>::Ptr albedo_transformed;
		field::TransformField<V3d>::Ptr extinction_transformed;

		V3d           extinction_min; // used for residual ratio tracking
		V3d           extinction_max; // used for distance sampling
		Transformd    localToWorld; // transform
		BoundingBox3d bboxLS; // bounding box in local space
		BoundingBox3d bboxWS; // bounding box in world space
		double        stepSize; // raymarching stepsize
	};




	Volume::Ptr C60()
	{
		std::string basePath = "c:/projects/visus/data";

		// volume ----
		FieldVolume::Ptr volume = std::make_shared<FieldVolume>();
		double stepSize;
		Transformd localToWorld;
		Fieldf::Ptr density = field::bgeo<float>(basePath + "/datasets/C60Large.vol.bgeo", &localToWorld, &stepSize);
		volume->setLocalToWorld(localToWorld);
		volume->setStepsize(stepSize);

		// use pure density
		//volume.setExtinction( field::scalarToVector<double, float>( density ) );
		// use transfer function
		volume->setExtinction( field::transferFunction<double, float>( density,
																	  basePath + "/datasets/C60.tf.extinction", 100.0) );
		// fixed albedo
		volume->setAlbedo( field::constant(V3d(0.8)) );
		// derive albedo from transferfunction on extinction
		//volume->setAlbedo( field::transferFunction<double, float>( density, basePath + "/datasets/C60.tf.albedo") );

		return volume;
	}

	Volume::Ptr nebulae()
	{
		std::string basePath = "c:/projects/visus/data";

		// volume ----
		FieldVolume::Ptr volume = std::make_shared<FieldVolume>();
		double stepSize;
		Transformd localToWorld;
		Fieldf::Ptr density = field::bgeo<float>(basePath + "/datasets/nebulae200.bgeo", &localToWorld, &stepSize);
		//Fieldf::Ptr density = field::bgeo<float>("C:/projects/epfl/temp/grid.bgeo", &localToWorld, &stepSize);

		volume->setLocalToWorld(localToWorld);
		volume->setStepsize(stepSize);

		// use pure density
		volume->setExtinction( field::scalar_to_vector<double, float>( density ) );

		// fixed albedo
		volume->setAlbedo( field::constant(V3d(0.8)) );

		return volume;
	}

/*
	Scene rbf()
	{
		Scene ds;
		ds.id = "rbf";

		RBFSField::Ptr field_rbf = std::make_shared<RBFSField>();
		// add rbf
		P3d center(0.5, 0.5, 0.5);
		V3d scale(0.15, 0.15, 0.03);
		V3d rotation(0.989496, 0.52273, 1.1019); // angle axis
		double weight = 1.0;
		field_rbf->m_rbfs.push_back( RBFSField::RBF( center, scale, rotation, weight ) );

		ds.volume = std::make_shared<Volume>();
		ds.volume->setAlbedo( field::constant(V3d(1.0, 1.0, 1.0)) );
		ds.volume->setExtinction( field::scalar_to_vector<double, double>(field_rbf) );

		ds.bound = ds.volume->bboxWS;

		return ds;
	}


	Scene engine()
	{
		Scene ds;

		ds.id = "engine";
		std::string basePath = "c:/projects/visus/data";

		ds.volume = std::make_shared<Volume>();
		double stepSize;
		Transformd localToWorld;
		Fieldf::Ptr density = field::bgeo<float>(basePath + "/datasets/engine.bgeo", &localToWorld, &stepSize);
		ds.volume->setLocalToWorld(localToWorld);
		ds.volume->setStepsize(stepSize);

		// setup transfer function from Steffens tfa files ----
		field::TransferFunctionField3<double, float>::Ptr tf;
		//tf = field::transferFunction<double, float>( density, basePath + "/datasets/engine.tf.extinction", 1.0);

		//
		{
			double scale = 1.0;
			tf = std::make_shared<field::TransferFunctionField3<double, float>>(density);

			// load samples
			std::vector<double> samples;
			int numRows, numCols;
			readSamples2<double>( basePath + "/datasets/engine.tfa", samples, numRows, numCols );
			//std::cout << "transferFunction: numRows=" << numRows << " numCols=" << numCols << std::endl;

			float min, max;
			density->getValueRange(min, max);

			// add mapping
			tf->m_mapping.clear();
			for( int j=0;j<numRows;++j )
			{
				double t = double(j)/double(numRows-1);
				double in =(max-min)*t + min;
				//double r = samples[j*numCols+0];
				//double g = samples[j*numCols+1];
				//double b = samples[j*numCols+2];
				double a = samples[j*numCols+3]/255.0;
				TVector<double, 3>  temp(a*scale,
									a*scale,
									a*scale);
				//std::cout << "transferFunction: in=" << in << "temp=" << temp.toString() << std::endl;

				//std::cout << "tf: " << in << " " << temp.toString() << std::endl;

				tf->m_mapping.addSample(in, temp);
			}
		}

		// use pure density
		//ds.volume->setExtinction( field::scalar_to_vector<double, float>( density ) );
		// use transfer function
		ds.volume->setExtinction( tf );
		// fixed albedo
		//volume.setAlbedo( field::constant(V3d(0.8)) );
		// derive albedo from transferfunction on extinction
		ds.volume->setAlbedo( field::transferFunction<double, float>( density, basePath + "/datasets/chameleon.tf.albedo") );

		//ds.bound = ds.volume->bboxWS;
		ds.bound = Box3d( P3d(14.5, 14.5, -9), P3d(227.5, 227.5, 119) );

		return ds;
	}



	Scene chameleon()
	{
		Scene ds;
		ds.id = "chameleon";
		std::string basePath = "c:/projects/visus/data";
		ds.volume = std::make_shared<Volume>();
		double stepSize;
		Transformd localToWorld;
		Fieldf::Ptr density = field::bgeo<float>(basePath + "/datasets/chameleon_rotated.bgeo", &localToWorld, &stepSize);
		ds.volume->setLocalToWorld(localToWorld);
		ds.volume->setStepsize(stepSize);

		// setup transfer function from Steffens tfa files ----
		field::TransferFunctionField3<double, float>::Ptr tf;
		//tf = field::transferFunction<double, float>( density, basePath + "/datasets/chameleon.tf.extinction", 100.0);


		// use pure density
		//volume.setExtinction( field::scalarToVector<double, float>( density ) );
		// use transfer function
		ds.volume->setExtinction( tf );
		// fixed albedo
		//volume.setAlbedo( field::constant(V3d(0.8)) );
		// derive albedo from transferfunction on extinction
		ds.volume->setAlbedo( field::transferFunction<double, float>( density, basePath + "/datasets/chameleon.tf.albedo") );

		ds.bound = ds.volume->bboxWS;

		return ds;
	}
	*/
}
