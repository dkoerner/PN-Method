#include <volume.h>


#include <volumes/hgridfield.h>

namespace volumes
{


	struct FieldVolume : public Volume
	{
		typedef std::shared_ptr<FieldVolume> Ptr;

		FieldVolume():
			Volume(),
			extinction_min(std::numeric_limits<double>::max()),
			extinction_max(-std::numeric_limits<double>::max())
		{
			bboxLS.reset();
			bboxLS.min = P3d(0.0f,0.0f,0.0f);
			bboxLS.max = P3d(1.0f,1.0f,1.0f);

			phase_field = phasefield_sggx();
		}

		// Volume overrides ----
		virtual Color3f evalExtinction( const P3d& pWS, bool debug = false )const override
		{

			P3d pLS = worldToLocal*pWS;
			V3d value = extinction_field->eval(pLS, debug);
			return Color3f( value.x(), value.y(), value.z() );
		}

		virtual Color3f evalAlbedo( const P3d& pWS, bool debug = false )const override
		{
			P3d pLS = worldToLocal*pWS;
			V3d value = albedo_field->eval(pLS);
			return Color3f( value.x(), value.y(), value.z() );
		}

		virtual double evalPhase( const P3d& pWS, const V3d& wi, const V3d& wo )const override
		{
			P3d pLS = worldToLocal*pLS;
			return phase_field->eval(pWS, wi, wo);
		}

		virtual double samplePhase( const P3d& pWS, const V3d& wi, V3d& wo, double& pdf, RNGd& rng )const override
		{
			P3d pLS = worldToLocal*pLS;
			return phase_field->sample(pLS, wi, wo, pdf, rng);
		}

		virtual Box3d getBound()const override
		{
			return bboxWS;
		}

		virtual bool intersectBound(const Ray3d& rayWS, double& mint, double& maxt, bool debug = false)const override
		{
			// intersect ray in local space to allow rotations etc.
			Ray3d rayLS = worldToLocal * rayWS;

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
			extinction_field = _extinction;
			_extinction->getValueRange(extinction_min, extinction_max);
		}

		void setAlbedo( Field3d::Ptr _albedo )
		{
			//_albedo->getValueRange(extinction_min, extinction_max);
			albedo_field = _albedo;
		}


		void setLocalToWorld( const Transformd& _localToWorld )
		{
			localToWorld = _localToWorld;
			worldToLocal = localToWorld.inverse();

			phase_field->setLocalToWorld(_localToWorld);

			// compute worldspace boundingbox
			bboxWS.reset();
			for( int i=0;i<8;++i )
				bboxWS.expandBy( localToWorld*bboxLS.getCorner(i) );
		}

		void setBound( Box3d boundWS )
		{
			setLocalToWorld( Transformd::from_aabb(boundWS) );
		}

		virtual Transformd getLocalToWorld()const override
		{
			return localToWorld;
		}


		// rgb dependent albedo and extinction values
		Field3d::Ptr extinction_field;
		Field3d::Ptr albedo_field;
		PhaseFunctionField::Ptr phase_field;

		V3d           extinction_min; // used for residual ratio tracking
		V3d           extinction_max; // used for distance sampling
		Transformd    localToWorld; // transform
		Transformd    worldToLocal; // transform
		BoundingBox3d bboxLS; // bounding box in local space
		BoundingBox3d bboxWS; // bounding box in world space
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
		//volume->setStepsize(stepSize);

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

		volume->setLocalToWorld(localToWorld);

		// use pure density
		volume->setExtinction( field::scalar_to_vector<double, float>( density ) );

		// fixed albedo
		volume->setAlbedo( field::constant(V3d(0.8)) );

		return volume;
	}

	Volume::Ptr homogeneous( const V3d& extinction, const V3d& albedo, Transformd localToWorld )
	{
		// volume ----
		FieldVolume::Ptr volume = std::make_shared<FieldVolume>();

		//volume->setLocalToWorld(localToWorld);
		volume->setBound( Box3d(P3d(-0.5), P3d(0.5)) );

		// use pure density
		volume->setExtinction(field::constant<V3d>(extinction));

		// fixed albedo
		volume->setAlbedo(field::constant(albedo));

		return volume;
	}

	Volume::Ptr scarf()
	{
		double scale = 100.0;
		///*
		FieldVolume::Ptr volume = std::make_shared<FieldVolume>();

		field::HGridField<double>::Ptr density = field::hgrid<double>( "c:/projects/epfl/data/scarf/data/volume_description.vol",
																	   "c:/projects/epfl/data/scarf/data/volume_$BX_$BY_$BZ-density.vol" );
		volume->setLocalToWorld(density->getLocalToWorld());
		volume->setExtinction( field::scalar_to_vector<double, double>( field::multiply_const<double, double>(density, scale) ) );
		//volume->setExtinction( field::scalar_to_vector<double, double>( density ) );
		volume->setAlbedo( field::constant(V3d(0.8)) );

		return volume;
		//*/

/*
		std::string basePath = ".";

		// volume ----
		FieldVolume::Ptr volume = std::make_shared<FieldVolume>();
		double stepSize;
		Transformd localToWorld;
		Fieldf::Ptr density = field::bgeo<float>(basePath + "/scarf.bgeo", &localToWorld, &stepSize);

		volume->setLocalToWorld(localToWorld);

		// use pure density
		volume->setExtinction( field::scalar_to_vector<double, float>( field::multiply_const<float>(density, 100.0f) ) );

		// fixed albedo
		volume->setAlbedo( field::constant(V3d(0.8)) );
*/
		return volume;
	}



} // namespace volumes
