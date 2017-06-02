
#include <volume.h>



namespace volumes
{

	struct Isotropic : public PhaseFunction
	{
		typedef std::shared_ptr<Isotropic> Ptr;

		virtual double sample( const V3d& wi, V3d& wo, double& pdf, RNGd& rng ) const override
		{
			wo = sampleSphere<double>(rng);
			pdf = INV_FOURPI;
			return 1.0; // eval/pdf = 1.0
		}

		virtual double eval( const V3d& wi, const V3d& wo )const override
		{
			return INV_FOURPI;
		}
	};


	PhaseFunction::Ptr phase_isotropic()
	{
		return std::make_shared<Isotropic>();
	}


} // namespace volumes
