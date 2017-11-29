#pragma once



#include <integrator.h>

#include <scene.h>
#include <math/frame.h>

#include <volume.h>


#include <PNSolution.h>


// intersects the volume bound and starts volume path tracing within the volume
// will stop once the path exits the volume
struct DirectPN : public Integrator
{
	typedef std::shared_ptr<DirectPN> Ptr;

	DirectPN( PNSolution* pns,
	          bool doSingleScattering = true,
	          bool doMultipleScattering = true );

	V3d trace( TraceInfo& ti, RNGd& rng )const;

	virtual V3d Li( const Scene* scene, RadianceQuery& rq, RNGd& rng )const override;
	virtual V3d sample_transmittance( const Scene* scene, const Ray3d& ray, double maxt, RNGd& rng )const override
	{
		V3d sigma_t;
		double distance = delta_tracking( scene, ray, maxt, 0, rng, sigma_t );

		if( distance < maxt )
			// since we return 0, we dont need to produce the pdf in the denominator
			return V3d(0.0, 0.0, 0.0);

		//NB: transmittance sampling pdf cancels out with the transmittance term

		return V3d(1.0, 1.0, 1.0);
	}

	virtual std::string toString()const override
	{
		std::ostringstream ss;
		ss << "DirectPN singleScattering=" << m_doSingleScattering << " multipleScattering=" << m_doMultipleScattering << std::endl;
		return ss.str();
	}

private:
	//PNSolution::Ptr m_pns;
	PNSolution* m_pns;

	bool m_doSingleScattering;
	bool m_doMultipleScattering;
};
