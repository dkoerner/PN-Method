#pragma once

#include <memory>
#include <fields/Field.h>
//#include <util/data.h>
#include <math/frame.h>
#include <math/transform.h>
#include <math/bbox.h>
#include <math/rng.h>

#include <common/Domain.h>



struct PNVolume
{
	typedef std::shared_ptr<PNVolume> Ptr;

	PNVolume( const Domain& domain );

	void setExtinctionAlbedo( Field3d::Ptr extinction, Field3d::Ptr albedo );
	//void setAbsorptionScattering( Field3d::Ptr absorption, Field3d::Ptr scattering );

	void setEmission(int l, int m, Field3d::Ptr field);
	void setPhase( int l, int m, Field3d::Ptr field);


	V3d evalExtinction( const P3d& pWS, bool debug = false )const;
	V3d evalAbsorption( const P3d& pWS, bool debug = false )const;
	V3d evalScattering( const P3d& pWS, bool debug = false )const;
	V3d evalAlbedo( const P3d& pWS, bool debug = false )const;
	V3d evalEmission(int l, int m , const P3d& pWS)const;
	V3d evalPhase(int l, int m , const P3d& pWS)const;

	/*
	P3d localToWorld(const P3d& pLS)const;
	P3d worldToLocal(const P3d& pWS)const;

	void setLocalToWorld( const Transformd& localToWorld );
	void setBound( Box3d boundWS );
	const Box3d& getBound()const;
	*/
	Domain& getDomain();
private:
	Domain m_domain;
	// rgb dependent albedo and extinction values (evaluated in local space)
	Field3d::Ptr       m_field_extinction;
	Field3d::Ptr       m_field_albedo;
	std::vector<Field3d::Ptr> m_field_q; // field for each SH-coefficient of q
	std::vector<Field3d::Ptr> m_field_p; // field for each SH-coefficient of p (phase function)

	/*
	//PhaseFunction::Ptr m_phaseFunction;
	//PhaseFunctionField::Ptr phase_field;

	Transformd    m_localToWorld; // transform
	Transformd    m_worldToLocal; // transform
	Box3d         m_bboxLS; // bounding box in local space (0,0,0->1,1,1)
	Box3d         m_bboxWS; // bounding box in world space
	*/
};
