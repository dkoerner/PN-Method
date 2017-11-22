#include <PNVolume.h>

#include <math/sph.h>









PNVolume::PNVolume(const Domain &domain):
	m_domain(domain)
	//m_worldToLocal(),
	//m_localToWorld()//,
	//m_bboxLS( P3d(0.0f,0.0f,0.0f), P3d(1.0f,1.0f,1.0f) ),
	//m_bboxWS( P3d(0.0f,0.0f,0.0f), P3d(1.0f,1.0f,1.0f) ),
	//m_phaseFunction( std::make_shared<IsotropicPhase>() )
{
	setExtinctionAlbedo( std::make_shared<ConstantField3d>(V3d(1.0, 1.0, 1.0)),
						 std::make_shared<ConstantField3d>(V3d(0.9, 0.9, 0.9)) );
}
void PNVolume::setExtinctionAlbedo( Field3d::Ptr extinction, Field3d::Ptr albedo )
{
	//std::pair<V3d, V3d> valueRange = extinction->getValueRange();
	//m_extinction_max = std::get<1>(valueRange);

	/*
	if( std::dynamic_pointer_cast<VoxelGridField3d>(extinction) )
	{
		VoxelGridField3d::Ptr vg = std::dynamic_pointer_cast<VoxelGridField3d>(extinction);
		V3d* data = vg->getData();
		int numElements = vg->getResolution()[0]*vg->getResolution()[1]*vg->getResolution()[2];
		m_extinction_max = V3d(0.0, 0.0, 0.0);
		for( int i=0;i<numElements;++i )
		{
			for( int j=0;j<3;++j )
				m_extinction_max[j] = std::max(m_extinction_max[j], data[i][j]);
		}
	}else
		m_extinction_max = extinction->getMaxValue();
	*/

	m_field_extinction = extinction;
	m_field_albedo = albedo;
}

void PNVolume::setEmission( int l, int m, Field3d::Ptr field)
{
	int index = sph::index(l, m);
	if( index>=m_field_q.size() )
		m_field_q.resize(index+1);
	m_field_q[index] = field;
}

void PNVolume::setPhase(int l, int m, Field3d::Ptr field)
{
	int index = sph::index(l, m);
	if( index>=m_field_p.size() )
		m_field_p.resize(index+1);
	m_field_p[sph::index(l, m)] = field;
}


Domain& PNVolume::getDomain()
{
	return m_domain;
}

V3d PNVolume::evalExtinction( const P3d& pWS, bool debug )const
{
	return m_field_extinction->eval(m_domain.worldToLocal(pWS));
}
V3d PNVolume::evalAbsorption( const P3d& pWS, bool debug )const
{
	V3d sigma_t = evalExtinction(pWS, debug);
	V3d albedo = evalAlbedo(pWS, debug);
	return V3d( sigma_t[0]*(1.0-albedo[0]),
				sigma_t[1]*(1.0-albedo[1]),
				sigma_t[2]*(1.0-albedo[2]) );
}
V3d PNVolume::evalScattering( const P3d& pWS, bool debug )const
{
	V3d sigma_t = evalExtinction(pWS, debug);
	V3d albedo = evalAlbedo(pWS, debug);
	return V3d( sigma_t[0]*albedo[0],
				sigma_t[1]*albedo[1],
				sigma_t[2]*albedo[2] );
}
V3d PNVolume::evalAlbedo( const P3d& pWS, bool debug )const
{
	return m_field_albedo->eval(m_domain.worldToLocal(pWS));
}

V3d PNVolume::evalEmission( int l, int m, const P3d& pWS )const
{
	int index = sph::index(l, m);
	if( index < m_field_q.size() && m_field_q[index] )
		return m_field_q[index]->eval(m_domain.worldToLocal(pWS));
	return V3d(0.0);
}

V3d PNVolume::evalPhase( int l, int m, const P3d& pWS )const
{
	int index = sph::index(l, m);
	if( index < m_field_p.size() && m_field_p[index] )
		return m_field_p[index]->eval(m_domain.worldToLocal(pWS));
	return V3d(0.0);
}
