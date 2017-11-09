#include <volume.h>






Volume::Volume():
	m_worldToLocal(),
	m_localToWorld(),
	m_bboxLS( P3d(0.0f,0.0f,0.0f), P3d(1.0f,1.0f,1.0f) ),
	m_bboxWS( P3d(0.0f,0.0f,0.0f), P3d(1.0f,1.0f,1.0f) ),
	m_phaseFunction( std::make_shared<IsotropicPhase>() )
{
	setExtinctionAlbedo( std::make_shared<ConstantField3d>(V3d(1.0, 1.0, 1.0)),
						 std::make_shared<ConstantField3d>(V3d(0.9, 0.9, 0.9)) );
}
void Volume::setExtinctionAlbedo( Field3d::Ptr extinction, Field3d::Ptr albedo )
{
	std::pair<V3d, V3d> valueRange = extinction->getValueRange();
	m_extinction_max = std::get<1>(valueRange);

	m_field_extinction = extinction;
	m_field_albedo = albedo;
}
void Volume::setLocalToWorld( const Transformd& localToWorld )
{
	m_localToWorld = localToWorld;
	m_worldToLocal = m_localToWorld.inverse();

	//phase_field->setLocalToWorld(_localToWorld);

	// compute worldspace boundingbox
	m_bboxWS.reset();
	for( int i=0;i<8;++i )
		m_bboxWS.expandBy( m_localToWorld*m_bboxLS.getCorner(i) );
}
P3d Volume::localToWorld(const P3d& pLS)const
{
	return m_localToWorld*pLS;
}

P3d Volume::worldToLocal(const P3d& pWS)const
{
	return m_worldToLocal*pWS;
}

void Volume::setBound( Box3d boundWS )
{
	setLocalToWorld( Transformd::from_aabb(boundWS) );
}
const Box3d& Volume::getBound()const
{
	return this->m_bboxWS;
}
V3d Volume::evalExtinction( const P3d& pWS, bool debug )const
{
	return m_field_extinction->eval(worldToLocal(pWS));
}
V3d Volume::evalAlbedo( const P3d& pWS, bool debug )const
{
	return m_field_albedo->eval(worldToLocal(pWS));
}
double Volume::evalPhase( const P3d& pWS, const V3d& wi, const V3d& wo )const
{
	return m_phaseFunction->eval(wi, wo);
}
double Volume::samplePhase( const P3d& pWS, const V3d& wi, V3d& wo, double& pdf, RNGd& rng )const
{
	return m_phaseFunction->sample(wi, wo, pdf, rng);
}
bool Volume::intersectBound(const Ray3d& rayWS, double& mint, double& maxt, bool debug)const
{
	// intersect ray in local space to allow rotations etc.
	Ray3d rayLS = m_worldToLocal * rayWS;

	if(debug)
	{
		std::cout << "FieldVolume::intersectBound: rayLS=" << rayLS.toString() << std::endl;
	}

	double tnearLS, tfarLS;
	if( m_bboxLS.rayIntersect(rayLS, tnearLS, tfarLS) )
	{
		mint = (m_localToWorld*rayLS(tnearLS) - rayWS.o).norm()*sign(tnearLS);
		maxt = (m_localToWorld*rayLS(tfarLS) - rayWS.o).norm()*sign(tfarLS);


		if(debug)
		{
			std::cout << "FieldVolume::intersectBound: mint=" << mint << std::endl;
			std::cout << "FieldVolume::intersectBound: maxt=" << maxt << std::endl;
			std::cout << "FieldVolume::intersectBound: bboxLS=" << m_bboxLS.toString() << std::endl;
		}


		mint = std::max( mint, rayWS.mint );
		maxt = std::min( maxt, rayWS.maxt );
		if(maxt-mint>0)
			return true;
	}
	return false;
}
V3d Volume::getMaxExtinction()const
{
	return m_extinction_max;
}
std::string Volume::toString()const
{
	std::ostringstream ss;
	ss << "Volume bound(worldspace)=" << m_bboxWS.min.toString() << " " << m_bboxWS.max.toString() << std::endl;
	//ss << "extinction field=" << m_field_extinction->toString() << std::endl;
	//ss << "albedo field=" << m_field_albedo->toString() << std::endl;
	return ss.str();
}
