#include <PNVolume.h>

#include <math/sph.h>









PNVolume::PNVolume(const Domain &domain):
	m_domain(domain),
	m_extinction_minimum_threshold(0.0)
	//m_worldToLocal(),
	//m_localToWorld()//,
	//m_bboxLS( P3d(0.0f,0.0f,0.0f), P3d(1.0f,1.0f,1.0f) ),
	//m_bboxWS( P3d(0.0f,0.0f,0.0f), P3d(1.0f,1.0f,1.0f) ),
	//m_phaseFunction( std::make_shared<IsotropicPhase>() )
{
	setExtinctionAlbedo( std::make_shared<ConstantField3d>(V3d(1.0, 1.0, 1.0)),
						 std::make_shared<ConstantField3d>(V3d(0.9, 0.9, 0.9)) );
}
void PNVolume::setExtinctionMinimumThreshold(double min_threshold)
{
	m_extinction_minimum_threshold = min_threshold;
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
	setEmission(sph::index(l, m), field);
}

void PNVolume::setEmission(int sh_index, Field3d::Ptr field)
{
	if( sh_index>=m_field_q.size() )
		m_field_q.resize(sh_index+1);
	m_field_q[sh_index] = field;
}

void PNVolume::setPhase(int l, int m, Field3d::Ptr field)
{
	setPhase(sph::index(l, m), field);
}

void PNVolume::setPhase( int sh_index, Field3d::Ptr field)
{
	if( sh_index>=m_field_p.size() )
		m_field_p.resize(sh_index+1);
	m_field_p[sh_index] = field;
}


Domain& PNVolume::getDomain()
{
	return m_domain;
}

//const Domain& PNVolume::getDomain()const
//{
//	return m_domain;
//}

V3d PNVolume::evalExtinction( const P3d& pWS, bool debug )const
{
	//return m_field_extinction->eval(m_domain.worldToLocal(pWS));

	///*
	// for testing we apply here a minimum threshold
	V3d ext = m_field_extinction->eval(m_domain.worldToLocal(pWS));

	double threshold = m_extinction_minimum_threshold;
	ext[0] = std::max(ext[0], threshold);
	ext[1] = std::max(ext[1], threshold);
	ext[2] = std::max(ext[2], threshold);
	return ext;
	//*/
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

// this creates the next coarser mipmap level of this problem
// used for multigrid solver
PNVolume::Ptr PNVolume::downsample()
{
	// downsample domain and check if we can actually do it
	V3i res_fine = m_domain.getResolution();

	bool is2D = res_fine[2] == 1;

	// for multigrid, we require the resolution to be even,
	// so that we can do restriction and interpolation straigh forwardly
	if( (res_fine[0]%2!=0)||(res_fine[1]%2!=0)||(!is2D && (res_fine[2]%2!=0)))
		return PNVolume::Ptr();

	V3i res_coarse( res_fine[0]/2, res_fine[1]/2, is2D ? 1:res_fine[2]/2 );

	// we stop if any dimension reduces to less than 2.
	if( (res_coarse[0]<2)||(res_coarse[1]<2)||(!is2D && (res_coarse[2]<2)))
		return PNVolume::Ptr();

	// create the coarse domain
	Domain domain_coarse(   m_domain.getBound().getExtents(),
							res_coarse,
							m_domain.getBound().min );

	PNVolume::Ptr problem_coarse = std::make_shared<PNVolume>(domain_coarse);

	problem_coarse->setExtinctionAlbedo( m_field_extinction->downsample(),
										 m_field_albedo->downsample() );
	for( int i=0;i<m_field_q.size();++i )
		if(m_field_q[i])
			problem_coarse->setEmission(i, m_field_q[i]->downsample());
	for( int i=0;i<m_field_p.size();++i )
		if(m_field_p[i])
			problem_coarse->setPhase(i, m_field_p[i]->downsample());

	return problem_coarse;
}

bool PNVolume::is2D()const
{
	return m_domain.getResolution()[2] == 1;
}
