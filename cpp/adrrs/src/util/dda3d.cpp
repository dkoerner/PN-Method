#include <util/dda3d.h>

bool DDA3D::initialize( const V3i& resolution, const Transformd& localToWorld, const Ray3d& rayWS, bool debug)
{
	res = resolution;

	// intersect with bounding box
	Ray3d rayLS = localToWorld.inverse() * rayWS;

	BoundingBox3d bboxLS;
	bboxLS.reset();
	bboxLS.min = P3d(0.0f,0.0f,0.0f);
	bboxLS.max = P3d(1.0f,1.0f,1.0f);

	double tnearLS, tfarLS;
	double maxt;
	P3d minP, maxP, minP_ls, maxP_ls;
	ratio = 1.0;
	if( bboxLS.rayIntersect(rayLS, tnearLS, tfarLS) )
	{
		// here we add Epsilon in order to make sure that the intersection point is within its first voxel
		minP_ls = rayLS(tnearLS+Epsilon);
		maxP_ls = rayLS(tfarLS);
		minP = localToWorld*minP_ls;
		maxP = localToWorld*maxP_ls;
		m_tstart = (minP - rayWS.o).norm()*sign(tnearLS);
		maxt = (maxP - rayWS.o).norm()*sign(tfarLS);

		m_tstart = std::max( m_tstart, rayWS.mint );
		maxt = std::min( maxt, rayWS.maxt );
		if(maxt-m_tstart<=0)
		{
			if(debug)
				std::cout << "DDA3D::initialize maxt-m_tstart<=0\n";

			// no intersection
			return false;
		}

		ratio = (maxt-m_tstart)/(tfarLS-tnearLS+Epsilon);
		m_tend = maxt;
	}else
	{
		if(debug)
			std::cout << "DDA3D::initialize ray does not intersect grid bound\n";
		return false;
	}


	// convert localspace position of entry point to voxelspace
	vox.x() = std::min(res.x(), std::max(int(minP_ls.x()*res.x()), 0));
	vox.y() = std::min(res.y(), std::max(int(minP_ls.y()*res.y()), 0));
	vox.z() = std::min(res.z(), std::max(int(minP_ls.z()*res.z()), 0));

	for( int i=0;i<3;++i )
	{
		if( rayLS.d[i] > 0.0 )
		{
			m_step[i] = 1;
			end[i] = res[i];
		}else
		{
			m_step[i] = -1;
			end[i] = -1;
		}
	}

	// delta
	delta.x() = m_step.x()/rayLS.d.x();
	delta.y() = m_step.y()/rayLS.d.y();
	delta.z() = m_step.z()/rayLS.d.z();

	// tmax
	tmax.x()=diff_distance( minP_ls.x()*res.x(), rayLS.d.x() );
	tmax.y()=diff_distance( minP_ls.y()*res.y(), rayLS.d.y() );
	tmax.z()=diff_distance( minP_ls.z()*res.z(), rayLS.d.z() );

	if(std::isinf(delta.x()))
		tmax.x() = std::numeric_limits<double>::infinity();
	if(std::isinf(delta.y()))
		tmax.y() = std::numeric_limits<double>::infinity();
	if(std::isinf(delta.z()))
		tmax.z() = std::numeric_limits<double>::infinity();
	if( debug )
	{
		std::cout << "DDA3D::initialize delta=" << delta.toString() << std::endl;
		std::cout << "DDA3D::initialize tmax=" << tmax.toString() << std::endl;
		std::cout << "DDA3D::initialize rayLS.d=" << rayLS.d.toString() << std::endl;
	}

	m_lastmax = m_tstart;

	if( tmax.x() < tmax.y() )
	{
		if( tmax.x() < tmax.z() )
			axis = 0;
		else
			axis = 2;
	}else
	{
		if( tmax.y() < tmax.z() )
			axis = 1;
		else
			axis = 2;
	}

	return true;
}

bool DDA3D::step( VoxelStep& info, bool debug )
{
	// check if we have left the grid
	if(vox[axis] == end[axis])
		return false;

	// do next step
	if( tmax.x() < tmax.y() )
	{
		if( tmax.x() < tmax.z() )
			axis = 0;
		else
			axis = 2;
	}else
	{
		if( tmax.y() < tmax.z() )
			axis = 1;
		else
			axis = 2;
	}

	info.min = m_lastmax;
	info.max = m_tstart + tmax[axis]*ratio/res[axis];
	m_lastmax = info.max;
	if(debug)
	{
		std::cout << "DDA3D info.min=" << info.min << std::endl;
		std::cout << "DDA3D m_lastmax=" << m_lastmax << std::endl;
		std::cout << "DDA3D m_tstart=" << m_tstart << std::endl;
		std::cout << "DDA3D axis=" << axis << std::endl;
		std::cout << "DDA3D tmax[axis]=" << tmax[axis] << std::endl;
		std::cout << "DDA3D ratio=" << ratio << std::endl;
		std::cout << "DDA3D res[axis]=" << res[axis] << std::endl;
	}
	info.coord = vox;

	// update
	vox[axis] += m_step[axis];
	tmax[axis] += delta[axis];
	return true;
}

