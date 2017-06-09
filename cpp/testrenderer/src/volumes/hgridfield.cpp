#include <volumes/hgridfield.h>






namespace field
{



	template<>
	V3d HGridField<V3d>::Block::eval( const P3d& pLS, bool debug )const
	{
		return V3d(0.0);
	}

	template<>
	double HGridField<double>::Block::eval( const P3d& pLS, bool debug )const
	{
		P3d pVS = localToVoxel(pLS) - V3d(0.5, 0.5, 0.5);

		// lower left corner
		V3i c1;
		c1[0] = (int)floor(pVS.x());
		c1[1] = (int)floor(pVS.y());
		c1[2] = (int)floor(pVS.z());

		// upper right corner
		V3i c2 = c1+Vector3i(1);

		// clamp the indexing coordinates
		c1[0] = std::max(0, std::min(c1.x(), m_resolution.x()-1));
		c2[0] = std::max(0, std::min(c2.x(), m_resolution.x()-1));
		c1[1] = std::max(0, std::min(c1.y(), m_resolution.y()-1));
		c2[1] = std::max(0, std::min(c2.y(), m_resolution.y()-1));
		c1[2] = std::max(0, std::min(c1.z(), m_resolution.z()-1));
		c2[2] = std::max(0, std::min(c2.z(), m_resolution.z()-1));

		double tx = pVS.x() - c1[0];
		double ty = pVS.y() - c1[1];
		double tz = pVS.z() - c1[2];

		switch (m_dataType)
		{
			case EUInt8:
			{
				const double
				d000 = m_data[(c1[2]*m_resolution.y() + c1[1])*m_resolution.x() + c1[0]]/255.0,
				d001 = m_data[(c1[2]*m_resolution.y() + c1[1])*m_resolution.x() + c2[0]]/255.0,
				d010 = m_data[(c1[2]*m_resolution.y() + c2[1])*m_resolution.x() + c1[0]]/255.0,
				d011 = m_data[(c1[2]*m_resolution.y() + c2[1])*m_resolution.x() + c2[0]]/255.0,
				d100 = m_data[(c2[2]*m_resolution.y() + c1[1])*m_resolution.x() + c1[0]]/255.0,
				d101 = m_data[(c2[2]*m_resolution.y() + c1[1])*m_resolution.x() + c2[0]]/255.0,
				d110 = m_data[(c2[2]*m_resolution.y() + c2[1])*m_resolution.x() + c1[0]]/255.0,
				d111 = m_data[(c2[2]*m_resolution.y() + c2[1])*m_resolution.x() + c2[0]]/255.0;
				return ((d000*(1.0-tx) + d001*tx)*(1.0-ty) +
						(d010*(1.0-tx) + d011*tx)*ty)*(1.0-tz) +
					   ((d100*(1.0-tx) + d101*tx)*(1.0-ty) +
						(d110*(1.0-tx) + d111*tx)*ty)*tz;
			}
		};
		return 0.0;
	}



} // namespace field
