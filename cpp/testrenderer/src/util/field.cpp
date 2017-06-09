
#include <util/field.h>
#include <math/bbox.h>














namespace field
{




	/*
	void write( const std::string& filename, Field<double>::Ptr field, const Box3d& bound, const V3i& res )
	{
		houio::math::V3i hres( res.x(), res.y(), res.z() );
		houio::math::Box3f hbound( houio::math::V3f(bound.min.x(), bound.min.y(), bound.min.z()),
								   houio::math::V3f(bound.max.x(), bound.max.y(), bound.max.z()));
		houio::ScalarField::Ptr hfield = houio::ScalarField::create( hres, hbound );
		houio::math::V3f voxelSize( (hbound.maxPoint.x-hbound.minPoint.x)/float(res.x()),
									(hbound.maxPoint.y-hbound.minPoint.y)/float(res.y()),
									(hbound.maxPoint.z-hbound.minPoint.z)/float(res.z()));

		for( int z=0;z<res.z();++z )
			for( int y=0;y<res.y();++y )
				for( int x=0;x<res.x();++x )
				{
					houio::math::V3f pWS = hbound.minPoint + houio::math::V3f((x+0.5)*voxelSize.x, (y+0.5)*voxelSize.y, (z+0.5)*voxelSize.z );
					hfield->lvalue( x, y, z ) = float(field->eval( P3d( pWS.x, pWS.y, pWS.z ) ));
				}

		float min, max;
		houio::field_range( *hfield, min, max );
		std::cout << "writing field to " << filename << " min="  << min << " max=" << max << std::endl;

		houio::HouGeoIO::xport( filename, hfield );
	}
	*/

	void write(const std::string& filename, Field<double>* field, const V3i& res, const Transformd &localToWorld, const Box3d &bound_ls)
	{
		houio::math::V3i hres( res.x(), res.y(), res.z() );

		const Eigen::Matrix<double, 4, 4>& transform = localToWorld.getMatrix();
		houio::math::M44f hlocalToWorld( transform(0,0), transform(1,0), transform(2,0), transform(3,0),
										 transform(0,1), transform(1,1), transform(2,1), transform(3,1),
										 transform(0,2), transform(1,2), transform(2,2), transform(3,2),
										 transform(0,3), transform(1,3), transform(2,3), transform(3,3) );


		// voxelsize in local space
		V3d bound_ls_extends = bound_ls.getExtents();
		houio::math::V3f bound_ls_min(bound_ls.min.x(), bound_ls.min.y(), bound_ls.min.z());
		houio::math::V3f voxelSize( float(bound_ls_extends.x())/float(res.x()),
									float(bound_ls_extends.y())/float(res.y()),
									float(bound_ls_extends.z())/float(res.z()));

		houio::math::M44f scale = houio::math::M44f::ScaleMatrix( float(bound_ls_extends.x()),
																  float(bound_ls_extends.y()),
																  float(bound_ls_extends.z()));
		houio::math::M44f offset = houio::math::M44f::TranslationMatrix(bound_ls_min);
		houio::ScalarField::Ptr hfield = houio::ScalarField::create( hres, scale*offset*hlocalToWorld );

		for( int z=0;z<res.z();++z )
			for( int y=0;y<res.y();++y )
				for( int x=0;x<res.x();++x )
				{
					// we create voxel centers in local space (0-1)
					houio::math::V3f pLS = bound_ls_min +houio::math::V3f((x+0.5)*voxelSize.x, (y+0.5)*voxelSize.y, (z+0.5)*voxelSize.z );

					// here we apply the transform to transform the evaluation positions into worldspace
					//hfield->lvalue( x, y, z ) = float(field->eval( localToWorld*P3d( pLS.x, pLS.y, pLS.z ) ));
					hfield->lvalue( x, y, z ) = float(field->eval( P3d( pLS.x, pLS.y, pLS.z ) ));
				}

		float min, max;
		houio::field_range( *hfield, min, max );
		std::cout << "writing field to " << filename << " min="  << min << " max=" << max << " res=" << res.toString() << std::endl;

		houio::HouGeoIO::xport( filename, hfield );
	}

	void write(const std::string& filename, Field<float>* field, const V3i& res, const Transformd &localToWorld, const Box3d &bound_ls)
	{
		houio::math::V3i hres( res.x(), res.y(), res.z() );

		const Eigen::Matrix<double, 4, 4>& transform = localToWorld.getMatrix();
		houio::math::M44f hlocalToWorld( transform(0,0), transform(1,0), transform(2,0), transform(3,0),
										 transform(0,1), transform(1,1), transform(2,1), transform(3,1),
										 transform(0,2), transform(1,2), transform(2,2), transform(3,2),
										 transform(0,3), transform(1,3), transform(2,3), transform(3,3) );


		// voxelsize in local space
		V3d bound_ls_extends = bound_ls.getExtents();
		houio::math::V3f bound_ls_min(bound_ls.min.x(), bound_ls.min.y(), bound_ls.min.z());
		houio::math::V3f voxelSize( float(bound_ls_extends.x())/float(res.x()),
									float(bound_ls_extends.y())/float(res.y()),
									float(bound_ls_extends.z())/float(res.z()));

		houio::math::M44f scale = houio::math::M44f::ScaleMatrix( float(bound_ls_extends.x()),
																  float(bound_ls_extends.y()),
																  float(bound_ls_extends.z()));
		houio::math::M44f offset = houio::math::M44f::TranslationMatrix(bound_ls_min);
		houio::ScalarField::Ptr hfield = houio::ScalarField::create( hres, scale*offset*hlocalToWorld );

		for( int z=0;z<res.z();++z )
			for( int y=0;y<res.y();++y )
				for( int x=0;x<res.x();++x )
				{
					// we create voxel centers in local space (0-1)
					houio::math::V3f pLS = bound_ls_min +houio::math::V3f((x+0.5)*voxelSize.x, (y+0.5)*voxelSize.y, (z+0.5)*voxelSize.z );
					// here we apply the transform to transform the evaluation positions into worldspace
					hfield->lvalue( x, y, z ) = float(field->eval( localToWorld*P3d( pLS.x, pLS.y, pLS.z ) ));
				}

		float min, max;
		houio::field_range( *hfield, min, max );
		std::cout << "writing field to " << filename << " min="  << min << " max=" << max << std::endl;

		houio::HouGeoIO::xport( filename, hfield );
	}
	/*
	void write( const std::string& filename, Field<double>::Ptr field, const Box3d& bound, const V3i& res, bool debug  )
	{
		houio::math::V3i hres( res.x(), res.y(), res.z() );
		houio::math::Box3f hbound( houio::math::V3f(bound.min.x(), bound.min.y(), bound.min.z()),
								   houio::math::V3f(bound.max.x(), bound.max.y(), bound.max.z()));
		houio::ScalarField::Ptr hfield = houio::ScalarField::create( hres, hbound );
		if(debug)
		{
			std::cout << "write:\n";
			std::cout << "voxelSize=" << hfield->getVoxelSize().x << " " << hfield->getVoxelSize().y << " " << hfield->getVoxelSize().z << std::endl;
			std::cout << "res=" << hfield->getResolution().x << " " << hfield->getResolution().x << " " << hfield->getResolution().x << std::endl;
		}

		houio::math::V3f voxelSize( (hbound.maxPoint.x-hbound.minPoint.x)/float(res.x()),
									(hbound.maxPoint.y-hbound.minPoint.y)/float(res.y()),
									(hbound.maxPoint.z-hbound.minPoint.z)/float(res.z()));
		//houio::math::V3f voxelSize = hfield->getVoxelSize();
		houio::Geometry::Ptr geo = houio::Geometry::createPointGeometry();
		houio::Attribute::Ptr pattr = geo->getAttr("P");
		houio::Attribute::Ptr vattr = houio::Attribute::createFloat();

		for( int z=0;z<res.z();++z )
			for( int y=0;y<res.y();++y )
				for( int x=0;x<res.x();++x )
				{
					houio::math::V3f pWS = hbound.minPoint + houio::math::V3f((x+0.5)*voxelSize.x, (y+0.5)*voxelSize.y, (z+0.5)*voxelSize.z );
					float value = float(field->eval( P3d( pWS.x, pWS.y, pWS.z ) ));
					//hfield->lvalue( x, y, z ) =

					int index = pattr->appendElement<houio::math::V3f>(pWS);
					vattr->appendElement<float>(value);
					geo->addPoint(index);
				}

		//float min, max;
		//houio::field_range( *hfield, min, max );
		//std::cout << "writing field to " << filename << " min="  << min << " max=" << max << std::endl;

		geo->setAttr("value", vattr);

		houio::HouGeoIO::xport( filename, geo );
	}
	*/
}
