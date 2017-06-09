#pragma once


#include <util/field.h>
#include <util/voxelgrid.h>


namespace field
{


	template<typename T>
	struct HGridField : public Field<T>
	{
		typedef std::shared_ptr<HGridField> Ptr;

		struct Block// : public Field<T>
		{
			typedef std::shared_ptr<Block> Ptr;

			enum EDataType
			{
				EFloat32 = 1,
				EFloat16 = 2,
				EUInt8 = 3,
				EQuantizedDirections = 4
			};

			Block( const std::string& filename )
				//:Field<T>()
			{
				load( filename );
			}


			void load( const std::string& filename )
			{
				std::ifstream in(filename, std::ios::binary);

				char header[3];
				in.read(header, sizeof(char)*3);
				if (header[0] != 'V' || header[1] != 'O' || header[2] != 'L')
					std::cout << "HGridField::Block::load: wrong volume format\n";

				uint8_t version;
				in.read( (char*)&version, sizeof(uint8_t) );

				{
					int type;
					in.read( (char*)&type, sizeof(int) );
					m_dataType = EDataType(type);
				}

				in.read( (char*) &m_resolution, sizeof(int)*3 );

				int channels;
				in.read( (char*)&channels, sizeof(int) );

				Box3f aabb;
				in.read( (char*)&aabb.min, sizeof(float)*3 );
				in.read( (char*)&aabb.max, sizeof(float)*3 );

				m_localToWorld = Transformd::from_aabb(Box3d(aabb.min.cast<double>(), aabb.max.cast<double>()));
				m_worldToLocal = m_localToWorld.inverse();

				int data_size = typeSize(m_dataType)*channels*m_resolution.x()*m_resolution.y()*m_resolution.z();
				m_data.resize( data_size );

				in.read( (char*) m_data.data(), data_size );
			}

			std::string typeToString( EDataType type )const
			{
				switch (type)
				{
				case EFloat32:
					return "EFloat32";
				case EFloat16:
					return "EFloat16";
				case EUInt8:
					return "EUInt8";
				case EQuantizedDirections:
					return "EQuantizedDirections";
				};
				return "unknown EVolumeType";
			}

			int typeSize( EDataType type )const
			{
				switch (type)
				{
				case EFloat32:
					return sizeof(float);
				case EFloat16:
					return 0;
				case EUInt8:
					return sizeof(uint8_t);
				case EQuantizedDirections:
					return 0;
				};
				return 0;
			}

			//virtual T eval( const P3d& p, bool debug = false )const override;
			T eval( const P3d& p, bool debug = false )const;

			//virtual void getValueRange( T& min, T& max )const override
			void getValueRange( T& min, T& max )const
			{
				//TODO
			}

			P3d localToVoxel( const P3d& pLS )const
			{
				return P3d( pLS.x()*m_resolution.x(),
							pLS.y()*m_resolution.y(),
							pLS.z()*m_resolution.z());
			}
			P3d voxelToLocal( const P3d& pVS )const
			{
				return P3d( pVS.x()/m_resolution.x(),
							pVS.y()/m_resolution.y(),
							pVS.z()/m_resolution.z());
			}

			std::vector<uint8_t> m_data;
			V3i                  m_resolution; // resolution of single block
			EDataType            m_dataType;
			Transformd           m_localToWorld;
			Transformd           m_worldToLocal;
		};

		HGridField( const std::string& filename_dictionary, const std::string& filename_blocks ):
			Field<T>()
		{
			load(filename_dictionary, filename_blocks);

			// now lets rasterize all blocks to bgeo files
			/*
			for( int k=0;k<m_resolution.z();++k )
				for( int j=0;j<m_resolution.y();++j )
					for( int i=0;i<m_resolution.x();++i )
					{
						int block_index = (m_resolution.y() * k + j) * m_resolution.x() + i;
						Block* block = m_blocks[block_index];
						if(block)
						{
							std::string block_filename = "scarf_bgeo/scarf_$BX_$BY_$BZ.bgeo";
							block_filename = replace( block_filename, "$BX", zeroPadNumber(i, 3) );
							block_filename = replace( block_filename, "$BY", zeroPadNumber(j, 3) );
							block_filename = replace( block_filename, "$BZ", zeroPadNumber(k, 3) );

							//field::write( block_filename, block, block->m_resolution, block->m_localToWorld );
							field::write( block_filename, this, block->m_resolution, getLocalToWorld(), field::voxelBound(V3i(i, j, k), m_resolution) );
						}
					}
			*/

			// lets rasterize the scarf
			//field::write( "scarf.bgeo", this, m_resolution*20, getLocalToWorld() );


			//TODO
		}

		~HGridField()
		{
			for( auto block : m_blocks )
				delete block;
		}

		//p is given in localspace
		virtual T eval( const P3d& pLS, bool debug = false )const override
		{
			// find block
			P3d pBS = localToBlock(pLS);

			// lower left corner
			P3i block_coord;
			block_coord[0] = (int)floor(pBS.x());
			block_coord[1] = (int)floor(pBS.y());
			block_coord[2] = (int)floor(pBS.z());

			// check if block is within resolution bound
			if( !((block_coord.x() >= 0)&&(block_coord.x() < m_resolution.x())&&
				  (block_coord.y() >= 0)&&(block_coord.y() < m_resolution.y())&&
				  (block_coord.z() >= 0)&&(block_coord.z() < m_resolution.z())))
				return T(0.0);


			P3d block_pLS(  pBS.x() - block_coord[0],
							pBS.y() - block_coord[1],
							pBS.z() - block_coord[2]);

			int block_index = (m_resolution.y() * block_coord.z() + block_coord.y()) * m_resolution.x() + block_coord.x();

			Block* block = m_blocks[block_index];
			if( !block )
				return T(0.0);

			return block->eval(block_pLS, debug);
		}

		virtual void getValueRange( T& min, T& max )const override
		{
			min = 0.0;
			max = 1.0;
		}

		const Transformd& getLocalToWorld()const
		{
			return m_localToWorld;
		}

	private:
		P3d localToBlock( const P3d& pLS )const
		{
			return P3d( pLS.x()*m_resolution.x(),
						pLS.y()*m_resolution.y(),
						pLS.z()*m_resolution.z());
		}
		P3d blockToLocal( const P3d& pBS )const
		{
			return P3d( pBS.x()/m_resolution.x(),
						pBS.y()/m_resolution.y(),
						pBS.z()/m_resolution.z());
		}

		void load( const std::string& filename_dictionary, const std::string& block_filename_template )
		{
			// load dictionary
			{
				// NB: Windows x86/x64 is always little endian
				std::ifstream in(filename_dictionary, std::ios::binary);

				Box3f aabb;
				in.read( (char*)&aabb.min, sizeof(float)*3 );
				in.read( (char*)&aabb.max, sizeof(float)*3 );

				m_localToWorld = Transformd::from_aabb( Box3d(aabb.min.cast<double>(), aabb.max.cast<double>()) );

				in.read( (char*)&m_resolution, sizeof(int)*3 );
				int numBlocks = m_resolution.x()*m_resolution.y()*m_resolution.z();
				m_blocks.resize( numBlocks, 0 );

				while( in.good() )
				{
					V3i block_coord;
					in.read( (char*)&block_coord, sizeof(int)*3 );

					std::string block_filename = block_filename_template;
					block_filename = replace( block_filename, "$BX", zeroPadNumber(block_coord.x(), 3) );
					block_filename = replace( block_filename, "$BY", zeroPadNumber(block_coord.y(), 3) );
					block_filename = replace( block_filename, "$BZ", zeroPadNumber(block_coord.z(), 3) );

					if( !file_exists(block_filename) )
						continue;

					//std::cout << "block:" << block_coord.toString() << std::endl;
					Block* block = new Block(block_filename);

					int block_index = (m_resolution.y() * block_coord.z() + block_coord.y()) * m_resolution.x() + block_coord.x();

					m_blocks[block_index] = block;
				}

				//std::cout << "aabb=" << aabb.min.toString() << " " << aabb.max.toString() << std::endl;
				//std::cout << "aabb extend =" << aabb.getExtents().toString() << std::endl;
				//std::cout << "aabb center=" << aabb.getCenter().toString() << std::endl;
				//std::cout << "res=" << resolution << std::endl;


			}

		}

		V3i                 m_resolution;
		Transformd          m_localToWorld;
		std::vector<Block*> m_blocks;
	};


	template<>
	V3d HGridField<V3d>::Block::eval( const P3d& pBLS, bool debug )const;
	template<>
	double HGridField<double>::Block::eval( const P3d& pBLS, bool debug )const;


	template<typename T>
	typename HGridField<T>::Ptr hgrid( const std::string& filename_dictionary, const std::string& filename_blocks )
	{
		return std::make_shared<HGridField<T>>(filename_dictionary, filename_blocks);
	}


} // namespace field

