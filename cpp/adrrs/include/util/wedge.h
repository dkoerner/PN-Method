#pragma once
#include <vector>
#include <map>

#include <util/string.h>


struct Wedge
{
	struct Parm
	{
		Parm( const std::string& name, int index ):
			name(name),
			index(index)
		{
		}


		struct Value
		{
			enum EType
			{
				EFloat,
				EInt
			};

			static Value make_int( int value, int index = -1 )
			{
				Value v;
				v.type = EInt;
				v.i = value;
				v.index = index;
				return v;
			}

			static Value make_float( float value, int index = -1 )
			{
				Value v;
				v.type = EFloat;
				v.f = value;
				v.index = index;
				return v;
			}

			std::string asString()
			{
				switch(type)
				{
				case EFloat:return toString<float>(f);break;
				case EInt:return toString<int>(i);break;
				};
				return "";
			}

			// -------
			EType type;
			union
			{
				float f;
				int i;
			};
			// if a parameter has N values associated, then this index will tell to which one of those N it is
			// this is used for string expansion
			int index;

		};


		std::string name;
		int index;
		std::vector<Value> values;
	};

	struct Iteration
	{
		Iteration()
		{
		}
		Iteration( const Iteration& other ):
			m_parms(other.m_parms),
			m_parm_values(other.m_parm_values)
		{
		}

		void print()
		{
			/*
			for( auto it:m_parm_values )
			{
				int parm_index = it.first;
				std::string parm_name = (*m_parms)[parm_index].name;
				int parm_value = m_parm_values[parm_index];
				std::cout << parm_name << "=" << parm_value << " ";
			}
			std::cout << std::endl;
			*/
		}

		std::string expand_value( const std::string& input )
		{
			std::string result = input;
			for( int i=0;i<m_parms->size();++i )
			{
				Parm::Value* value = getValue(getParm(i));
				result = replace( result, "$"+toString<int>(i), value->asString() );
				//result = replace( result, "$"+toString<int>(i), toString<int>(value->index) );
			}
			return result;
		}

		std::string expand_index( const std::string& input )
		{
			std::string result = input;
			for( int i=0;i<m_parms->size();++i )
			{
				Parm::Value* value = getValue(getParm(i));
				result = replace( result, "$"+toString<int>(i), toString<int>(value->index) );
			}
			return result;
		}

		Parm* getParm(int index)
		{
			return &(*m_parms)[index];
		}

		Parm* getParm(const std::string& parm_name)
		{
			Parm* parm = 0;
			for( auto& parm_item : *m_parms )
			{
				if( parm_item.name == parm_name )
				{
					parm = &parm_item;
					break;
				}
			}
			return parm;
		}

		Parm::Value* getValue(Parm* parm)
		{
			if(!parm)
				throw std::runtime_error("Iteration::getValue: parm does not exist");
			return &m_parm_values[parm->index];
		}

		int getValueIndex( int parameter_index )
		{
			return getValue(getParm(parameter_index))->index;
		}



		int getInt( const std::string& parm_name )
		{
			Parm::Value* v = getValue(getParm(parm_name));
			if(v && v->type == Parm::Value::EInt)
				return v->i;
			throw std::runtime_error("Iteration::getInt: failed to retrieve integer value from parameter.");
			return 0;
		}

		float getFloat( const std::string& parm_name )
		{
			Parm::Value* v = getValue(getParm(parm_name));
			if(v && v->type == Parm::Value::EFloat)
				return v->f;
			throw std::runtime_error("Iteration::getFloat: failed to retrieve integer value from parameter.");
			return 0.0;
		}

		std::map<int, Parm::Value> m_parm_values;
		std::vector<Parm>* m_parms;
	};




	void addParm( const std::string& name, const std::vector<int>& values )
	{
		Parm& parm = addParm(name);
		for( int i=0;i<values.size();++i )
			parm.values.push_back(Parm::Value::make_int(values[i], i));
	}


	void addParm( const std::string& name, const std::vector<float>& values )
	{
		Parm& parm = addParm( name );
		for( int i=0;i<values.size();++i )
			parm.values.push_back(Parm::Value::make_float(values[i], i));
	}


	void build()
	{
		m_iterations.clear();
		curIter = Iteration();
		curIter.m_parms = &m_parms;
		build_parm(0);
	}


	std::vector<Iteration>& iterations()
	{
		return m_iterations;
	}

	std::vector<Iteration> find_iterations( std::vector<std::pair<int, int>>& fixed )
	{
		std::vector<Iteration> result;

		for( auto& iteration:m_iterations )
		{
			int match = true;
			for( auto& f : fixed )
			{
				int parameter_index = f.first;
				int parameter_value_index = f.second;

				Parm* parm = &m_parms[parameter_index];
				Parm::Value* v = iteration.getValue( parm );
				if(!v)
					throw std::runtime_error("WedgeBuilder::find_iterations no value associated for parm");
				if(v->index != parameter_value_index)
				{
					match = false;
					break;
				}
			}

			// add iteration
			if(match)
				result.push_back(iteration);
		}

		return result;
	}

	// indices is a list which contains the value index to select for each parameter
	Iteration* find_iteration( std::vector<int> indices)
	{
		for( auto& it:m_iterations )
		{
			bool match = true;
			for( int i=0;i<indices.size();++i )
			{
				Parm* parm = &m_parms[i];
				Parm::Value* v = it.getValue( parm );
				if(!v)
					throw std::runtime_error("WedgeBuilder::find_iteration no value associated for parm");
				if(v->index != indices[i])
					match = false;
			}

			if(match)
				return &it;
		}
		return 0;
	}


//private:
	void build_parm( int index )
	{
		Parm& parm = m_parms[index];
		bool isLast = index == m_parms.size()-1;

		for( auto& value:parm.values )
		{
			curIter.m_parm_values[index] = value;

			if( isLast )
			{
				// commit iteration
				m_iterations.push_back(Iteration(curIter));
			}else
				// recurse
				build_parm(index+1);
		}
	}

	Parm& addParm( const std::string& name )
	{
		int index = m_parms.size();
		m_parms.push_back(Parm(name, index));
		return m_parms.back();
	}

	Iteration curIter;
	std::vector<Iteration> m_iterations;
	std::vector<Parm> m_parms;
};
