#pragma once
#include <vector>
#include <map>
#include <string>

#include <util/string.h>











struct Wedge
{
	struct Value
	{
		enum EType
		{
			EFloat,
			EInt
		};

		static Value make_int( int value )
		{
			Value v;
			v.type = EInt;
			v.i = value;
			return v;
		}

		static Value make_float( float value )
		{
			Value v;
			v.type = EFloat;
			v.f = value;
			return v;
		}


		std::string toString()
		{
			switch(type)
			{
			case EFloat:return ::toString<float>(f);break;
			case EInt:return ::toString<int>(i);break;
			};
			return "";
		}

		int asInt()
		{
			return i;
		}


		float asFloat()
		{
			return f;
		}

		// -------
		EType type;
		union
		{
			float f;
			int i;
		};
	};

	struct Iterator
	{
		Iterator( int currentIteration, Wedge* wedge):
			m_currentIteration(currentIteration),
			m_indices(wedge->getNumParameters(), 0),
			m_wedge(wedge)
		{
		}

		int getIndex( const std::string& name )
		{
			return m_indices[m_wedge->m_parameter_name_to_index[name]];
		}

		Value& getValue( const std::string& name )
		{
			int parm_index = m_wedge->getParameterIndex(name);
			int parm_value_index = m_indices[parm_index];
			return m_wedge->m_parameter_values[parm_index][parm_value_index];
		}

		int getInt( const std::string& name )
		{
			return getValue(name).asInt();
		}

		float getFloat( const std::string& name )
		{
			return getValue(name).asFloat();
		}

		void advance()
		{
			++m_currentIteration;

			int increment = 1;
			for( int i=0;i<m_indices.size();++i )
			{
				m_indices[i] += increment;

				// check if index is smaller than dimension
				if(m_indices[i] < m_wedge->m_parameter_value_count[i])
					return;

				// otherwise we need to handle carry
				increment = m_indices[i]/m_wedge->m_parameter_value_count[i];
				m_indices[i] = m_indices[i] % m_wedge->m_parameter_value_count[i];
			}

			return;
		}


		Iterator& operator ++ ()
		{
			advance();
			return *this;
		}
		Iterator operator ++ (int)
		{
			Iterator it = *this;
			it.advance();
			return it;
		}

		bool operator==(const Iterator& other)
		{
			int numIndices = m_indices.size();
			if( other.m_indices.size() != numIndices)
				return false;
			return m_currentIteration == other.m_currentIteration;
		}

		bool operator!=(const Iterator& other)
		{
			int numIndices = m_indices.size();
			if( other.m_indices.size() != numIndices)
				return true;
			return m_currentIteration != other.m_currentIteration;
		}

		std::string toString()
		{
			std::string result = "";
			for( auto& index:m_indices )
			{
				result += ::toString(index) + " ";
			}
			return result;
		}


		std::string expand_index( const std::string& input )
		{
			std::string result = input;
			int numParms = m_wedge->getNumParameters();
			for( int i=0;i<numParms;++i )
			{
				int currentIndex = m_indices[i];
				result = replace( result, "$"+::toString<int>(i), ::toString<int>(currentIndex) );
			}
			return result;
		}


		int m_currentIteration;
		std::vector<int> m_indices;
		Wedge* m_wedge;
	};

	Iterator begin()
	{
		return Iterator(0, this);
	}

	Iterator end()
	{
		int numIterations = 0;
		if(m_parameter_value_count.size() > 0)
		{
			numIterations = 1;
			for( auto c:m_parameter_value_count )
				numIterations *= c;
		}

		return Iterator(numIterations, this);
	}


	int getNumParameters()
	{
		return m_parameter_names.size();
	}

	int getParameterIndex( const std::string& name )
	{
		return m_parameter_name_to_index[name];
	}


	void addParm( const std::string& name, const std::vector<int>& values )
	{
		int parm_index = m_parameter_names.size();
		m_parameter_names.push_back(name);
		m_parameter_value_count.push_back(values.size());
		m_parameter_values[parm_index] = std::vector<Value>();
		m_parameter_name_to_index[name] = parm_index;
		std::vector<Value>& parameter_values = m_parameter_values[parm_index];
		for( int i=0;i<values.size();++i )
			parameter_values.push_back(Value::make_int(values[i]));
	}


	void addParm( const std::string& name, const std::vector<float>& values )
	{
		int parm_index = m_parameter_names.size();
		m_parameter_names.push_back(name);
		m_parameter_value_count.push_back(values.size());
		m_parameter_values[parm_index] = std::vector<Value>();
		m_parameter_name_to_index[name] = parm_index;
		std::vector<Value>& parameter_values = m_parameter_values[parm_index];
		for( int i=0;i<values.size();++i )
			parameter_values.push_back(Value::make_float(values[i]));
	}



	std::vector<std::string> m_parameter_names;
	std::vector<int> m_parameter_value_count;
	std::map<int, std::vector<Value>> m_parameter_values;
	std::map<std::string, int> m_parameter_name_to_index;


};











/*
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

			int asInt()
			{
				return i;
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
//			for( auto it:m_parm_values )
//			{
//				int parm_index = it.first;
//				std::string parm_name = (*m_parms)[parm_index].name;
//				int parm_value = m_parm_values[parm_index];
//				std::cout << parm_name << "=" << parm_value << " ";
//			}
//			std::cout << std::endl;
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


	int getInt( const std::string& name, int value_index )
	{
		for( auto& p: m_parms)
			if(p.name == name)
				return p.values[value_index].asInt();
		throw std::runtime_error("getInt");
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

*/
