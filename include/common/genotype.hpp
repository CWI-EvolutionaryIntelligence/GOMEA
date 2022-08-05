#pragma once

#include "common/gomea_defs.hpp"

namespace gomea{

template<class T>
class genotype_t 
{
	public:
		virtual T& operator[](std::size_t idx) = 0;
    	virtual const T& operator[](std::size_t idx) const = 0;
		virtual int getNumberOfVariables() = 0;
};

}
