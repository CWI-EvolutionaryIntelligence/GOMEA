# pragma once

namespace gomea{
namespace realvalued{

class fitness_buffer_t
{
	
};

class double_buffer_t : public fitness_buffer_t
{
	public:
		double val;
		void update( double x );
};

}}
