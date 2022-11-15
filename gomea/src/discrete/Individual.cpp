#include "gomea/src/discrete/Individual.hpp"

namespace gomea{
namespace discrete{

ostream & operator << (ostream &out, const Individual &individual)
{
	for (size_t i = 0; i < individual.getNumberOfVariables(); ++i)
		out << +individual.variables[i];
	out << " | " << individual.getObjectiveValue();
	return out;
}

}}