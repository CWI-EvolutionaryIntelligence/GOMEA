#include "discrete_gomea/Individual.hpp"

ostream & operator << (ostream &out, const Individual &individual)
{
	for (size_t i = 0; i < individual.numberOfVariables; ++i)
		out << +individual.genotype[i];
	out << " | " << individual.fitness;
	return out;
}