#pragma once

#include <iostream> 
#include <vector>
#include <random>
#include "gomea/src/common/solution.hpp"
using namespace std;

namespace gomea{
namespace discrete{

class Individual : public gomea::solution_t<char>
{
public:
    size_t alphabetSize;

    Individual(size_t numberOfVariables_, size_t alphabetSize_ ): solution_t(numberOfVariables_), alphabetSize(alphabetSize_)
    {
        fill(variables.begin(), variables.end(), 0);
    }

    void randomInit(mt19937 *rng)
    {
        for (size_t i = 0; i < getNumberOfVariables(); ++i)
        {
            variables[i] = (*rng)() % alphabetSize;
        }
    }

    friend ostream & operator << (ostream &out, const Individual &individual);

    bool operator==(const Individual& solutionB)
    {
        for (size_t i = 0; i < getNumberOfVariables(); ++i)
        {
            if (this->variables[i] != solutionB.variables[i])
                return false;
        }
        return true;
    }
};

}}