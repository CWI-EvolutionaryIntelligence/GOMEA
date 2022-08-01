#pragma once

#include <iostream> 
#include <vector>
#include <random>
using namespace std;

class Individual
{
public:
    size_t numberOfVariables;
    size_t alphabetSize;
    vector<char> genotype;
    double fitness;

    Individual() {};

    Individual(size_t numberOfVariables_, size_t alphabetSize_): numberOfVariables(numberOfVariables_), alphabetSize(alphabetSize_)
    {
        genotype.resize(numberOfVariables_);
        fill(genotype.begin(), genotype.end(), 0);
    }

    Individual(vector<char> &genotype_, double fitness_): fitness(fitness_)
    {
        numberOfVariables = genotype_.size();
        genotype.resize(numberOfVariables);
        copy(genotype_.begin(), genotype_.end(), genotype.begin());
    }

    void randomInit(mt19937 *rng)
    {
        for (size_t i = 0; i < numberOfVariables; ++i)
        {
            genotype[i] = (*rng)() % alphabetSize;
        }
    }

    friend ostream & operator << (ostream &out, const Individual &individual);

    Individual& operator=(const Individual& other)
    {
        alphabetSize = other.alphabetSize;
        numberOfVariables = other.numberOfVariables;

        genotype = other.genotype;
        fitness = other.fitness;

        return *this;
    }

    bool operator==(const Individual& solutionB)
    {
        for (size_t i = 0; i < numberOfVariables; ++i)
        {
            if (this->genotype[i] != solutionB.genotype[i])
                return false;
        }
        return true;
    }
};

