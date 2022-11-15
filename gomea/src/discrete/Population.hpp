#pragma once

#include <cmath>
#include <iostream> 
#include <vector>
using namespace std;

#include "gomea/src/discrete/Individual.hpp"
#include "gomea/src/discrete/Config.hpp"
#include "gomea/src/discrete/shared.hpp"
#include "gomea/src/discrete/FOS.hpp"
#include "gomea/src/fitness/fitness.hpp"
#include "gomea/src/fitness/benchmarks-discrete.hpp"
#include "gomea/src/common/solution.hpp"
#include "gomea/src/common/partial_solution.hpp"

namespace gomea{
namespace discrete{

class Population
{
public:
    Config *config;
    fitness_t *problemInstance;
    sharedInformation *sharedInformationPointer;
    size_t GOMEAIndex;
    size_t populationSize;

    vector<Individual*> population;
    vector<Individual*> offspringPopulation;
    vector<int> noImprovementStretches;

    bool terminated;
    double averageFitness;
    size_t numberOfGenerations;
    
    FOS_t FOSInstance = NULL;

    Population(Config *config_, fitness_t *problemInstance_, sharedInformation *sharedInformationPointer_, size_t GOMEAIndex_, size_t populationSize_, FOS_t FOSInstance_ = NULL );
    ~Population();


    friend ostream & operator << (ostream &out, const Individual &individual);

    void calculateAverageFitness();
    bool allSolutionsAreEqual();
    void makeOffspring();
    void copyOffspringToPopulation();
    void generateOffspring();
    void evaluateSolution(Individual *solution);
    void evaluateSolution(Individual *solution, Individual *solutionBefore, vector<int> &touchedGenes, double fitnessBefore);
    bool GOM(size_t offspringIndex, vector<int> FOSOrder );
    bool FI(size_t offspringIndex, vector<int> FOSOrder );
    void updateElitistAndCheckVTR(Individual *solution);
    void checkTimeLimit();
};

}}