#pragma once

#include <cmath>
#include <iostream> 
#include <vector>

#include "gomea/src/discrete/Config.hpp"
#include "gomea/src/discrete/shared.hpp"
#include "gomea/src/fitness/fitness.hpp"
#include "gomea/src/fitness/benchmarks-discrete.hpp"
#include "gomea/src/common/solution.hpp"
#include "gomea/src/common/partial_solution.hpp"
#include "gomea/src/common/linkage_model.hpp"

namespace gomea{
namespace discrete{

class Population
{
public:
    Config *config;
    fitness_pt<char> problemInstance;
    sharedInformation *sharedInformationPointer;
    size_t GOMEAIndex;
    size_t populationSize;

    vec_t<solution_t<char>*> population;
    vec_t<solution_t<char>*> offspringPopulation;
    vec_t<int> noImprovementStretches;

    bool terminated;
    double averageFitness;
    size_t numberOfGenerations;
    
    linkage_model_pt FOSInstance = NULL;

    Population(Config *config_, fitness_pt<char> problemInstance_, sharedInformation *sharedInformationPointer_, size_t GOMEAIndex_, size_t populationSize_, linkage_model_pt FOSInstance_ = NULL );
    ~Population();

    friend std::ostream & operator << (std::ostream &out, const Population &populationInstance);

    void calculateAverageFitness();
    double getFitnessMean();
    double getFitnessVariance();
    double getConstraintValueMean();
    double getConstraintValueVariance();
    solution_t<char> *getBestSolution();
    solution_t<char> *getWorstSolution();
    bool allSolutionsAreEqual();
    void makeOffspring();
    void copyOffspringToPopulation();
    void generateOffspring();
    void evaluateSolution(solution_t<char> *solution);
    void evaluateSolution(solution_t<char> *solution, solution_t<char> *solutionBefore, vec_t<int> &touchedGenes, double fitnessBefore);
    bool GOM(size_t offspringIndex);
    bool FI(size_t offspringIndex);
    void updateElitistAndCheckVTR(solution_t<char> *solution);
    void checkTimeLimit();
};

}}