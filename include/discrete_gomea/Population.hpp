#pragma once

#include <cmath>
#include <iostream> 
#include <vector>
using namespace std;

#include "Individual.hpp"
#include "Config.hpp"
#include "shared.hpp"
#include "problems.hpp"
#include "FOS.hpp"

class Population
{
public:
    Config *config;
    Problem *problemInstance;
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

    Population(Config *config_, Problem *problemInstance_, sharedInformation *sharedInformationPointer_, size_t GOMEAIndex_, size_t populationSize_, FOS_t FOSInstance_ = NULL );
    ~Population();


    friend ostream & operator << (ostream &out, const Individual &individual);

    void calculateAverageFitness(); 
    void makeOffspring();
    void copyOffspringToPopulation();
    void generateOffspring();
    void evaluateSolution(Individual *solution);
    void evaluateSolution(Individual *solution, Individual *solutionBefore, vector<int> &touchedGenes, double fitnessBefore);
    bool GOM(size_t offspringIndex, Individual *backup, vector<int> FOSOrder );
    bool FI(size_t offspringIndex, Individual *backup, vector<int> FOSOrder );
    void updateElitistAndCheckVTR(Individual *solution);
    void checkTimeLimit();
};
