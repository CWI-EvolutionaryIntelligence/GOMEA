#pragma once

#include <vector>
#include <unordered_map>
#include <string>
#include <fstream>
#include <deque>
#include <algorithm>
#include <unordered_set>
#include <set>
#include <cassert>
using namespace std;

#include "gomea/src/discrete/Individual.hpp"
#include "gomea/src/discrete/utils.hpp"
#include "gomea/src/discrete/time.hpp"
#include "gomea/src/discrete/FOS.hpp"

class Config;
#include "gomea/src/discrete/Config.hpp"

class Problem
{
public:
    int numberOfVariables;
    int usePartialEvaluations;
    
    Problem(){};
    virtual ~Problem(){};
    virtual void initializeProblem(int numberOfVariables)=0;
    virtual double calculateFitness(Individual *solution)=0;
    virtual double calculateFitnessPartialEvaluations(Individual *solution, Individual *solutionBefore, vector<int> &touchedGenes, double fitnessBefore);
	virtual vector<vector<double>> getMIMatrix();
	virtual vector<vector<int>> getVIG();
};

class oneMax:public Problem
{
public:
    oneMax(){}
    void initializeProblem(int numberOfVariables_)
    {
        numberOfVariables = numberOfVariables_;
    };
    double calculateFitness(Individual *solution);
    double calculateFitnessPartialEvaluations(Individual *solution, Individual *solutionBefore, vector<int> &touchedGenes, double fitnessBefore);
};

class concatenatedDeceptiveTrap:public Problem
{
    int k, s;
    bool bimodal;
    vector<vector<int> > trapsForVariable; //determines traps which each variables belongs to
public:
    concatenatedDeceptiveTrap(int k_, int s_, bool bimodal_): k(k_), s(s_), bimodal(bimodal_)
    {
        if (not bimodal_)
		{}
		else
        {
            if (k != 10)
            {
                cout << "Bimodal trap with k=" << k << " not implemented!" << endl;
                exit(0);
            }
        }
    }
    void initializeProblem(int numberOfVariables_);
    double calculateFitness(Individual *solution);
    double calculateFitnessPartialEvaluations(Individual *solution, Individual *solutionBefore, vector<int> &touchedGenes, double fitnessBefore);
};

struct NKSubfunction
{
    vector<int> variablesPositions;
    vector<double> valuesTable;
};

class ADF:public Problem
{
    string problemInstanceFilename;
    vector<NKSubfunction> subfunctions;
    vector<vector<int> > subfunctionsForVariable;

public:
    ADF(string problemInstanceFilename_): problemInstanceFilename(problemInstanceFilename_)
    {
    }

    void initializeProblem(int numberOfVariables_);

    double calculateFitness(Individual *solution);
    double calculateFitnessPartialEvaluations(Individual *solution, Individual *solutionBefore, vector<int> &touchedGenes, double fitnessBefore);
};

class hierarchialDeceptiveTrap:public Problem
{
    int k;
    vector<int> transform;
public:
    hierarchialDeceptiveTrap(int k_): k(k_)
    {
    }

    void initializeProblem(int numberOfVariables_)
    {
        numberOfVariables = numberOfVariables_;
        if (!isPowerOfK(numberOfVariables, k))
        {
            cerr << "Number of bits should be a power of k! " << numberOfVariables << " is not a power of " << k << endl;
            exit(0);
        }
        transform.resize(numberOfVariables);        
    };

    double generalTrap(int unitation, double leftPeak, double rightPeak);
    double calculateFitness(Individual *solution);
    //double calculateFitnessPartialEvaluations(Individual *solution, Individual *solutionBefore, vector<double> &touchedGenes, double fitnessBefore);
};

class hierarchialIfAndOnlyIf:public Problem
{
public:
    hierarchialIfAndOnlyIf()
    {
    }
    void initializeProblem(int numberOfVariables_)
    {
        numberOfVariables = numberOfVariables_;
        if (!isPowerOfK(numberOfVariables, 2))
        {
            cerr << "Number of bits should be a power of 2! " << numberOfVariables<< " is not a power of 2" << endl;
            exit(0);
        }
    };

    double calculateFitness(Individual *solution);
};

class maxCut:public Problem
{
    string problemInstanceFilename;
    vector<pair<pair<int, int>, double > > edges;
    vector<vector<int> > edgesForVariable;
    vector<vector<int> > neighborsForVariable;

public:
    maxCut(string problemInstanceFilename_): problemInstanceFilename(problemInstanceFilename_)
    {
    }
    void initializeProblem(int numberOfVariables_);
    double calculateFitness(Individual *solution);
    double calculateFitnessPartialEvaluations(Individual *solution, Individual *solutionBefore, vector<int> &touchedGenes, double fitnessBefore);
	vector<vector<double>> getMIMatrix();
	vector<vector<int>> getVIG();

};

class leadingOnes:public Problem
{
public:
    leadingOnes()
    {
    }
    void initializeProblem(int numberOfVariables_)
    {
        numberOfVariables = numberOfVariables_;
    };

    double calculateFitness(Individual *solution);
    //double calculateFitnessPartialEvaluations(Individual *solution, Individual *solutionBefore, vector<double> &touchedGenes, double fitnessBefore);
};

class Clustering:public Problem
{
    string problemInstanceFilename;
    vector<vector<double> > points;
    vector<vector<double> > distances;
    
    int Dim;
public:
    Clustering(string problemInstanceFilename_): problemInstanceFilename(problemInstanceFilename_)
    {
    }
    void initializeProblem(int numberOfVariables_);
    double calculateFitness(Individual *solution);
};


double deceptiveTrap(int unitation, int k);
double bimodalDeceptiveTrap(int unitation, int k);

void createProblemInstance(int problemIndex, int numberOfVariables, Config *config, Problem **problemInstance, string &instancePath, int k = 1, int s = 1);
bool problemNameByIndex(Config *config, string &problemName);

