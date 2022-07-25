#pragma once

#include <vector>
#include <set>
#include <string>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <deque>
#include <random>
#include <memory>
using namespace std;

#include "Individual.hpp"
#include "utils.hpp"
#include "time.hpp"

template<typename T>
struct varMat {
	T *data;
	int *rowStartIndices;
	int numRows;
};

class FOS
{   
public:
    vector<vector<int> > FOSStructure;  
	varMat<int> FOSStructure_d;
	size_t numberOfVariables;
    size_t alphabetSize;
	vector<vector<int>> parallelFOSGroups;
    
    vector<int> improvementCounters;
    vector<int> usageCounters;

    FOS(size_t numberOfVariables_, size_t alphabetSize_): numberOfVariables(numberOfVariables_), alphabetSize(alphabetSize_)
    {}

    virtual ~FOS(){};

    size_t FOSSize()
    {
        return FOSStructure.size();
    }

    size_t FOSElementSize(int i)
    {
        return FOSStructure[i].size();
    }
    
	virtual void learnFOS(vector<Individual*> &population, mt19937 *rng = NULL) = 0;
	virtual void learnFOS( vector<vector<double>> MI_matrix, mt19937 *rng = NULL) = 0;
    void writeToFileFOS(string folder, int populationIndex, int generation);
    void writeFOSStatistics(string folder, int populationIndex, int generation);
    void setCountersToZero();
    void shuffleFOS(vector<int> &indices, mt19937 *rng);
	void determineParallelFOSOrder(vector<int> &indices, vector<vector<int>> VIG, mt19937 *rng);
	vector<int> graphColoring( vector<vector<int>> &VIG );
    
    virtual void writeMIMatrixToFile(vector<vector<double>> MI_Matrix, string folder, int populationIndex, int generation){};

private:
	vector<int> colors;
};


class LTFOS: public FOS
{
private:
    vector<vector<double> > S_Matrix;
    bool filtered;
	int maximumSetSize;
    int similarityMeasure;
    int determineNearestNeighbour(size_t index, vector< vector< int > > &mpm ); 
    vector<vector<double>> computeMIMatrix(vector<Individual*> &population);
    vector<vector<double>> computeNMIMatrix(vector<Individual*> &population);
    void estimateParametersForSingleBinaryMarginal(vector<Individual*> &population, vector<size_t> &indices, size_t  &factorSize, vector<double> &result);

public: 
    LTFOS(size_t numberOfVariables_, size_t alphabetSize_, int similarityMeasure, bool filtered=false, int maximumSetSize = -1);
    ~LTFOS(){};

	void learnFOS(vector<Individual*> &population, mt19937 *rng = NULL);
	void learnFOS( vector<vector<double>> MI_matrix, mt19937 *rng = NULL);
    void writeMIMatrixToFile(vector<vector<double>> MI_Matrix, string folder, int populationIndex, int generation);
};

typedef std::shared_ptr<FOS> FOS_t;
bool FOSNameByIndex(size_t FOSIndex, string &FOSName);
FOS_t createFOSInstance(size_t FOSIndex, size_t numberOfVariables, size_t alphabetSize, int similarityMeasure, int maximumFOSSetSize = -1);

