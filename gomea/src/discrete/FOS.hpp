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

#include "gomea/src/discrete/Individual.hpp"
#include "gomea/src/discrete/utils.hpp"
#include "gomea/src/utils/time.hpp"

namespace gomea{
namespace discrete{

template<typename T>
struct varMat {
	T *data;
	int *rowStartIndices;
	int numRows;
};

class FOS
{   
public:
    vec_t<vec_t<int> > FOSStructure;  
	varMat<int> FOSStructure_d;
	size_t numberOfVariables;
    size_t alphabetSize;
	vec_t<vec_t<int>> parallelFOSGroups;
    
    vec_t<int> improvementCounters;
    vec_t<int> usageCounters;

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
    
	virtual void learnFOS(vec_t<Individual*> &population, mt19937 *rng = NULL) = 0;
	virtual void learnFOS( vec_t<vec_t<double>> MI_matrix, mt19937 *rng = NULL) = 0;
    void writeToFileFOS(string folder, int populationIndex, int generation);
    void writeFOSStatistics(string folder, int populationIndex, int generation);
    void setCountersToZero();
    void shuffleFOS(vec_t<int> &indices, mt19937 *rng);
	void determineParallelFOSOrder(vec_t<int> &indices, vec_t<vec_t<int>> VIG, mt19937 *rng);
    void determineParallelFOSOrder(vec_t<int> &indices, std::map<int, std::set<int>> VIG, mt19937 *rng);
    vec_t<int> graphColoring(std::map<int, std::set<int>> &VIG);

    virtual void writeMIMatrixToFile(vec_t<vec_t<double>> MI_Matrix, string folder, int populationIndex, int generation){};

private:
	vec_t<int> colors;
};


class LTFOS: public FOS
{
private:
    vec_t<vec_t<double> > S_Matrix;
    bool filtered;
	int maximumSetSize;
    int similarityMeasure;
    int determineNearestNeighbour(size_t index, vec_t< vec_t< int > > &mpm ); 
    vec_t<vec_t<double>> computeMIMatrix(vec_t<Individual*> &population);
    vec_t<vec_t<double>> computeNMIMatrix(vec_t<Individual*> &population);
    void estimateParametersForSingleBinaryMarginal(vec_t<Individual*> &population, vec_t<size_t> &indices, size_t  &factorSize, vec_t<double> &result);

public: 
    LTFOS(size_t numberOfVariables_, size_t alphabetSize_, int similarityMeasure, bool filtered=false, int maximumSetSize = -1);
    ~LTFOS(){};

	void learnFOS(vec_t<Individual*> &population, mt19937 *rng = NULL);
	void learnFOS( vec_t<vec_t<double>> MI_matrix, mt19937 *rng = NULL);
    void writeMIMatrixToFile(vec_t<vec_t<double>> MI_Matrix, string folder, int populationIndex, int generation);
};

typedef std::shared_ptr<FOS> FOS_t;
bool FOSNameByIndex(size_t FOSIndex, string &FOSName);
FOS_t createFOSInstance(size_t FOSIndex, size_t numberOfVariables, size_t alphabetSize, int similarityMeasure, int maximumFOSSetSize = -1);

}}