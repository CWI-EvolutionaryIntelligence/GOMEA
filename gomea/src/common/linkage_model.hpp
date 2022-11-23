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

#include "gomea/src/common/solution.hpp"
#include "gomea/src/utils/time.hpp"
#include "gomea/src/utils/tools.hpp"

namespace gomea{

typedef enum{
	UNIVARIATE,
    FULL,
    MPM,
	LINKAGE_TREE,
    CUSTOM_LM
} linkage_model_type;

class linkage_model_t
{   
public:
    vec_t<vec_t<int>> FOSStructure;  
	size_t numberOfVariables;
	vec_t<vec_t<int>> parallelFOSGroups;
    vec_t<int> FOSorder;
    
    vec_t<int> improvementCounters;
    vec_t<int> usageCounters;

    linkage_model_type type = CUSTOM_LM;

    static std::shared_ptr<linkage_model_t> univariate(size_t numberOfvariables_);
    static std::shared_ptr<linkage_model_t> linkage_tree(size_t numberOfVariables_, int similarityMeasure_ = 0, bool filtered_ = false, int maximumSetSize_ = -1);
    static std::shared_ptr<linkage_model_t> conditional(size_t numberOfvariables_);
    static std::shared_ptr<linkage_model_t> marginal_product_model(size_t numberOfvariables_, size_t block_size);
    static std::shared_ptr<linkage_model_t> custom_fos(size_t numberOfvariables_, const vec_t<vec_t<int>> &FOS);
    static std::shared_ptr<linkage_model_t> from_file( FILE *file );

    linkage_model_t( linkage_model_type type, size_t numberOfVariables_ ); 
    virtual ~linkage_model_t(){};

    size_t size()
    {
        return FOSStructure.size();
    }

    size_t elementSize(int i)
    {
        return FOSStructure[i].size();
    }

    void addGroup(int var_index);
    void addGroup(const std::set<int> &group);
    virtual void addGroup( vec_t<int> group ); 

    void writeToFileFOS(string folder, int populationIndex, int generation);
    void writeFOSStatistics(string folder, int populationIndex, int generation);
    void setCountersToZero();
    void shuffleFOS();
	void determineParallelFOSOrder(vec_t<vec_t<int>> VIG, mt19937 *rng);
    void determineParallelFOSOrder(std::map<int, std::set<int>> VIG, mt19937 *rng);
    vec_t<int> graphColoring(std::map<int, std::set<int>> &VIG);

    void writeMIMatrixToFile(vec_t<vec_t<double>> MI_Matrix, string folder, int populationIndex, int generation);
    
    void learnLinkageTreeFOS( vec_t<solution_t<char>*> &population, size_t alphabetSize, mt19937 *rng = NULL );
	void learnLinkageTreeFOS( vec_t<vec_t<double>> MI_matrix, mt19937 *rng = NULL);

protected:
    vec_t<int> colors;
    vec_t<vec_t<double>> S_Matrix;
    bool filtered;
	int maximumSetSize;
    int similarityMeasure;
    
    linkage_model_t(size_t numberOfVariables_): numberOfVariables(numberOfVariables_){}
    linkage_model_t( size_t numberOfVariables_, size_t block_size ); 
    linkage_model_t( size_t numberOfVariables_, const vec_t<vec_t<int>> &FOS );
    linkage_model_t(size_t numberOfVariables_, int similarityMeasure, bool filtered=false, int maximumSetSize = -1);
    linkage_model_t( FILE *file );
	
    int determineNearestNeighbour(size_t index, vec_t< vec_t< int > > &mpm ); 
    vec_t<vec_t<double>> computeMIMatrix(vec_t<solution_t<char>*> &population, size_t alphabetSize);
    vec_t<vec_t<double>> computeNMIMatrix(vec_t<solution_t<char>*> &population, size_t alphabetSize);
    void estimateParametersForSingleBinaryMarginal(vec_t<solution_t<char>*> &population, size_t alphabetSize, vec_t<size_t> &indices, size_t &factorSize, vec_t<double> &result);
};

typedef std::shared_ptr<linkage_model_t> linkage_model_pt;
bool FOSNameByIndex(size_t FOSIndex, string &FOSName);

linkage_model_pt createFOSInstance(size_t FOSIndex, size_t numberOfVariables, int similarityMeasure, int maximumFOSSetSize = -1);

}