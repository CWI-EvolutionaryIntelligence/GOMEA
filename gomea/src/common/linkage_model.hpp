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

#include "gomea/src/common/linkage_config.hpp"
#include "gomea/src/common/solution.hpp"
#include "gomea/src/fitness/fitness.hpp"
#include "gomea/src/utils/time.hpp"
#include "gomea/src/utils/tools.hpp"

namespace gomea{

class linkage_model_t
{   
public:
    vec_t<vec_t<int>> FOSStructure;  
	size_t numberOfVariables;
	vec_t<vec_t<int>> parallelFOSGroups;
    vec_t<int> FOSorder;
    
    vec_t<int> improvementCounters;
    vec_t<int> usageCounters;

    linkage::linkage_model_type type = linkage::CUSTOM_LM;
    bool is_static = true;
	int maximumSetSize;
    double MI_truncation_factor = 1.0; // Calculate MI on this fraction of the (best solutions in the) population -- 1.0 means all solutions, 0.5 means the 50% best solutions, etc.

    virtual ~linkage_model_t(){};

    static std::shared_ptr<linkage_model_t> createFOSInstance( const linkage_config_t &config, size_t numberOfVariables = 0 );
    static std::shared_ptr<linkage_model_t> createLinkageTreeFOSInstance(size_t FOSIndex, size_t numberOfVariables, linkage::similarity_measure_type similarityMeasure, int maximumFOSSetSize, bool is_static = false );
    static std::string getTypeName( linkage::linkage_model_type type );
    
    static std::shared_ptr<linkage_model_t> univariate(size_t numberOfVariables_);
    static std::shared_ptr<linkage_model_t> full(size_t numberOfVariables_);
    static std::shared_ptr<linkage_model_t> linkage_tree(size_t numberOfVariables_, linkage::similarity_measure_type similarityMeasure_, bool filtered_,  int maximumSetSize_, bool is_static_ );
    static std::shared_ptr<linkage_model_t> marginal_product_model(size_t numberOfVariables_, size_t block_size);
    static std::shared_ptr<linkage_model_t> custom_fos(size_t numberOfVariables_, const vec_t<vec_t<int>> &FOS);
    static std::shared_ptr<linkage_model_t> from_file( std::string filename );

    size_t size()
    {
        return FOSStructure.size();
    }

    size_t elementSize(int i)
    {
        return FOSStructure[i].size();
    }
	
    vec_t<int> getSet( int element_index )
    {
        return FOSStructure[element_index];
    }

    linkage::similarity_measure_type getSimilarityMeasure();

    void addGroup(int var_index);
    void addGroup(const std::set<int> &group);
    virtual void addGroup( vec_t<int> group ); 

    void writeToFileFOS(std::string folder, int populationIndex, int generation);
    void writeFOSStatistics(std::string folder, int populationIndex, int generation);
    void setCountersToZero();
    void shuffleFOS();
	void determineParallelFOSOrder(vec_t<vec_t<int>> VIG );
    void determineParallelFOSOrder(std::map<int, std::set<int>> VIG );
    vec_t<int> graphColoring(std::map<int, std::set<int>> &VIG);
    void initializeDependentSubfunctions( std::map<int,std::set<int>> &subfunction_dependency_map );
    std::set<int> getDependentSubfunctions( int linkage_set_index );

    void writeMIMatrixToFile(vec_t<vec_t<double>> MI_Matrix, std::string folder, int populationIndex, int generation);
    
    void learnLinkageTreeFOS( fitness::fitness_t<char> *problemInstance, vec_t<solution_t<char>*> &population );
	void learnLinkageTreeFOS( vec_t<vec_t<double>> similarity_matrix, bool include_full_fos_element );

    void printFOS();

protected:
    vec_t<int> colors;
    vec_t<vec_t<double>> S_Matrix;
    vec_t<std::set<int>> dependent_subfunctions;
    bool filtered;
    linkage::similarity_measure_type similarityMeasure;
    
    linkage_model_t( size_t numberOfVariables_ ): numberOfVariables(numberOfVariables_), maximumSetSize(numberOfVariables_) {}
    linkage_model_t( size_t numberOfVariables_, size_t block_size ); 
    linkage_model_t( size_t numberOfVariables_, const vec_t<vec_t<int>> &FOS );
    linkage_model_t( size_t numberOfVariables_, linkage::similarity_measure_type similarityMeasure, bool filtered, int maximumSetSize, bool is_static );
    linkage_model_t( std::string filename );
	
    int determineNearestNeighbour(size_t index, const vec_t< vec_t< int > > &mpm );
    vec_t<vec_t<double>> computeMIMatrix(vec_t<solution_t<char>*> &population, size_t alphabetSize);
    vec_t<vec_t<double>> computeNMIMatrix(vec_t<solution_t<char>*> &population, size_t alphabetSize);
    vec_t<vec_t<double>> computeHammingDistanceSimilarityMatrix( vec_t<solution_t<char>*> &population );
    void estimateParametersForSingleBinaryMarginal(vec_t<solution_t<char>*> &population, size_t alphabetSize, vec_t<size_t> &indices, size_t &factorSize, vec_t<double> &result, double truncation_factor = 1.0);
};

typedef std::shared_ptr<linkage_model_t> linkage_model_pt;

}