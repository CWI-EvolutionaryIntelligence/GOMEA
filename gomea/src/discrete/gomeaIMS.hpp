#pragma once

#include <vector>

#include "gomea/src/discrete/Config.hpp"
#include "gomea/src/discrete/Population.hpp"
#include "gomea/src/discrete/shared.hpp"
#include "gomea/src/discrete/gomea.hpp"
#include "gomea/src/fitness/benchmarks-discrete.hpp"
#include "gomea/src/utils/time.hpp"
#include "gomea/src/common/output_statistics.hpp"

namespace gomea{
namespace discrete{

class gomeaIMS: public GOMEA
{
public:
    int maximumNumberOfGOMEAs;
    int IMSsubgenerationFactor,
		basePopulationSize,
		numberOfGOMEAs = 0,
		numberOfGenerationsIMS = 0,
		minimumGOMEAIndex = 0,
		currentGOMEAIndex = 0,
        numberOfStatisticsWrites = 0;
	bool isInitialized = false,
        hasTerminated = false;
    vec_t<Population*> GOMEAs;
    fitness_t<char> *problemInstance = NULL;
    sharedInformation *sharedInformationInstance = NULL;

	gomeaIMS();
	//gomeaIMS(int _problemIndex, int _numberOfVariables, int _maximumNumberOfGOMEAs, int _IMSsubgenerationFactor, int _basePopulationSize, int _maxArchiveSize, string _folder );
    gomeaIMS(Config *config_);
    ~gomeaIMS();
   
   	void initialize();
    void initializeNewGOMEA();
	void initializeGPU( void );
    bool checkTermination();
    bool checkEvaluationLimitTerminationCriterion();
	bool checkTimeLimitTerminationCriterion();
    double getProgressUntilTermination();
	void generationalStepAllGOMEAs();
    bool checkTerminationGOMEA(int GOMEAIndex);
    void GOMEAGenerationalStepAllGOMEAsRecursiveFold(int indexSmallest, int indexBiggest);
    void run();
	void runGeneration();
	void runGeneration( int GOMEAIndex );
    void writeStatistics( Population *population );
    void writeStatistics( int population_index );
    void ezilaitini();
};

}}