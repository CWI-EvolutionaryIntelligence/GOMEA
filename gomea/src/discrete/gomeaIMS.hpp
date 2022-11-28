#pragma once

#include <vector>
#include <unordered_map>
using namespace std;

#include "gomea/src/discrete/Config.hpp"
#include "gomea/src/discrete/Population.hpp"
#include "gomea/src/discrete/shared.hpp"
#include "gomea/src/discrete/gomea.hpp"
#include "gomea/src/fitness/benchmarks-discrete.hpp"
#include "gomea/src/utils/time.hpp"

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
		currentGOMEAIndex = 0;
	bool isInitialized = false,
        hasTerminated = false;
    time_t start_time;

    Config *config;
    vector<Population*> GOMEAs;
    fitness_t *problemInstance = NULL;
    sharedInformation *sharedInformationInstance = NULL;

	gomeaIMS();
	//gomeaIMS(int _problemIndex, int _numberOfVariables, int _maximumNumberOfGOMEAs, int _IMSsubgenerationFactor, int _basePopulationSize, int _maxArchiveSize, string _folder );
    gomeaIMS(Config *config_);
    ~gomeaIMS();
   
   	void initialize();
    void initializeNewGOMEA();
	void initializeGPU( void );
    bool checkTermination();
	bool checkTimeLimitTerminationCriterion();
    double getProgressUntilTermination();
	void generationalStepAllGOMEAs();
    bool checkTerminationGOMEA(int GOMEAIndex);
    void GOMEAGenerationalStepAllGOMEAsRecursiveFold(int indexSmallest, int indexBiggest);
    void run();
	void runGeneration();
	void runGeneration( int GOMEAIndex );
};

}}