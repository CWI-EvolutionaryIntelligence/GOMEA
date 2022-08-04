#pragma once

#include <vector>
#include <unordered_map>
using namespace std;

#include "discrete_gomea/Config.hpp"
#include "discrete_gomea/Population.hpp"
#include "discrete_gomea/problems.hpp"
#include "discrete_gomea/shared.hpp"
#include "discrete_gomea/gomea.hpp"

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
	bool isInitialized = false;

    Config *config;
    vector<Population*> GOMEAs;
    Problem *problemInstance = NULL;
    sharedInformation *sharedInformationInstance = NULL;

	gomeaIMS();
	gomeaIMS(int _problemIndex, int _numberOfVariables, int _maximumNumberOfGOMEAs, int _IMSsubgenerationFactor, int _basePopulationSize, int _maxArchiveSize, string _folder );
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
