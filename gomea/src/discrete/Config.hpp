#pragma once

#include <string>
#include <cstdlib>
#include <iostream>
#include <random>
#include <chrono>
#include <cassert>
#include <unistd.h>
#include <cxxopts.hpp>

using namespace std;

#include "gomea/src/common/linkage_model.hpp"
#include "gomea/src/fitness/fitness.hpp"
#include "gomea/src/fitness/benchmarks-discrete.hpp"

namespace gomea{
namespace discrete{

typedef gomea::fitness::fitness_t<char> fitness_t;

class Config
{
	void splitString(const string &str, vector<string> &splitted, char delim);
    bool isNumber(const string &s);

public:
	Config();

    bool parseCommandLine(int argc, char **argv);
    void checkOptions();
    void printUsage();
    void printOverview();
    
	fitness_t *fitness;
	int usePartialEvaluations              = 0,                  
		useParallelGOM		               = 1,                  
		useParallelFOSOrder	               = 0,
		popUpdatesDuringGOM				   = 0,
		fixFOSOrderForPopulation		   = 0,
        AnalyzeFOS                         = 0,
        writeElitists					   = 0,
        saveEvaluations                    = 0,
        useForcedImprovements              = 1,
        printHelp                          = 0,
		maximumNumberOfEvaluations		   = -1,
		maximumNumberOfGenerations		   = -1;
	double maximumNumberOfSeconds = -1;
    double vtr = 1e+308;
    size_t k = 1, s = 1,   
        FOSIndex = 0;
	int GPUIndex = -1;
	int maximumFOSSetSize = -1;

    string folder = "discrete_gomea_output";
    //string problemName,
    string FOSName;
    string problemInstancePath = "";

    //long long timelimitMilliseconds = -1,
    bool fix_seed = false;
    long long randomSeed;
    
    size_t alphabetSize = 2;
    size_t maxArchiveSize = 1000000;
    int maximumNumberOfGOMEAs  = 100,
        IMSsubgenerationFactor = 4,
        basePopulationSize     = 2;
    linkage_config_t *linkage_config;

    private:
        cxxopts::Options options = cxxopts::Options("DiscreteGOMEA", "GOMEA for discrete optimization");
        int problemIndex = 0, numberOfVariables = 10;
};

}}