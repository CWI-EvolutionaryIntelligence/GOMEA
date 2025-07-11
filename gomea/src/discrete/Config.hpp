#pragma once

#include <string>
#include <cstdlib>
#include <iostream>
#include <random>
#include <chrono>
#include <cassert>
#include <cxxopts.hpp>

#include "gomea/src/common/linkage_model.hpp"
#include "gomea/src/common/linkage_config.hpp"
#include "gomea/src/fitness/fitness.hpp"
#include "gomea/src/fitness/benchmarks-discrete.hpp"

namespace gomea{
namespace discrete{

using gomea::fitness::fitness_t;

class Config
{
	void splitString(const std::string &str, vec_t<std::string> &splitted, char delim);
    bool isNumber(const std::string &s);

public:
	Config();
    virtual ~Config();

    bool parseCommandLine(int argc, char **argv);
    void checkOptions();
    void printOverview();
    
	fitness_t<char> *fitness;
    bool gene_invariant = false, 
        generational_statistics = true,
        generational_solution = false,
        write_generational_population = false,
        verbose = false;
	int usePartialEvaluations              = 1,                  
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
    output_frequency_t output_frequency = output_frequency_t::GEN;

    std::string folder = "discrete_gomea_output";
    //std::string problemName,
    std::string FOSName;
    std::string problemInstancePath = "";

    //long long timelimitMilliseconds = -1,
    bool fix_seed = false;
    long long randomSeed;
    
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
