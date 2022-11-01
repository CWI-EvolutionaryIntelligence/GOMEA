#pragma once

#include <string>
#include <cstdlib>
#include <iostream>
#include <random>
#include <chrono>
#include <cassert>
#include <unistd.h>
#include <getopt.h>

using namespace std;

#include "gomea/src/discrete/problems.hpp"
#include "gomea/src/discrete/FOS.hpp"

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
    size_t problemIndex = 0, k = 1, s = 1,   
        FOSIndex = 0,              
        numberOfVariables = 10;
	int GPUIndex = -1;
    int similarityMeasure = 0;
	int maximumFOSSetSize = -1;

    string folder = "test";
    string problemName,
           FOSName;
    string problemInstancePath = "";

    //long long timelimitMilliseconds = -1,
    long long randomSeed,
			  cudaRandstateOffset;
    
    size_t alphabetSize = 2;
    size_t maxArchiveSize = 1000000;
    int maximumNumberOfGOMEAs  = 100,
        IMSsubgenerationFactor = 4,
        basePopulationSize     = 2;
    
    mt19937 rng;
};
