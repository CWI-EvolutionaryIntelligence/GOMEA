#include "gomea/src/discrete/Config.hpp"

namespace gomea{
namespace discrete{

Config::Config(){}

/**
 * Parses the command line.
 * For options, see printUsage.
 */
void Config::splitString(const std::string &str, vec_t<std::string> &splitted, char delim)
{
    size_t current, previous = 0;
    current = str.find(delim);
    while (current != std::string::npos)
    {
        splitted.push_back(str.substr(previous, current - previous));
        previous = current + 1;
        current = str.find(delim, previous);
    }
    splitted.push_back(str.substr(previous, current - previous));
}

bool Config::isNumber(const std::string &str)
{
    return !str.empty() && all_of(str.begin(), str.end(), ::isdigit);
}

bool Config::parseCommandLine(int argc, char **argv)
{
    // Default parameters for options are listed in header file.
    options.add_options()
            ("h,help", "Prints out usage information")
            ("g,partial", "Whether to use partial evaluations. Default: no")
            ("w,analyzeFOS", "Whether to write FOS statistics to file. Default: no")
            ("e,writeElitists", "Whether to write the genotype of the elitist solution to a file each time it is updated. Default: no")
            ("s,saveEvals", "Whether to cache evaluations. Default: no")
            ("X,parallelFOSOrder", "Whether to order the linkage model to maximize parallelizability. Default: no")
            ("Y,fixFOSOrderForPopulation", "Whether to use the same FOS order for each individual in the population. Default: no")
            //("p, useParallelGOM", "Whether to apply GOM in parallel to multiple independent linkage sets. Default: no")
            ("P,problem", "Index of optimization problem to be solved (maximization). Default: 0 (oneMax)", cxxopts::value<int>())
            ("L,numVariables", "Number of variables. Default: 10", cxxopts::value<int>())
            ("F,FOS", "FOS type, 0 - Linkage Tree, 1 - Filtered Linkage Tree. Default: Linkage Tree", cxxopts::value<int>())
            ("m,maximumFOSSetSize", "The maximum size of a linkage set. Default: no limit", cxxopts::value<int>())
            ("S,seed", "Random seed. Default: Random timestamp", cxxopts::value<long>())
            ("I,instance", "Path to problem instance.", cxxopts::value<std::string>())
            ("V,vtr", "Value to reach. Default: inf", cxxopts::value<double>())
            ("T,time", "Time limit for run in seconds. Default: no limit", cxxopts::value<double>())
            ("O,folder", "Folder where GOMEA runs. Default: \"test\"", cxxopts::value<std::string>())
            ("n,basePopulationSize", "Population size of initial population in IMS. Default: 2", cxxopts::value<int>())
            ("Z,similarityMeasure", "Measure to build Linkage Tree, 0 - Mutual Information, 1 - Normalized Mutual Information, 2 - Problem specific. Default: Mutual Information", cxxopts::value<int>())
            ("f,useForcedImprovements", "Whether to enable the Forced Improvement (FI) procedure (1: yes, 0: no). Default: 1", cxxopts::value<int>());

    auto result = options.parse(argc, argv);
    for( auto param : result.arguments())
    {
        std::string str_param = param.as<std::string>();
        char c = param.as<char>();
        switch (c)
        {
        case 'g':
            usePartialEvaluations = 1;
            break;
        case 'X':
            useParallelFOSOrder = 1;
            break;
        case 'Y':
            fixFOSOrderForPopulation = 1;
            break;
        case 'w':
            AnalyzeFOS = 1;
            break;
        case 'e':
            writeElitists = 1;
            break;
        case 's':
            saveEvaluations = 1;
            break;
        case 'h':
            printHelp = 1;
            break;
        case 'f':
            useForcedImprovements = result[str_param].as<int>();
            break;
        case 'n':
            basePopulationSize = result[str_param].as<int>();
            break;
        case 'P':
        {
            const std::string optarg_str = result[str_param].as<std::string>();
            if (isNumber(optarg_str))
                problemIndex = result[str_param].as<int>();
            else
            {
                vec_t<std::string> tmp;
                splitString(optarg_str, tmp, '_');
                std::cout << tmp[0] << " " << tmp[1] << " " << tmp[2] << std::endl;

                problemIndex = atoi(tmp[0].c_str());
                k = atoi(tmp[1].c_str());
                s = atoi(tmp[2].c_str());
                std::cout << problemIndex << std::endl;
            }
        }
        break;
        case 'F':
            FOSIndex = result[str_param].as<int>();
            break;
        case 'm':
            maximumFOSSetSize = result[str_param].as<int>();
            break;
        case 'L':
            numberOfVariables = result[str_param].as<int>();
            break;
        case 'O':
            folder = result[str_param].as<std::string>();
            break;
        case 'T':
            maximumNumberOfSeconds = result[str_param].as<double>();
            break;
        case 'V':
            vtr = result[str_param].as<double>();
            break;
        case 'S':
        {
            fix_seed = true;
            randomSeed = result[str_param].as<long long>();
        }
        break;
        case 'I':
            problemInstancePath = result[str_param].as<std::string>();
            break;
        case 'Z':
            linkage_config->lt_similarity_measure = result[str_param].as<int>();
            break;
        }
    }

    if (printHelp)
    {
        printUsage();
        exit(0);
    }

    // TODO
    // fitness = fitness::getFitnessClassDiscrete(problemIndex,numberOfVariables);
    assert(0);

    return 1;
}

void Config::printUsage()
{
  std::cout << options.help() << std::endl;
}

void Config::printOverview()
{
  std::cout << "### Settings ######################################\n";
  std::cout << "#\n";
  std::cout << "# Use partial evaluations : " << (usePartialEvaluations ? "enabled" : "disabled")  << std::endl;
  std::cout << "# Write FOS to file : " << (AnalyzeFOS ? "enabled" : "disabled") << std::endl;
  std::cout << "# Save all evaluations : " << (saveEvaluations ? "enabled" : "disabled") << std::endl;
  std::cout << "# similarity measure : " << (linkage_config->lt_similarity_measure ? "normalized MI" : "MI") << std::endl;
  
  std::cout << "#\n";
  std::cout << "###################################################\n";
  std::cout << "#\n";
  std::cout << "# Problem                      = " << fitness->name << std::endl;
  std::cout << "# Problem Instance Filename    = " << problemInstancePath << std::endl;
  std::cout << "# FOS                          = " << FOSName << std::endl;
  std::cout << "# Number of variables          = " << numberOfVariables << std::endl;
  std::cout << "# Time Limit (seconds)         = " << maximumNumberOfSeconds << std::endl;
  std::cout << "# VTR                          = " << ((vtr < 1e+308) ? std::to_string(vtr) : "not set") << std::endl;
  std::cout << "# Random seed                  = " << randomSeed << std::endl;
  std::cout << "# Folder                       = " << folder << std::endl;
  std::cout << "#\n";
  std::cout << "### Settings ######################################\n";
}

void Config::checkOptions()
{
}

}}