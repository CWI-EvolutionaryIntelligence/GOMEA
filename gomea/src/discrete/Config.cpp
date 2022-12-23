#include "gomea/src/discrete/Config.hpp"

namespace gomea{
namespace discrete{

Config::Config(){}

/**
 * Parses the command line.
 * For options, see printUsage.
 */
void Config::splitString(const string &str, vector<string> &splitted, char delim)
{
    size_t current, previous = 0;
    current = str.find(delim);
    while (current != string::npos)
    {
        splitted.push_back(str.substr(previous, current - previous));
        previous = current + 1;
        current = str.find(delim, previous);
    }
    splitted.push_back(str.substr(previous, current - previous));
}

bool Config::isNumber(const string &str)
{
    return !str.empty() && all_of(str.begin(), str.end(), ::isdigit);
}

bool Config::parseCommandLine(int argc, char **argv)
{
  const struct option longopts[] =
  {
    {"help",        no_argument,         0, 'h'},    
    {"partial",     no_argument,         0, 'g'},
    {"analyzeFOS",  no_argument,         0, 'w'},
    {"writeElitists",no_argument,        0, 'e'},
    {"saveEvals",   no_argument,         0, 's'},  
    {"parallelFOSOrder",no_argument,     0, 'X'},  
    {"fixFOSOrderForPopulation",no_argument, 0, 'Y'},  
    {"popUpdatesDuringGOM", no_argument, 0, 'Q'},  
    {"useParallelGOM", required_argument,0, 'p'},
    {"problem",     required_argument,   0, 'P'},  
    {"L",           required_argument,   0, 'L'},  
    {"FOS",         required_argument,   0, 'F'},
    {"maximumFOSSetSize",required_argument, 0, 'm'},
    {"seed",        required_argument,   0, 'S'},
    {"instance",    required_argument,   0, 'I'},
    {"vtr",         required_argument,   0, 'V'},
    {"time",        required_argument,   0, 'T'},
    {"folder",      required_argument,   0, 'O'},
    {"basePopulationSize", required_argument,   0, 'n'}, 
    {"similarityMeasure", required_argument,   0, 'Z'}, 
    {"useForcedImprovements", required_argument,   0, 'f'}, 
    {"GPUIndex", required_argument,   0, 'G'}, 
               
    {0,             0,                   0,  0 }
  };


  int c, index;
  while ((c = getopt_long(argc, argv, "h::n::p::X::Y::Q::g::w::e::s::f::P::F::m::L::O::T::S::V::I::B::Z::G::", longopts, &index)) != -1)
  {
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
		case 'Q':
			popUpdatesDuringGOM = 1;
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
            useForcedImprovements = atoi(optarg);
            break;
        case 'n':
            basePopulationSize = atoi(optarg);
            break;
        case 'P':
            {
                const string optarg_str = string(optarg);
                if (isNumber(optarg_str))   
                    problemIndex = atoi(optarg);
                else
                {
                    vector<string> tmp;
                    splitString(optarg_str, tmp, '_');
                    cout << tmp[0] << " " << tmp[1] << " " << tmp[2] << endl;   
                
                    problemIndex = atoi(tmp[0].c_str());
                    k = atoi(tmp[1].c_str());
                    s = atoi(tmp[2].c_str());
                    cout << problemIndex << endl;
                }
            }
            break;
        case 'p':
            useParallelGOM = atoi(optarg);
            break;
        case 'F':
            FOSIndex = atoi(optarg);
            break;
        case 'm':
            maximumFOSSetSize = atoi(optarg);
            break;
        case 'L':
            numberOfVariables = atoi(optarg);
            break;
        case 'O':
            folder= string(optarg);
            break;
        case 'T':
            maximumNumberOfSeconds = atof(optarg);
            break;
        case 'V':
            vtr = atof(optarg);
            break;
        case 'S':
			{
                fix_seed = true;
				randomSeed = atoll(optarg);
			}
            break;
        case 'I':
            problemInstancePath = string(optarg);
            break;
        case 'Z':
            linkage_config->lt_similarity_measure = atoi(optarg);
            break;
        case 'G':
            GPUIndex = atoi(optarg);
            break;
        default:
            abort();
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
  cout << "  -h: Prints out this usage information.\n";
  cout << "  --partial: Whether to use partial evaluations. Default: no\n";
  cout << "  --analyzeFOS: Whether to write FOS statistics to file. Default: no\n";
  cout << "  --writeElitists: Whether to write the genotype of the elitist solution to a file each time it is updated. Default: no\n";
  cout << "  --saveEvals: Whether to cache evaluations. Default: no.\n";
  cout << endl;
  cout << "  --problem: Index of optimization problem to be solved (maximization). Default: 0 (oneMax)\n";

  cout << "  List of available problems:\n";
  cout << "    0: OneMax\n";
  cout << "    1: Concatenated Deceptive Trap, specify k and s by writing 1_k_s\n";
  cout << "    2: ADF, specify k and s by writing 2_k_s\n";
  cout << "    3: MaxCut\n";
  cout << "    4: Hierarhical If-And-Only-If\n";
  cout << "    5: Leading Ones\n";        
  cout << "    6: Hierarhical Deceptive Trap-3\n";
  cout << "    7: Concatenated Bimodal Deceptive Trap, specify k and s by writing 7_k_s\n";
  cout << endl;
  
  cout << "  --L: Number of variables. Default: 10\n";
  cout << "  --FOS: FOS type, 0 - Linkage Tree, 1 - Filtered Linkage Tree. Default: Linkage Tree\n";    
  cout << "  --similarityMeasure: Measure to build Linkage Tree, 0 - Mutual Information, 1 - Normalized Mutual Information, 2 - Problem specific. Default: Mutual Information\n";      
  cout << "  --folder: Folder where GOMEA runs. Default: \"test\"\n";
  cout << "  --instance: Path to problem instance. Default: \"\"\n";
  cout << "  --vtr: Value to reach. Default: inf\n";  
  cout << "  --time: time limit for run in seconds. Default: 1\n";
  cout << "  --seed: Random seed. Default: Random timestamp\n";  
}

void Config::printOverview()
{
  cout << "### Settings ######################################\n";
  cout << "#\n";
  cout << "# Use partial evaluations : " << (usePartialEvaluations ? "enabled" : "disabled")  << endl;
  cout << "# Write FOS to file : " << (AnalyzeFOS ? "enabled" : "disabled") << endl;
  cout << "# Save all evaluations : " << (saveEvaluations ? "enabled" : "disabled") << endl;
  cout << "# similarity measure : " << (linkage_config->lt_similarity_measure ? "normalized MI" : "MI") << endl;
  
  cout << "#\n";
  cout << "###################################################\n";
  cout << "#\n";
  cout << "# Problem                      = " << fitness->name << endl;
  cout << "# Problem Instance Filename    = " << problemInstancePath << endl;
  cout << "# FOS                          = " << FOSName << endl;
  cout << "# Number of variables          = " << numberOfVariables << endl;
  cout << "# Time Limit (seconds)         = " << maximumNumberOfSeconds << endl;
  cout << "# VTR                          = " << ((vtr < 1e+308) ? to_string(vtr) : "not set") << endl;
  cout << "# Random seed                  = " << randomSeed << endl;
  cout << "# Folder                       = " << folder << endl;
  cout << "#\n";
  cout << "### Settings ######################################\n";
}

void Config::checkOptions()
{
}

}}