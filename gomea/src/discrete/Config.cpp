#include "gomea/src/discrete/Config.hpp"

namespace gomea{
namespace discrete{

Config::Config(){};

Config::~Config(){};

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
            ("P,problems", "Prints out all installed problems")
            ("problem_index", "Index of optimization problem to be solved (maximization).", cxxopts::value<int>())
            ("L,numVariables", "Number of variables.", cxxopts::value<int>())
            ("g,partial", "Whether to use partial evaluations. Default: yes")
            ("S,random_seed", "Random seed. Default: Random timestamp", cxxopts::value<long>())
            ("instance", "Path to problem instance.", cxxopts::value<std::string>())
            //("V,vtr", "Value to reach. Default: inf", cxxopts::value<double>())
            ("VTRfile", "Path to VTR file.", cxxopts::value<std::string>())
            ("O,folder", "Folder where output is saved to. Default: \"output_discrete_gomea\"", cxxopts::value<std::string>())
            ("F,max_number_of_evaluations", "Maximum number of function evaluations. Default: inf", cxxopts::value<int>())
            ("T,max_number_of_seconds", "Time limit for run in seconds. Default: inf", cxxopts::value<double>())
            ("n,base_population_size", "Population size of initial population in IMS. Default: 2", cxxopts::value<int>())
            ("max_number_of_populations", "Maximum number of populations in IMS. Default: 25", cxxopts::value<int>())
            ("max_number_of_generations", "Maximum number of generations per population. Default: inf", cxxopts::value<int>())
            ("IMS_subgeneration_factor", "Factor in which larger populations are started in IMS. Default: 4", cxxopts::value<int>())
            // Linkage config
            ("FOS_index", "Index of the FOS to be used.", cxxopts::value<int>())
            ("mpm_block_size", "Block size of marginal product FOS.", cxxopts::value<int>())
            ("lt_similarity_measure", "Similarity measure for linkage tree FOS.", cxxopts::value<int>())
            ("filtered_lt", "Enable for filtered linkage tree.")
            ("static_lt", "Enable for static linkage tree.")
            ("cond_max_clique_size", "Maximum clique size for conditional linkage models.", cxxopts::value<int>())
            ("cond_include_cliques_as_fos_elements", "Include cliques as FOS elements in conditional linkage models.")
            ("cond_include_full_fos_element", "Include full FOS element in conditional linkage models.")
            ("filename", "Filename for FOS matrix.", cxxopts::value<std::string>());
    options.parse_positional({"problem_index", "numVariables", "FOS_index"});
    options.positional_help("[problem_index] [numVariables] [FOS_index]").show_positional_help();

    if (argc == 1 )
    {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
        std::cout << options.help() << std::endl;
        exit(0);
    }
    if (result.count("problems"))
    {
        printf("TODO - print all installed problems\n");
        //printAllInstalledProblems();
        exit(0);
    }

    int problem_index = result["problem_index"].as<int>();
    int num_variables = result["numVariables"].as<int>();
    int FOSIndex = result["FOS_index"].as<int>();
    std::string vtrFilePath = "";
    if(result.count("partial"))
    {
        usePartialEvaluations = result["partial"].as<bool>() ? 1 : 0;
    }
    if(result.count("random_seed"))
    {
        fix_seed = true;
        randomSeed = result["random_seed"].as<long>();
    }
    if(result.count("instance"))
    {
        problemInstancePath = result["instance"].as<std::string>();
    }
    /*if(result.count("vtr"))
    {
        vtr = result["vtr"].as<double>();
    }*/
    if(result.count("VTRfile"))
    {
        vtrFilePath = result["VTRfile"].as<std::string>();
    }
    if(result.count("folder"))
    {
        folder = result["folder"].as<std::string>();
    }
    if(result.count("max_number_of_evaluations"))
    {
        maximumNumberOfEvaluations = result["max_number_of_evaluations"].as<int>();
    }
    if(result.count("max_number_of_seconds"))
    {
        maximumNumberOfSeconds = result["max_number_of_seconds"].as<double>();
    }
    if(result.count("base_population_size"))
    {
        basePopulationSize = result["base_population_size"].as<int>();
    }
    if(result.count("max_number_of_populations"))
    {
        maximumNumberOfGOMEAs = result["max_number_of_populations"].as<int>();
    }
    if(result.count("max_number_of_generations"))
    {
        maximumNumberOfGenerations = result["max_number_of_generations"].as<int>();
    }
    if(result.count("IMS_subgeneration_factor"))
    {
        IMSsubgenerationFactor = result["IMS_subgeneration_factor"].as<int>();
    }

    // Linkage model parameters
    int mpm_block_size = -1, lt_similarity_measure = 0, cond_max_clique_size = 1, lt_max_set_size = -1;
    bool filtered_lt, static_lt, cond_include_cliques_as_fos_elements, cond_include_full_fos_element;
    std::string filename;
    if(result.count("FOS_index"))
    {
        FOSIndex = result["FOS_index"].as<int>();
    }
    if(result.count("mpm_block_size"))
    {
        mpm_block_size = result["mpm_block_size"].as<int>();
    }
    if(result.count("lt_similarity_measure"))
    {
        lt_similarity_measure = result["lt_similarity_measure"].as<int>();
    }
    if(result.count("lt_max_set_size"))
    {
        lt_max_set_size = result["lt_max_set_size"].as<int>();
    }
    if(result.count("filtered_lt"))
    {
        filtered_lt = true;
    }
    if(result.count("static_lt"))
    {
        static_lt = true;
    }
    if(result.count("cond_max_clique_size"))
    {
        cond_max_clique_size = result["cond_max_clique_size"].as<int>();
    }
    if(result.count("cond_include_cliques_as_fos_elements"))
    {
        cond_include_cliques_as_fos_elements = true;
    }
    if(result.count("cond_include_full_fos_element"))
    {
        cond_include_full_fos_element = true;
    }
    if(result.count("filename"))
    {
        filename = result["filename"].as<std::string>();
    }

    switch(FOSIndex){
        case linkage::linkage_model_type::UNIVARIATE:
            linkage_config = new linkage_config_t();
            break;
        case linkage::linkage_model_type::FULL:
            throw std::invalid_argument("Full FOS is invalid for discrete optimization.");
            break;
        case linkage::linkage_model_type::MPM:
            linkage_config = new linkage_config_t(true, mpm_block_size);
            break;
        case linkage::linkage_model_type::LINKAGE_TREE:
            linkage_config = new linkage_config_t(lt_similarity_measure, filtered_lt, lt_max_set_size, static_lt);
            break;
        case linkage::linkage_model_type::CONDITIONAL:
            linkage_config = new linkage_config_t(cond_max_clique_size, cond_include_cliques_as_fos_elements, cond_include_full_fos_element);
            break;
        case linkage::linkage_model_type::FROM_FILE:
            linkage_config = new linkage_config_t(filename);
            break;
        default:
            throw std::invalid_argument("Invalid FOS index.");
    }
    
    if(!usePartialEvaluations){
        switch(problem_index){
            case 0:
                fitness = new fitness::oneMaxBBO_t(num_variables);
                break;
            case 1:
                fitness = new fitness::deceptiveTrapBBO_t(num_variables, 5);
                break;
            case 2:
                fitness = new fitness::maxCutBBO_t(problemInstancePath, vtrFilePath);
                break;
            default:
                throw std::invalid_argument("Invalid problem index.");
        }
    }
    else{
        switch(problem_index){
            case 0:
                fitness = new fitness::oneMax_t(num_variables);
                break;
            case 1:
                fitness = new fitness::deceptiveTrap_t(num_variables, 5);
                break;
            case 2:
                fitness = new fitness::maxCut_t(problemInstancePath, vtrFilePath);
                break;
            default:
                throw std::invalid_argument("Invalid problem index.");
        }
    }

    return 1;
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
  std::cout << "# Number of variables          = " << fitness->number_of_variables << std::endl;
  std::cout << "# Time Limit (seconds)         = " << maximumNumberOfSeconds << std::endl;
  std::cout << "# VTR                          = " << ((fitness->vtr < 1e+308) ? std::to_string(fitness->vtr) : "not set") << std::endl;
  std::cout << "# Random seed                  = " << randomSeed << std::endl;
  std::cout << "# Folder                       = " << folder << std::endl;
  std::cout << "#\n";
  std::cout << "### Settings ######################################\n";
}

void Config::checkOptions()
{
}

}}