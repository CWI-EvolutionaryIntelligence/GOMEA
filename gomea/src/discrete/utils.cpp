#include "gomea/src/discrete/utils.hpp"

namespace gomea{
namespace discrete{

void prepareFolder(std::string &folder)
{
    if (!std::filesystem::exists(folder))
    {
		std::filesystem::create_directories(folder);
    }
	std::filesystem::create_directories(folder + "/fos");
	std::filesystem::create_directories(folder + "/output");
}

void initStatisticsFile(std::string &folder)
{
    std::ofstream outFile(folder + "/statistics.txt", std::ofstream::out);
    if (outFile.fail())
    {
        std::cerr << "Problems with opening file " << folder + "/elitists.txt!\n";
        exit(0);
    }
    outFile << "#Evaluations " << "Time,sec. " << "Fitness " << std::endl;
    outFile.close();
}

void initElitistFile(std::string &folder)
{
    std::ofstream outFile(folder + "/elitists.txt", std::ofstream::out);
    if (outFile.fail())
    {
        std::cerr << "Problems with opening file " << folder + "/elitists.txt!\n";
        exit(0);
    }
    outFile << "#Evaluations " << "Time,sec. " << "Fitness " << "Solution" << std::endl;
    outFile.close();
}

void writeStatisticsToFile(std::string &folder, long long numberOfEvaluations, long long time, solution_t<char> *solution)
{
    std::ofstream outFile(folder + "/statistics.txt", std::ofstream::app);
    if (outFile.fail())
    {
        std::cerr << "Problems with opening file " << folder + "/statistics.txt!\n";
        exit(0);
    }

    outFile << (int)numberOfEvaluations << " " << std::fixed << std::setprecision(3) << time/1000.0 << " " <<  std::setprecision(6) << solution->getObjectiveValue();
    outFile << std::endl;

    outFile.close();
}

void writeElitistSolutionToFile(std::string &folder, long long numberOfEvaluations, long long time, solution_t<char> *solution)
{
    std::ofstream outFile(folder + "/elitists.txt", std::ofstream::app);
    if (outFile.fail())
    {
        std::cerr << "Problems with opening file " << folder + "/elitists.txt!\n";
        exit(0);
    }

    outFile << (int)numberOfEvaluations << " " << std::fixed << std::setprecision(3) << time/1000.0 << " " <<  std::setprecision(6) << solution->getObjectiveValue() << " ";
    for (int i = 0; i < solution->getNumberOfVariables(); ++i)
        outFile << +solution->variables[i];
    outFile << std::endl;

    outFile.close();
}

void solutionsArchive::checkAlreadyEvaluated(vec_t<char> &genotype, archiveRecord *result)
{
    result->isFound = false;

    std::unordered_map<vec_t<char>, double, hashVector >::iterator it = archive.find(genotype);
    if (it != archive.end())
    {
        result->isFound = true;
        result->value = it->second;
    }
}

void solutionsArchive::insertSolution(vec_t<char> &genotype, double fitness)
{
    // #if DEBUG
    //  std::cout << "Inserting solution ";
    //  for (size_t i = 0; i < solution.size(); ++i)
    //      std::cout << solution[i];
    // #endif
    if (archive.size() >= maxArchiveSize)
        return;
    archive.insert(std::pair<vec_t<char>, double> (genotype, fitness));
}


}}