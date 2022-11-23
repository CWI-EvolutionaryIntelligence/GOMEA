#include "gomea/src/discrete/utils.hpp"

namespace gomea{
namespace discrete{

void prepareFolder(string &folder)
{
    if (!filesystem::exists(folder))
    {
		filesystem::create_directories(folder);
    }
	filesystem::create_directories(folder + "/fos");
	filesystem::create_directories(folder + "/output");
}

void initStatisticsFile(string &folder)
{
    ofstream outFile(folder + "/statistics.txt", ofstream::out);
    if (outFile.fail())
    {
        cerr << "Problems with opening file " << folder + "/elitists.txt!\n";
        exit(0);
    }
    outFile << "#Evaluations " << "Time,sec. " << "Fitness " << endl;
    outFile.close();
}

void initElitistFile(string &folder)
{
    ofstream outFile(folder + "/elitists.txt", ofstream::out);
    if (outFile.fail())
    {
        cerr << "Problems with opening file " << folder + "/elitists.txt!\n";
        exit(0);
    }
    outFile << "#Evaluations " << "Time,sec. " << "Fitness " << "Solution" << endl;
    outFile.close();
}

void writeStatisticsToFile(string &folder, long long numberOfEvaluations, long long time, solution_t<char> *solution)
{
    ofstream outFile(folder + "/statistics.txt", ofstream::app);
    if (outFile.fail())
    {
        cerr << "Problems with opening file " << folder + "/statistics.txt!\n";
        exit(0);
    }

    outFile << (int)numberOfEvaluations << " " << fixed << setprecision(3) << time/1000.0 << " " <<  setprecision(6) << solution->getObjectiveValue();
    outFile << endl;

    outFile.close();
}

void writeElitistSolutionToFile(string &folder, long long numberOfEvaluations, long long time, solution_t<char> *solution)
{
    ofstream outFile(folder + "/elitists.txt", ofstream::app);
    if (outFile.fail())
    {
        cerr << "Problems with opening file " << folder + "/elitists.txt!\n";
        exit(0);
    }

    outFile << (int)numberOfEvaluations << " " << fixed << setprecision(3) << time/1000.0 << " " <<  setprecision(6) << solution->getObjectiveValue() << " ";
    for (int i = 0; i < solution->getNumberOfVariables(); ++i)
        outFile << +solution->variables[i];
    outFile << endl;

    outFile.close();
}

void solutionsArchive::checkAlreadyEvaluated(vector<char> &genotype, archiveRecord *result)
{
    result->isFound = false;

    unordered_map<vector<char>, double, hashVector >::iterator it = archive.find(genotype);
    if (it != archive.end())
    {
        result->isFound = true;
        result->value = it->second;
    }
}

void solutionsArchive::insertSolution(vector<char> &genotype, double fitness)
{
    // #if DEBUG
    //  cout << "Inserting solution ";
    //  for (size_t i = 0; i < solution.size(); ++i)
    //      cout << solution[i];
    // #endif
    if (archive.size() >= maxArchiveSize)
        return;
    archive.insert(pair<vector<char>, double> (genotype, fitness));
}

std::chrono::high_resolution_clock::time_point start_time;
double getElapsedTime()
{
	auto end = std::chrono::high_resolution_clock::now();
	auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - start_time);
	double seconds = diff.count() / 1000.0;
	return( seconds );
}

void startTimer()
{
	start_time = std::chrono::high_resolution_clock::now();
	getElapsedTime();
}

}}