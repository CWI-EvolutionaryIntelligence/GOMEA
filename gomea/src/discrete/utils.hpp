#pragma once

#include <sys/stat.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <exception>
#include <cassert>

using namespace std;

#include "gomea/src/common/solution.hpp"
#include "gomea/src/common/gomea_defs.hpp"

namespace gomea{
namespace discrete{

typedef std::chrono::time_point<std::chrono::steady_clock> chtime;

struct archiveRecord
{
  bool isFound = false;
  double value = 0.0;
};

struct hashVector
{ 
    size_t operator()(const vector<char> &vec) const
    { 
        hash <char> hashChar; 
        size_t hash_value = 0;
        for (size_t i = 0; i < vec.size(); ++i) 
            hash_value = hash_value*31 + hashChar(vec[i]); 
        return hash_value; 
    } 
}; 

class solutionsArchive
{
    size_t maxArchiveSize;
public:
    solutionsArchive(size_t maxArchiveSize_): maxArchiveSize(maxArchiveSize_){};
    unordered_map<vector<char>, double, hashVector > archive;
    void checkAlreadyEvaluated(vector<char> &genotype, archiveRecord *result);
    void insertSolution(vector<char> &genotype, double fitness);
};

void prepareFolder(string &folder);
void initElitistFile(string &folder);
void writeStatisticsToFile(string &folder, long long numberOfEvaluations, long long time, solution_t<char> *solution);
void writeElitistSolutionToFile(string &folder, long long numberOfEvaluations, long long time, solution_t<char> *solution);

}}