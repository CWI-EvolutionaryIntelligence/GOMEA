#pragma once

#include <vector>
#include <unordered_map>
using namespace std;

#include "utils.hpp"
#include "time.hpp"

class sharedInformation
{
	public:
		double numberOfEvaluations;
		long long startTimeMilliseconds;
		double elitistSolutionHittingTimeMilliseconds,
			   elitistSolutionHittingTimeEvaluations;

		solutionsArchive *evaluatedSolutions;
		bool firstEvaluationEver;
		double elitistFitness;
		Individual elitist;

    sharedInformation(int maxArchiveSize)
    {
        numberOfEvaluations = 0;
        startTimeMilliseconds = getTimestamp();
        firstEvaluationEver = true;
        evaluatedSolutions = new solutionsArchive(maxArchiveSize);
    }

    ~sharedInformation()
    {
        delete evaluatedSolutions;
    }
};
