#pragma once

#include <vector>
#include <unordered_map>
using namespace std;

#include "gomea/src/discrete/utils.hpp"
#include "gomea/src/utils/time.hpp"

namespace gomea{
namespace discrete{

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
		Individual elitist = Individual(1,2);

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

}}