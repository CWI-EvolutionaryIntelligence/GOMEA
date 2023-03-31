#pragma once

#include <vector>
using namespace std;

#include "gomea/src/discrete/utils.hpp"
#include "gomea/src/utils/time.hpp"

namespace gomea{
namespace discrete{

class sharedInformation
{
	public:
		time_t startTime;
		double elitistSolutionHittingTimeMilliseconds,
			   elitistSolutionHittingTimeEvaluations;

		solutionsArchive *evaluatedSolutions;
		bool firstEvaluationEver;
		double elitistFitness;
		double elitistConstraintValue;
		solution_t<char> elitist = solution_t<char>(1,2);

    sharedInformation(int maxArchiveSize)
    {
        startTime = utils::getTimestamp();
        firstEvaluationEver = true;
        evaluatedSolutions = new solutionsArchive(maxArchiveSize);
        gomea::utils::clearTimers();
    }

    ~sharedInformation()
    {
        delete evaluatedSolutions;
    }
};

}}