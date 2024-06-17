#pragma once

#include <vector>

#include "gomea/src/discrete/utils.hpp"
#include "gomea/src/utils/time.hpp"

namespace gomea{
namespace discrete{

class sharedInformation
{
	public:
		double elitistSolutionHittingTimeMilliseconds,
			   elitistSolutionHittingTimeEvaluations;

		solutionsArchive *evaluatedSolutions;
		bool firstEvaluationEver;
		double elitistFitness;
		double elitistConstraintValue;
		solution_t<char> elitist = solution_t<char>(1,2);

    sharedInformation(int maxArchiveSize)
    {
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