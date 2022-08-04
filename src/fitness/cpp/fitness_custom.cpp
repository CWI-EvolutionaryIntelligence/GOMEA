/**
 *
 * RV-GOMEA
 *
 * If you use this software for any purpose, please cite the most recent publication:
 * A. Bouter, C. Witteveen, T. Alderliesten, P.A.N. Bosman. 2017.
 * Exploiting Linkage Information in Real-Valued Optimization with the Real-Valued
 * Gene-pool Optimal Mixing Evolutionary Algorithm. In Proceedings of the Genetic 
 * and Evolutionary Computation Conference (GECCO 2017).
 * DOI: 10.1145/3071178.3071272
 *
 * Copyright (c) 1998-2017 Peter A.N. Bosman
 *
 * The software in this file is the proprietary information of
 * Peter A.N. Bosman.
 *
 * IN NO EVENT WILL THE AUTHOR OF THIS SOFTWARE BE LIABLE TO YOU FOR ANY
 * DAMAGES, INCLUDING BUT NOT LIMITED TO LOST PROFITS, LOST SAVINGS, OR OTHER
 * INCIDENTIAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR THE INABILITY
 * TO USE SUCH PROGRAM, EVEN IF THE AUTHOR HAS BEEN ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGES, OR FOR ANY CLAIM BY ANY OTHER PARTY. THE AUTHOR MAKES NO
 * REPRESENTATIONS OR WARRANTIES ABOUT THE SUITABILITY OF THE SOFTWARE, EITHER
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR NON-INFRINGEMENT. THE
 * AUTHOR SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY ANYONE AS A RESULT OF
 * USING, MODIFYING OR DISTRIBUTING THIS SOFTWARE OR ITS DERIVATIVES.
 *
 * The software in this file is the result of (ongoing) scientific research.
 * The following people have been actively involved in this research over
 * the years:
 * - Peter A.N. Bosman
 * - Dirk Thierens
 * - Jörn Grahl
 * - Anton Bouter
 * 
 */

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "utils/embed.hpp"
#include "fitness/fitness_basic.hpp"
#include "fitness/fitness_custom.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{
namespace fitness{

customFitnessFunction_t::customFitnessFunction_t( int number_of_parameters, double vtr )
{
	this->name = "Sphere function";
	this->number_of_parameters = number_of_parameters;
	this->number_of_subfunctions = number_of_parameters;
	this->vtr = vtr;
	initializeFitnessFunction();
}
		
void customFitnessFunction_t::evaluationFunction( solution_t<double> *solution )
{
	double result = 0.0;
	for( int i = 0; i < number_of_subfunctions; i++ )
		result += subfunction( solution->variables[i] );

	solution->setObjectiveValue(result);
	solution->setConstraintValue(0);
	full_number_of_evaluations++;
	number_of_evaluations++;
}

void customFitnessFunction_t::partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution )
{
	double result = 0.0;
	for( int i = 0; i < solution->getNumberOfTouchedVariables(); i++ )
	{
		int ind = solution->touched_indices[i];
		result += subfunction( solution->touched_variables[i] );
		result -= subfunction( parent->variables[ind] );
	}
	
	solution->setObjectiveValue(parent->getObjectiveValue() + result);
	solution->setConstraintValue(parent->getConstraintValue());
	full_number_of_evaluations++;
	number_of_evaluations += solution->getNumberOfTouchedVariables() / (double) number_of_subfunctions;
}

double customFitnessFunction_t::subfunction( double x )
{
	return( x * x );
}

double customFitnessFunction_t::getLowerRangeBound( int dimension )
{
	return( -1e308 );
}
		
double customFitnessFunction_t::getUpperRangeBound( int dimension )
{
	return( 1e308 );
}

double customFitnessFunction_t::subfunctionEmbedded( int subfunction_index, double *variables )
{
	// TODO
	/*if (fitness_embedded() < 0) {
        PyErr_Print();
        fprintf(stderr, "Error in Python code, exception was printed.\n");
		gomea::utils::freePythonEmbedding();
		exit(1);
    }*/
	return 0;
}
		
}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
