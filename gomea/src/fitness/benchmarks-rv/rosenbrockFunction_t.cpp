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
#include "gomea/src/fitness/benchmarks-rv.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{
namespace fitness{

using namespace gomea;

rosenbrockFunction_t::rosenbrockFunction_t( int number_of_variables, double vtr ) : fitness_t(number_of_variables,vtr)
{
	this->name = "Rosenbrock function";
	if( !black_box_optimization )
		initializeVariableInteractionGraph();
}

int rosenbrockFunction_t::getNumberOfSubfunctions()
{
	return number_of_variables-1;
}

void rosenbrockFunction_t::evaluationFunction( solution_t<double> *solution )
{
	double result = 0.0;
	for( int i = 0; i < number_of_variables-1; i++ )
		result += subfunction( solution->variables[i], solution->variables[i+1] );

	solution->setObjectiveValue(result);
	solution->setConstraintValue(0);
	full_number_of_evaluations++;
	number_of_evaluations++;
}

void rosenbrockFunction_t::univariatePartialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution )
{
	assert( solution->getNumberOfTouchedVariables() == 1 );

	int num_subfunctions_evaluated = 0;	
	double result = 0.0;
	int ind = solution->touched_indices[0];
	if( ind > 0 )
	{
		result += subfunction( parent->variables[ind-1], solution->touched_variables[0] );
		result -= subfunction( parent->variables[ind-1], parent->variables[ind] );
		num_subfunctions_evaluated++;
	}
	if( ind < getNumberOfSubfunctions() )
	{
		result += subfunction( solution->touched_variables[0], parent->variables[ind+1] );
		result -= subfunction( parent->variables[ind], parent->variables[ind+1] );
		num_subfunctions_evaluated++;
	}

	solution->setObjectiveValue(parent->getObjectiveValue() + result);
	solution->setConstraintValue(parent->getConstraintValue());
	full_number_of_evaluations++;
	number_of_evaluations += num_subfunctions_evaluated / (double) getNumberOfSubfunctions();
}


void rosenbrockFunction_t::partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution )
{
	if( solution->getNumberOfTouchedVariables() == 1 )
	{
		univariatePartialEvaluationFunction( parent, solution );
	}
	else
	{
		int num_subfunctions_evaluated = 0;
		double result = 0.0;
		vec_t<int> order = gomea::utils::getSortedOrder( solution->touched_indices );
		for( int i = 0; i < solution->getNumberOfTouchedVariables(); i++ )
		{
			int ind = solution->touched_indices[order[i]];
			if( ind > 0 )
			{
				double y = parent->variables[ind-1];
				if( i > 0 && solution->touched_indices[i-1] == ind-1 )
					y = solution->touched_variables[i-1];

				result += subfunction( y, solution->touched_variables[i] );
				result -= subfunction( parent->variables[ind-1], parent->variables[ind] );
				num_subfunctions_evaluated++;
			}
			if( ind < getNumberOfSubfunctions() )
			{
				double y = parent->variables[ind+1];
				if( !(i < solution->getNumberOfTouchedVariables()-1 && solution->touched_indices[i+1] == ind+1) )
				{
					result += subfunction( solution->touched_variables[i], y );
					result -= subfunction( parent->variables[ind], parent->variables[ind+1] );
					num_subfunctions_evaluated++;
				}
			}
		}

		solution->setObjectiveValue(parent->getObjectiveValue() + result);
		solution->setConstraintValue(parent->getConstraintValue());
		full_number_of_evaluations++;
		number_of_evaluations += num_subfunctions_evaluated / (double) getNumberOfSubfunctions();
	}	
}

double rosenbrockFunction_t::subfunction( double x, double y )
{
	return( 100*(y-x*x)*(y-x*x) + (1.0-x)*(1.0-x) );
}

void rosenbrockFunction_t::initializeVariableInteractionGraph()
{
	for( int i = 0; i < number_of_variables; i++ )
	{
		std::set<int> dependent_vars; 
		if( i > 0 )
			dependent_vars.insert(i-1);
		if( i < number_of_variables-1 )
			dependent_vars.insert(i+1);
		variable_interaction_graph[i] = dependent_vars;
	}
}

}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
