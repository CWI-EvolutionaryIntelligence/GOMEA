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
#include "fitness/cpp_fitness_custom.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{
namespace fitness{

cpp_customFitnessFunction_t::cpp_customFitnessFunction_t( int number_of_parameters, int number_of_subfunctions, double vtr )
{
	this->name = "Custom fitness function (C++)";
	this->number_of_parameters = number_of_parameters;
	this->number_of_subfunctions = number_of_subfunctions;
	this->vtr = vtr;
	initializeFitnessFunction();
}
		
void cpp_customFitnessFunction_t::evaluationFunction( solution_t<double> *solution )
{
	double result = 0.0;
	for( int i = 0; i < number_of_subfunctions; i++ )
	{
		// TODO
		//double fsub = subfunction( i, solution->variables );
		//solution->setPartialObjectiveValue(i,fsub);
		//result += fsub;
	}

	solution->setObjectiveValue(result);
	solution->setConstraintValue(0);
	full_number_of_evaluations++;
	number_of_evaluations++;
}

void cpp_customFitnessFunction_t::partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution )
{
	std::set<int> touched_subfunctions;
	for( int ind : solution->touched_variables )
	{
		touched_subfunctions.insert(subfunction_dependency_graph[ind].begin(), subfunction_dependency_graph[ind].end());
	}

	//vec_t<double> partial_backup = parent->createPartialBackup( solution->touched_indices );
	
	double objective_value_delta = 0.0;
	for( int subfunction_index : touched_subfunctions ) 
	{
		// TODO: backup of partial objective values
		objective_value_delta -= parent->getPartialObjectiveValue(subfunction_index);
		
		//double subf_result = subfunction( subfunction_index, parent->variables );
		//parent->setPartialObjectiveValue( subfunction_index, subf_result );
		//objective_value_delta += subf_result;
	}

	//parent->insertPartialBackup(partial_backup);
	
	solution->setObjectiveValue(parent->getObjectiveValue() + objective_value_delta);
	solution->setConstraintValue(parent->getConstraintValue());
	full_number_of_evaluations++;
	number_of_evaluations += solution->getNumberOfTouchedVariables() / (double) number_of_subfunctions;
}

double cpp_customFitnessFunction_t::subfunction( int subfunction_index, genotype_t<double> &variables )
{
	return( variables[subfunction_index] * variables[subfunction_index] );
}

double cpp_customFitnessFunction_t::getLowerRangeBound( int dimension )
{
	return( -1e308 );
}
		
double cpp_customFitnessFunction_t::getUpperRangeBound( int dimension )
{
	return( 1e308 );
}

}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
