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
#include "fitness/benchmarks-rv.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{
namespace fitness{

using namespace gomea;

sphereFunction_t::sphereFunction_t( int number_of_parameters, double vtr )
{
	this->name = "Sphere function";
	this->number_of_parameters = number_of_parameters;
	this->number_of_subfunctions = number_of_parameters;
	this->vtr = vtr;
	initializeFitnessFunction();
}
		
void sphereFunction_t::evaluationFunction( solution_t<double> *solution )
{
	double result = 0.0;
	for( int i = 0; i < number_of_subfunctions; i++ )
		result += subfunction( solution->variables[i] );

	/*if( number_of_evaluations == 0 )
	{
		fitness_embedded_pub();
		printf("TEST - SPHERE\n");
	}*/

	solution->setObjectiveValue(result);
	solution->setConstraintValue(0);
	full_number_of_evaluations++;
	number_of_evaluations++;
}

void sphereFunction_t::partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution )
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

double sphereFunction_t::subfunction( double x )
{
	return( x * x );
}

double sphereFunction_t::getLowerRangeBound( int dimension )
{
	return( -1e308 );
}
		
double sphereFunction_t::getUpperRangeBound( int dimension )
{
	return( 1e308 );
}

rosenbrockFunction_t::rosenbrockFunction_t( int number_of_parameters, double vtr )
{
	this->name = "Rosenbrock function";
	this->number_of_parameters = number_of_parameters;
	this->vtr = vtr;
	this->number_of_subfunctions = number_of_parameters-1;
	initializeFitnessFunction();
	if( !black_box_optimization )
		initializeVariableInteractionGraph();
}

void rosenbrockFunction_t::evaluationFunction( solution_t<double> *solution )
{
	double result = 0.0;
	for( int i = 0; i < number_of_parameters-1; i++ )
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
	if( ind < number_of_subfunctions )
	{
		result += subfunction( solution->touched_variables[0], parent->variables[ind+1] );
		result -= subfunction( parent->variables[ind], parent->variables[ind+1] );
		num_subfunctions_evaluated++;
	}

	solution->setObjectiveValue(parent->getObjectiveValue() + result);
	solution->setConstraintValue(parent->getConstraintValue());
	full_number_of_evaluations++;
	number_of_evaluations += num_subfunctions_evaluated / (double) number_of_subfunctions;
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
			if( ind < number_of_subfunctions )
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
		number_of_evaluations += num_subfunctions_evaluated / (double) number_of_subfunctions;
	}	
}

double rosenbrockFunction_t::subfunction( double x, double y )
{
	return( 100*(y-x*x)*(y-x*x) + (1.0-x)*(1.0-x) );
}

double rosenbrockFunction_t::getLowerRangeBound( int dimension )
{
	return( -1e308 );
}
		
double rosenbrockFunction_t::getUpperRangeBound( int dimension )
{
	return( 1e308 );
}

void rosenbrockFunction_t::initializeVariableInteractionGraph()
{
	for( int i = 0; i < number_of_parameters; i++ )
	{
		std::set<int> dependent_vars; 
		if( i > 0 )
			dependent_vars.insert(i-1);
		if( i < number_of_parameters-1 )
			dependent_vars.insert(i+1);
		variable_interaction_graph[i] = dependent_vars;
	}
	/*for( auto p : variable_interaction_graph )
	{
		printf("[%d] ",p.first);
		for( int x : p.second )
			printf("%d ",x);
	}
	printf("\n");*/
}


sorebFunction_t::sorebFunction_t( int number_of_parameters, double vtr, double conditioning_number, double rotation_angle, int block_size, int overlap_size )
{
	this->name = "Sum of Rotated Ellipsoid Blocks function";
	this->number_of_parameters = number_of_parameters;
	this->vtr = vtr;
	this->rotation_angle = rotation_angle;
	this->conditioning_number = conditioning_number;
	this->rotation_block_size = block_size;
	this->overlap_size = overlap_size;
	this->number_of_subfunctions = 1+(number_of_parameters-block_size)/(block_size-overlap_size);

	assert( number_of_parameters >= block_size );
	assert( overlap_size < block_size );
	assert( (number_of_parameters-block_size) % (block_size-overlap_size) == 0 ); // total number of parameters matches an integer number of blocks

	initializeFitnessFunction();
	rotation_matrix = initializeObjectiveRotationMatrix( rotation_angle, rotation_block_size );
	if( !black_box_optimization )
		initializeVariableInteractionGraph();
}

int sorebFunction_t::getIndexOfFirstBlock( int var )
{
	int block_index = (var-overlap_size) / (rotation_block_size-overlap_size); // rounded down
	block_index = fmax( block_index, 0 );
	block_index = fmin( block_index, number_of_subfunctions-1 );
	return( block_index );
}

int sorebFunction_t::getStartingIndexOfBlock( int block_index )
{
	return( block_index * (rotation_block_size-overlap_size) );
}

void sorebFunction_t::evaluationFunction( solution_t<double> *solution )
{
	double result = 0.0;
	for( int i = 0; i < number_of_subfunctions; i++ )
		result += subfunction( &solution->variables[getStartingIndexOfBlock(i)], rotation_block_size );

	solution->setObjectiveValue(result);
	solution->setConstraintValue(0);
	full_number_of_evaluations++;
	number_of_evaluations++;
}

void sorebFunction_t::partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution )
{
	double result = 0.0;
	int last_evaluated_block = -1;
	double variables_copy[rotation_block_size];
	int num_subfunctions_evaluated = 0;

	/*printf("ORIGINAL VARS:\n");
	for( int i = 0; i < number_of_parameters; i++ )
		printf("%.3lf ",parent->variables[i]);
	printf("\n");
	evaluationFunction(parent);
	printf("fbefore = %10.3e\n",parent->getObjectiveValue());

	printf("MODIFIED VARS:\n");
	for( int i = 0; i < solution->getNumberOfTouchedVariables(); i++ )
		printf("%d ",solution->touched_indices[i]);
	printf("\n");
	for( int i = 0; i < solution->getNumberOfTouchedVariables(); i++ )
		printf("%.3lf ",solution->touched_variables[i]);
	printf("\n");*/
	for( int i = 0; i < solution->getNumberOfTouchedVariables(); i++ )
	{
		int ind = solution->touched_indices[i];
		int block_ind = getIndexOfFirstBlock(ind);
		if( i > 0 )
			assert( ind > solution->touched_indices[i-1] ); // this partial evaluation requires indices to be sorted
	
		while( block_ind < number_of_subfunctions && getStartingIndexOfBlock(block_ind) <= ind ) // while the touched variable is within the block to be evaluated
		{
			if( block_ind > last_evaluated_block )
			{
				int block_start = getStartingIndexOfBlock( block_ind );
				for( int j = 0; j < rotation_block_size; j++ )
					variables_copy[j] = parent->variables[block_start+j];
				result -= subfunction( variables_copy, rotation_block_size );

				int j = 0;
				while( i+j < solution->getNumberOfTouchedVariables() && solution->touched_indices[i+j]-block_start < rotation_block_size ) // find other touched variables in this block and put them in local array
				{
					int cur_ind = solution->touched_indices[i+j];
					variables_copy[cur_ind-block_start] = solution->touched_variables[i+j];
					j++;
				}
				result += subfunction( variables_copy, rotation_block_size );
				num_subfunctions_evaluated++;

				last_evaluated_block = block_ind;
			}

			block_ind++;
		}
	}
	
	solution->setObjectiveValue(parent->getObjectiveValue() + result);
	solution->setConstraintValue(parent->getConstraintValue());
	full_number_of_evaluations++;
	number_of_evaluations += num_subfunctions_evaluated / (double) number_of_subfunctions;
}

double sorebFunction_t::subfunction( double *vars, int num_vars )
{
	double *rotated_vars = vars;
   	if( rotation_angle != 0.0 )
		rotated_vars = rotateVariables( vars, num_vars, rotation_matrix );
	double result = 0.0;
	
	/*for( int i = 0; i < num_vars; i++ )
		rotated_vars[i] += 1.0;
	for( int i = 0; i < num_vars-1; i++ )
		result += 100*(rotated_vars[i+1]-rotated_vars[i]*rotated_vars[i])*(rotated_vars[i+1]-rotated_vars[i]*rotated_vars[i]) + (1.0-rotated_vars[i])*(1.0-rotated_vars[i]);*/

	//result += exp(rotated_vars[0]*rotated_vars[0] + rotated_vars[1]*rotated_vars[1] + rotated_vars[0]*rotated_vars[1]       ) - 1.0;

    	for(int i = 0; i < num_vars; i++ )
        	result += pow( 10.0, conditioning_number*(((double) (i))/((double) (rotation_block_size-1))) )*rotated_vars[i]*rotated_vars[i];
	
	if( rotation_angle != 0.0 )
		free( rotated_vars );
	return( result ); 
}

double sorebFunction_t::getLowerRangeBound( int dimension )
{
	return( -1e308 );
}
		
double sorebFunction_t::getUpperRangeBound( int dimension )
{
	return( 1e308 );
}

void sorebFunction_t::initializeVariableInteractionGraph()
{
	for( int i = 0; i < number_of_subfunctions; i++ )
	{
		int block_start = getStartingIndexOfBlock(i);
		for( int j = 0; j < rotation_block_size; j++ )
		{
			int ind = block_start + j;
			std::set<int> dependent_vars = variable_interaction_graph[ind];
			for( int k = 0; k < rotation_block_size; k++ )
				if( block_start+k != ind )
					dependent_vars.insert(block_start+k);
			variable_interaction_graph[ind] = dependent_vars;
		}
	}
	/*for( auto p : variable_interaction_graph )
	{
		printf("[%d] ",p.first);
		for( int x : p.second )
			printf("%d ",x);
	}
	printf("\n");*/
}
		
sorebFunction_t::~sorebFunction_t()
{
	ezilaitiniObjectiveRotationMatrix( rotation_matrix, rotation_angle, rotation_block_size );
}

osorebFunction_t::osorebFunction_t( int number_of_parameters, double vtr )
{
	this->name = "Overlapping Sum of Rotated Ellipsoid Blocks function";
	this->number_of_parameters = number_of_parameters;
	this->vtr = vtr;
	this->rotation_angle = 45;
	this->rotation_block_size = 5;
	this->number_of_large_rotated_blocks = (number_of_parameters+4)/5;
	this->number_of_small_rotated_blocks = number_of_large_rotated_blocks-1;
	this->number_of_subfunctions = number_of_large_rotated_blocks + number_of_small_rotated_blocks;
	initializeFitnessFunction();
	rotation_matrix_big = initializeObjectiveRotationMatrix( rotation_angle, rotation_block_size );
	rotation_matrix_small = initializeObjectiveRotationMatrix( rotation_angle, 2 );
	//exit(0); // TODO : check evaluation functions
}

void osorebFunction_t::evaluationFunction( solution_t<double> *solution )
{
	assert( black_box_optimization );
	
	double result = 0.0;
	for( int i = 0; i < number_of_large_rotated_blocks; i++ )
	{
		//printf("[%d-%d]\n",rotation_block_size*i,rotation_block_size*i + rotation_block_size - 1);
		result += subfunction( &solution->variables[rotation_block_size*i], rotation_block_size );
	}
	for( int i = 1; i < number_of_small_rotated_blocks; i++ )
	{
		//printf("[%d,%d]\n",rotation_block_size*i-1,rotation_block_size*i);
		result += subfunction( &solution->variables[rotation_block_size*i-1], 2 );
	}

	solution->setObjectiveValue(result);
	solution->setConstraintValue(0);
	full_number_of_evaluations++;
	number_of_evaluations++;
}

void osorebFunction_t::partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution )
{
	int num_subfunctions_evaluated = 0;
	double result = 0.0;
	for( int i = 0; i < solution->getNumberOfTouchedVariables(); i++ )
	{
		int ind = solution->touched_indices[i];
		int block_ind = ind / rotation_block_size;
		if( i > 0 )
		{
			assert( ind > solution->touched_indices[i-1] ); // assume indices are sorted
			int prev_block_ind = solution->touched_indices[i-1] / rotation_block_size;
			if( block_ind == prev_block_ind )
				continue;
		}
		
		int block_start = block_ind * rotation_block_size;
		double *variables_copy = new double[rotation_block_size];
		for( int j = 0; j < rotation_block_size; j++ )
			variables_copy[j] = parent->variables[block_start+j];
		result -= subfunction( variables_copy, rotation_block_size );
		
		double variables_copy_smallblock[2];
		for( int j = 0; j < 2; j++ )
			variables_copy_smallblock[j] = parent->variables[block_start+j-1];
		result -= subfunction( variables_copy_smallblock, 2 );
		
		int j = 0;
		while( i+j < solution->getNumberOfTouchedVariables() && block_ind == solution->touched_indices[i+j] / rotation_block_size )
		{
			int cur_ind = solution->touched_indices[i+j];
			variables_copy[cur_ind % rotation_block_size] = solution->touched_variables[i+j];
			j++;
		}
		result += subfunction( variables_copy, rotation_block_size );
		num_subfunctions_evaluated++;
		
		if( i > 0 )
		{
			if( solution->touched_indices[i-1] == ind-1 )
				variables_copy_smallblock[0] = solution->touched_variables[i-1];
			variables_copy_smallblock[1] = solution->touched_variables[i];
		}
		result += subfunction( variables_copy_smallblock, 2 );
		num_subfunctions_evaluated++;

		delete( variables_copy );
	}
	
	solution->setObjectiveValue(parent->getObjectiveValue() + result);
	solution->setConstraintValue(parent->getConstraintValue());
	full_number_of_evaluations++;
	number_of_evaluations += num_subfunctions_evaluated / (double) number_of_subfunctions;
}

double osorebFunction_t::subfunction( double *vars, int num_vars )
{
	double *rotated_vars;
	if( num_vars == rotation_block_size )
		rotated_vars = rotateVariables( vars, rotation_block_size, rotation_matrix_big );
	else if( num_vars == 2 )
		rotated_vars = rotateVariables( vars, 2, rotation_matrix_small );
	else{printf("Undefined operation\n");exit(0);}
	double result = 0.0;
    for(int i = 0; i < num_vars; i++ )
        result += pow( 10.0, 6.0*(((double) (i))/((double) (num_vars-1))) )*rotated_vars[i]*rotated_vars[i];
	free( rotated_vars );
	return( result ); 
}

double osorebFunction_t::getLowerRangeBound( int dimension )
{
	return( -1e308 );
}
		
double osorebFunction_t::getUpperRangeBound( int dimension )
{
	return( 1e308 );
}

osorebFunction_t::~osorebFunction_t()
{
	if( rotation_angle != 0.0 )
	{
		ezilaitiniObjectiveRotationMatrix( rotation_matrix_big, rotation_angle, rotation_block_size );
		ezilaitiniObjectiveRotationMatrix( rotation_matrix_small, rotation_angle, 2 );
	}
}

sorebChainFunction_t::sorebChainFunction_t( int number_of_parameters, double vtr, double conditioning_number, double rotation_angle, bool wrap_around )
{
	this->name = "Chain of Sum of Rotated Ellipsoid Blocks function";
	this->number_of_parameters = number_of_parameters;
	this->vtr = vtr;
	this->rotation_angle = rotation_angle;
	this->rotation_block_size = 2;
	this->number_of_subfunctions = number_of_parameters-1;
	this->conditioning_number = conditioning_number;
	this->wrap_around = wrap_around;
	initializeFitnessFunction();
	rotation_matrix = initializeObjectiveRotationMatrix( rotation_angle, rotation_block_size );
	if( !black_box_optimization )
		initializeVariableInteractionGraph();
}

void sorebChainFunction_t::evaluationFunction( solution_t<double> *solution )
{
	double result = 0.0;
	for( int i = 0; i < number_of_subfunctions; i++ )
		result += subfunction( &solution->variables[i], rotation_block_size );

	solution->setObjectiveValue(result);
	solution->setConstraintValue(0);
	full_number_of_evaluations++;
	number_of_evaluations++;
}

void sorebChainFunction_t::partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution )
{
	int num_subfunctions_evaluated = 0;
	double result = 0.0;
	for( int i = 0; i < solution->getNumberOfTouchedVariables(); i++ )
	{
		int ind = solution->touched_indices[i];
		double *variables_copy = new double[rotation_block_size];
		if( ind > 0 )
		{
			variables_copy[0] = parent->variables[ind-1];
			variables_copy[1] = parent->variables[ind];
			result -= subfunction( variables_copy, rotation_block_size );
			
			if( i > 0 && solution->touched_indices[i-1] == ind-1 )
				variables_copy[0] = solution->touched_variables[i-1];
			variables_copy[1] = solution->touched_variables[i];
			result += subfunction( variables_copy, rotation_block_size );
			num_subfunctions_evaluated++;
		}
		if( ind < number_of_parameters-1 && !(i < solution->getNumberOfTouchedVariables()-1 && solution->touched_indices[i+1] == ind+1) )
		{
			variables_copy[0] = parent->variables[ind];
			variables_copy[1] = parent->variables[ind+1];
			result -= subfunction( variables_copy, rotation_block_size );
			
			variables_copy[0] = solution->touched_variables[i];
			if( i+1 < solution->getNumberOfTouchedVariables() && solution->touched_indices[i+1] == ind+1 )
				variables_copy[1] = solution->touched_variables[i+1];
			result += subfunction( variables_copy, rotation_block_size );
			num_subfunctions_evaluated++;
		}
		delete[] variables_copy;
	}
	
	solution->setObjectiveValue(parent->getObjectiveValue() + result);
	solution->setConstraintValue(parent->getConstraintValue());
	full_number_of_evaluations++;
	number_of_evaluations += num_subfunctions_evaluated / (double) number_of_subfunctions;
}

double sorebChainFunction_t::subfunction( double *vars, int num_vars )
{
	double *rotated_vars = vars;
	if( rotation_angle != 0.0 )
		rotated_vars = rotateVariables( vars, num_vars, rotation_matrix );
	double result = 0.0;
    for(int i = 0; i < num_vars; i++ )
        result += pow( 10.0, conditioning_number*(((double) (i))/((double) (rotation_block_size-1))) )*rotated_vars[i]*rotated_vars[i];
	if( rotation_angle != 0.0 )
		free( rotated_vars );
	return( result ); 
}

double sorebChainFunction_t::getLowerRangeBound( int dimension )
{
	return( -1e308 );
}
		
double sorebChainFunction_t::getUpperRangeBound( int dimension )
{
	return( 1e308 );
}

void sorebChainFunction_t::initializeVariableInteractionGraph()
{
	for( int i = 0; i < number_of_parameters; i++ )
	{
		std::set<int> dependent_vars;
		if( i > 0 || wrap_around )
			dependent_vars.insert((i+number_of_parameters-1)%number_of_parameters);	
		if( i+1 < number_of_parameters || wrap_around )
			dependent_vars.insert((i+1)%number_of_parameters);	
		variable_interaction_graph[i] = dependent_vars;
	}
	/*for( auto p : variable_interaction_graph )
	{
		printf("[%d] ",p.first);
		for( int x : p.second )
			printf("%d ",x);
	}
	printf("\n");*/
}

sorebChainFunction_t::~sorebChainFunction_t()
{
	ezilaitiniObjectiveRotationMatrix( rotation_matrix, rotation_angle, rotation_block_size );
}

sorebGridFunction_t::sorebGridFunction_t( int number_of_parameters, double vtr, double conditioning_number, double rotation_angle, bool wrap_around_x, bool wrap_around_y )
{
	this->name = "Chain of Sum of Rotated Ellipsoid Blocks function";
	this->number_of_parameters = number_of_parameters;
	this->vtr = vtr;
	this->rotation_angle = rotation_angle;
	this->conditioning_number = conditioning_number;
	this->wrap_around_x = wrap_around_x;
	this->wrap_around_y = wrap_around_y;
	this->number_of_subfunctions = number_of_parameters;
	this->grid_width = round(sqrt(number_of_parameters));
	assert( grid_width * grid_width == number_of_parameters );
	//if( !wrap_around_x ) this->number_of_subfunctions -= grid_width;
	//if( !wrap_around_y ) this->number_of_subfunctions -= grid_width;
	initializeFitnessFunction();
	if( !black_box_optimization )
		initializeVariableInteractionGraph();
}

sorebGridFunction_t::~sorebGridFunction_t()
{
	for( auto it : rotation_matrices )
	{
		int n = it.first;
		double **rot_mat = it.second;
		for( int i = 0; i < n; i++ )
			free( rot_mat[i] );
		free( rot_mat );
	}
}

void sorebGridFunction_t::evaluationFunction( solution_t<double> *solution )
{
	double result = 0.0;
	double *var_tmp = new double[5];
	for( int y = 0; y < grid_width; y++ )
	{
		for( int x = 0; x < grid_width; x++ )
		{
			int num_vars = 0;
			int ind = y*grid_width + x;
			var_tmp[num_vars++] = solution->variables[ind];
			//result += solution->variables[ind] * solution->variables[ind];
			//printf("x%d**2 + ",ind);
			if( x+1 < grid_width || wrap_around_x )
			{
				int ind_right = y*grid_width + ((x+1)%grid_width);
				var_tmp[num_vars++] = solution->variables[ind_right];
				//result += fabs(2 * solution->variables[ind] * solution->variables[ind_right]);
				//printf("2*x%d*x%d + ",ind,ind_right);
			}
			if( y+1 < grid_width || wrap_around_y )
			{
				int ind_down = ((y+1)%grid_width)*grid_width + x;
				var_tmp[num_vars++] = solution->variables[ind_down];
				//result += fabs(2 * solution->variables[ind] * solution->variables[ind_down]);
				//printf("2*x%d*x%d + ",ind,ind_down);
			}
			if( x > 0 || wrap_around_x )
			{
				int ind_left = y*grid_width + ((x+grid_width-1)%grid_width);
				var_tmp[num_vars++] = solution->variables[ind_left];
				//result += solution->variables[ind] * solution->variables[ind_left];
				//printf("x%d * x%d + ",ind,ind_left);
			}
			if( y > 0 || wrap_around_y )
			{
				int ind_up = ((y+grid_width-1)%grid_width)*grid_width + x;
				var_tmp[num_vars++] = solution->variables[ind_up];
				//result += solution->variables[ind] * solution->variables[ind_up];
				//printf("x%d * x%d + ",ind,ind_up);
			}
			//printf("%d-D ellipsoid\n",num_vars);
			result += subfunction( var_tmp, num_vars );
		}
	}
	//exit(0);
	delete[] var_tmp;

	solution->setObjectiveValue(result);
	//printf("\n");
	//exit(0);
	solution->setConstraintValue(0);
	full_number_of_evaluations++;
	number_of_evaluations++;
}

void sorebGridFunction_t::partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution )
{
	/*double result = 0.0;
	double *variables_copy = new double[2];
	int num_subfunctions_evaluated = 0;
	for( int i = 0; i < solution->getNumberOfTouchedVariables(); i++ )
	{
		int ind = solution->touched_indices[i];
		int x = ind%grid_width;
		int y = ind/grid_width;
		
		int ind_left = y*grid_width + (x+grid_width-1)%grid_width;
		int ind_right = y*grid_width + (x+1)%grid_width;
		int ind_up = ((y+grid_width-1)%grid_width)*grid_width + x;
		int ind_down = ((y+1)%grid_width)*grid_width + x;
		int ind_touched_left = solution->getTouchedIndex(ind_left);
		int ind_touched_up = solution->getTouchedIndex(ind_up);
		int ind_touched_right = solution->getTouchedIndex(ind_right);
		int ind_touched_down = solution->getTouchedIndex(ind_down);

		if( ind_touched_left == -1 && (x > 0 || wrap_around_x) ) // ind_left was not touched; if it was, this subfunction was already recomputed
		{
			// update left subfunction
			variables_copy[0] = parent->variables[ind_left];
			variables_copy[1] = parent->variables[ind];
			result -= subfunction( variables_copy, 2 );
			
			variables_copy[1] = solution->touched_variables[i];
			result += subfunction( variables_copy, 2 );
			num_subfunctions_evaluated++;
		}
		if( ind_touched_up == -1 && (y > 0 || wrap_around_y) ) // ind_up was not touched; if it was, this subfunction was already recomputed
		{
			// update up subfunction
			variables_copy[0] = parent->variables[ind_up];
			variables_copy[1] = parent->variables[ind];
			result -= subfunction( variables_copy, 2 );
			
			variables_copy[1] = solution->touched_variables[i];
			result += subfunction( variables_copy, 2 );
			num_subfunctions_evaluated++;
		}

		if( x+1 < grid_width || wrap_around_x )
		{
			// subtract old right subfunction
			variables_copy[0] = parent->variables[ind];
			variables_copy[1] = parent->variables[ind_right];
			result -= subfunction( variables_copy, 2 );

			// add new right subfunction
			variables_copy[0] = solution->touched_variables[i];
			if( ind_touched_right == -1 )
				variables_copy[1] = parent->variables[ind_right];
			else
				variables_copy[1] = solution->touched_variables[ind_touched_right];
			result += subfunction( variables_copy, 2 );
			num_subfunctions_evaluated++;
		}

		if( y+1 < grid_width || wrap_around_y )
		{
			// subtract old down subfunction
			variables_copy[0] = parent->variables[ind];
			variables_copy[1] = parent->variables[ind_down];
			result -= subfunction( variables_copy, 2 );

			// add new down subfunction
			variables_copy[0] = solution->touched_variables[i];
			if( ind_touched_down  == -1 )
				variables_copy[1] = parent->variables[ind_down];
			else
				variables_copy[1] = solution->touched_variables[ind_touched_down];
			result += subfunction( variables_copy, 2 );
			num_subfunctions_evaluated++;
		}

		delete[] variables_copy;
	}
	solution->setObjectiveValue(parent->getObjectiveValue() + result);
	solution->setConstraintValue(parent->getConstraintValue());
	full_number_of_evaluations++;
	number_of_evaluations += num_subfunctions_evaluated / (double) number_of_subfunctions;*/
	
	std::set<int> subfunction_indices;
	for( int i = 0; i < solution->getNumberOfTouchedVariables(); i++ )
	{
		int ind = solution->touched_indices[i];
		subfunction_indices.insert(ind);
		for( int x : getNeighborsInGrid(ind) )
			subfunction_indices.insert(x);
	}

	int num_subfunctions_evaluated;
   	if( solution->getNumberOfTouchedVariables() == 1 )
	{
		num_subfunctions_evaluated = 1 + getNeighborsInGrid( solution->touched_indices[0] ).size();
		assert( num_subfunctions_evaluated == subfunction_indices.size() );
	}
	else
		num_subfunctions_evaluated = subfunction_indices.size();

	assert( num_subfunctions_evaluated <= number_of_subfunctions );
	evaluatePartialSolutionBlackBox( parent, solution );
	number_of_evaluations--;
	number_of_evaluations += num_subfunctions_evaluated / (double) number_of_subfunctions;
}

double sorebGridFunction_t::subfunction( double *vars, int num_vars )
{
	double *rotated_vars = vars;
	if( rotation_angle != 0.0 )
	{
		if( rotation_matrices.find(num_vars) == rotation_matrices.end() )
			rotation_matrices[num_vars] = initializeObjectiveRotationMatrix( rotation_angle, num_vars );
		double **rotation_matrix = rotation_matrices.find(num_vars)->second;
		rotated_vars = rotateVariables( vars, num_vars, rotation_matrix );
	}
	double result = 0.0;
	for(int i = 0; i < num_vars; i++ )
		result += pow( 10.0, conditioning_number*(((double) (i))/((double) (num_vars-1))) )*(rotated_vars[i])*(rotated_vars[i]);
	if( rotation_angle != 0.0 )
		free( rotated_vars );
	return( result ); 
}

double sorebGridFunction_t::getLowerRangeBound( int dimension )
{
	return( -1e308 );
}
		
double sorebGridFunction_t::getUpperRangeBound( int dimension )
{
	return( 1e308 );
}

std::set<int> sorebGridFunction_t::getNeighborsInGrid( int ind )
{
	int j = ind % grid_width;
	int i = ind / grid_width;
	int ind_up = (ind+number_of_parameters-grid_width)%number_of_parameters;
	int ind_left = (i*grid_width)+((j+grid_width-1)%grid_width);
	int ind_right = (i*grid_width)+((j+1)%grid_width);
	int ind_down = (ind+grid_width)%number_of_parameters;
	
	std::set<int> dependent_vars;
	if( ind_up < ind || wrap_around_y )
		dependent_vars.insert(ind_up);	
	if( ind_left < ind || wrap_around_x )
		dependent_vars.insert(ind_left);
	if( ind_down > ind || wrap_around_y )
		dependent_vars.insert(ind_down);
	if( ind_right > ind || wrap_around_x )
		dependent_vars.insert(ind_right);	
	return( dependent_vars );
}

void sorebGridFunction_t::initializeVariableInteractionGraph()
{
	for( int i = 0; i < grid_width; i++ )
	{
		for( int j = 0; j < grid_width; j++ )
		{
			int ind = i*grid_width+j;
			variable_interaction_graph[ind] = getNeighborsInGrid(ind);
		}
	}
	
	// Add an edge to neighbors of neighbors
	for( int i = 0; i < grid_width; i++ )
	{
		for( int j = 0; j < grid_width; j++ )
		{
			int ind = i*grid_width+j;
			std::set<int> dependent_vars = variable_interaction_graph[ind];
			std::set<int> neighbors = getNeighborsInGrid(ind);
			for( int x : neighbors )
			{
				for( int y : getNeighborsInGrid(x) )
					if( y != ind ) dependent_vars.insert(y);
			}				
			variable_interaction_graph[ind] = dependent_vars;
		}
	}
	/*for( auto p : variable_interaction_graph )
	{
		printf("[%d] ",p.first);
		for( int x : p.second )
			printf("%d ",x);
	}
	printf("\n");*/
	
}

sorebCubeFunction_t::sorebCubeFunction_t( int number_of_parameters, double vtr, double conditioning_number, double rotation_angle, bool wrap_around_x, bool wrap_around_y, bool wrap_around_z )
{
	this->name = "Chain of Sum of Rotated Ellipsoid Blocks function";
	this->number_of_parameters = number_of_parameters;
	this->vtr = vtr;
	this->rotation_angle = rotation_angle;
	this->conditioning_number = conditioning_number;
	this->wrap_around_x = wrap_around_x;
	this->wrap_around_y = wrap_around_y;
	this->wrap_around_z = wrap_around_z;
	this->number_of_subfunctions = number_of_parameters;
	this->cube_width = round(cbrt(number_of_parameters));
	assert( cube_width * cube_width * cube_width == number_of_parameters );
	//if( !wrap_around_x ) this->number_of_subfunctions -= cube_width;
	//if( !wrap_around_y ) this->number_of_subfunctions -= cube_width;
	//if( !wrap_around_z ) this->number_of_subfunctions -= cube_width;
	initializeFitnessFunction();
	rotation_matrix = initializeObjectiveRotationMatrix( rotation_angle, rotation_block_size );
	if( !black_box_optimization )
		initializeVariableInteractionGraph();
}

std::set<int> sorebCubeFunction_t::getNeighborsInGrid( int ind )
{
	int x = ind % cube_width;
	ind /= cube_width;
	int y = ind % cube_width;
	ind /= cube_width;
	int z = ind % cube_width;
	int ind_xnext = z*cube_width*cube_width + y*cube_width + (x+1)%cube_width;

	int ind_ynext = z*cube_width*cube_width + ((y+1)%cube_width)*cube_width + x;
	int ind_znext = ((z+1)%cube_width)*cube_width*cube_width + y*cube_width + x;
	int ind_xprev = z*cube_width*cube_width + y*cube_width + (x+cube_width-1)%cube_width;
	int ind_yprev = z*cube_width*cube_width + ((y+cube_width-1)%cube_width)*cube_width + x;
	int ind_zprev = ((z+cube_width-1)%cube_width)*cube_width*cube_width + y*cube_width + x;

	std::set<int> dependent_vars;
	if( x > 0 || wrap_around_x )
		dependent_vars.insert(ind_xprev);
	if( x+1 < cube_width || wrap_around_x )
		dependent_vars.insert(ind_xnext);
	if( y > 0 || wrap_around_y )
		dependent_vars.insert(ind_yprev);
	if( y+1 < cube_width || wrap_around_y )
		dependent_vars.insert(ind_ynext);
	if( z > 0 || wrap_around_z )
		dependent_vars.insert(ind_zprev);
	if( z+1 < cube_width || wrap_around_z )
		dependent_vars.insert(ind_znext);

	return( dependent_vars );
}

void sorebCubeFunction_t::initializeVariableInteractionGraph()
{
	for( int z = 0; z < cube_width; z++ )
	{
		for( int y = 0; y < cube_width; y++ )
		{
			for( int x = 0; x < cube_width; x++ )
			{
				int ind = z*cube_width*cube_width + y*cube_width + x;
				std::set<int> neighbors = getNeighborsInGrid( ind );
				std::set<int> dependent_vars;
				for( int x : neighbors )
				{
					dependent_vars.insert(x);
					for( int y : getNeighborsInGrid( x ) )
						if( y != ind ) dependent_vars.insert(y);
				}
				variable_interaction_graph[ind] = dependent_vars;
			}
		}
	}
	/*for( auto p : variable_interaction_graph )
	{
		printf("[%d] ",p.first);
		for( int x : p.second )
			printf("%d ",x);
	}
	printf("\n");*/
	
}

void sorebCubeFunction_t::evaluationFunction( solution_t<double> *solution )
{
	double result = 0.0;
	double *var_tmp = new double[7];
	for( int z = 0; z < cube_width; z++ )
	{
		for( int y = 0; y < cube_width; y++ )
		{
			for( int x = 0; x < cube_width; x++ )
			{
				int num_vars = 0;
				int ind = z*cube_width*cube_width + y*cube_width + x;
				var_tmp[num_vars++] = solution->variables[ind];

				int ind_xprev = z*cube_width*cube_width + y*cube_width + (x+cube_width-1)%cube_width;
				if( x > 0 || wrap_around_x )
					var_tmp[num_vars++] = solution->variables[ind_xprev];
				int ind_xnext = z*cube_width*cube_width + y*cube_width + (x+1)%cube_width;
				if( x+1 < cube_width || wrap_around_x )
					var_tmp[num_vars++] = solution->variables[ind_xnext];

				int ind_yprev = z*cube_width*cube_width + ((y+cube_width-1)%cube_width)*cube_width + x;
				if( y > 0 || wrap_around_y )
					var_tmp[num_vars++] = solution->variables[ind_yprev];
				int ind_ynext = z*cube_width*cube_width + ((y+1)%cube_width)*cube_width + x;
				if( y+1 < cube_width || wrap_around_y )
					var_tmp[num_vars++] = solution->variables[ind_ynext];

				int ind_zprev = ((z+cube_width-1)%cube_width)*cube_width*cube_width + y*cube_width + x;
				if( z > 0 || wrap_around_z )
					var_tmp[num_vars++] = solution->variables[ind_zprev];
				int ind_znext = ((z+1)%cube_width)*cube_width*cube_width + y*cube_width + x;
				if( z+1 < cube_width || wrap_around_z )
					var_tmp[num_vars++] = solution->variables[ind_znext];
				
				result += subfunction( var_tmp, num_vars );
			}
		}
	}
	delete[] var_tmp;

	solution->setObjectiveValue(result);
	solution->setConstraintValue(0);
	full_number_of_evaluations++;
	number_of_evaluations++;
}

void sorebCubeFunction_t::partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution )
{
	/*int num_subfunctions_evaluated = 0;
	double result = 0.0;
	double *variables_copy = new double[2];
	for( int i = 0; i < solution->getNumberOfTouchedVariables(); i++ )
	{
		int ind = solution->touched_indices[i];
		
		int indc = ind;
		int x = indc%cube_width;
		indc /= cube_width;
		int y = indc%cube_width;
		indc /= cube_width;
		int z = indc%cube_width;
		
		int ind_xnext = z*cube_width*cube_width + y*cube_width + (x+1)%cube_width;
		int ind_ynext = z*cube_width*cube_width + ((y+1)%cube_width)*cube_width + x;
		int ind_znext = ((z+1)%cube_width)*cube_width*cube_width + y*cube_width + x;
		int ind_xprev = z*cube_width*cube_width + y*cube_width + (x+cube_width-1)%cube_width;
		int ind_yprev = z*cube_width*cube_width + ((y+cube_width-1)%cube_width)*cube_width + x;
		int ind_zprev = ((z+cube_width-1)%cube_width)*cube_width*cube_width + y*cube_width + x;
		int ind_touched_xnext = solution->getTouchedIndex(ind_xnext);
		int ind_touched_ynext = solution->getTouchedIndex(ind_ynext);
		int ind_touched_znext = solution->getTouchedIndex(ind_znext);
		int ind_touched_xprev = solution->getTouchedIndex(ind_xprev);
		int ind_touched_yprev = solution->getTouchedIndex(ind_yprev);
		int ind_touched_zprev = solution->getTouchedIndex(ind_zprev);

		if( ind_touched_xprev == -1 && (x > 0 || wrap_around_x) ) // ind was not touched; if it was, this subfunction was already recomputed
		{
			variables_copy[0] = parent->variables[ind_xprev];
			variables_copy[1] = parent->variables[ind];
			result -= subfunction( variables_copy, 2 );
			
			variables_copy[1] = solution->touched_variables[i];
			result += subfunction( variables_copy, 2 );
			num_subfunctions_evaluated++;
		}
		if( ind_touched_yprev == -1 && (y > 0 || wrap_around_y) ) // ind_left was not touched; if it was, this subfunction was already recomputed
		{
			variables_copy[0] = parent->variables[ind_yprev];
			variables_copy[1] = parent->variables[ind];
			result -= subfunction( variables_copy, 2 );
			
			variables_copy[1] = solution->touched_variables[i];
			result += subfunction( variables_copy, 2 );
			num_subfunctions_evaluated++;
		}
		if( ind_touched_zprev == -1 && (z > 0 || wrap_around_z) ) // ind_left was not touched; if it was, this subfunction was already recomputed
		{
			variables_copy[0] = parent->variables[ind_zprev];
			variables_copy[1] = parent->variables[ind];
			result -= subfunction( variables_copy, 2 );
			
			variables_copy[1] = solution->touched_variables[i];
			result += subfunction( variables_copy, 2 );
			num_subfunctions_evaluated++;
		}

		if( x+1 < cube_width || wrap_around_x )
		{
			// subtract old subfunction for xnext
			variables_copy[0] = parent->variables[ind];
			variables_copy[1] = parent->variables[ind_xnext];
			result -= subfunction( variables_copy, 2 );

			// add new subfunction for xnext
			variables_copy[0] = solution->touched_variables[i];
			if( ind_touched_xnext == -1 )
				variables_copy[1] = parent->variables[ind_xnext];
			else
				variables_copy[1] = solution->touched_variables[ind_touched_xnext];
			result += subfunction( variables_copy, 2 );
			num_subfunctions_evaluated++;
		}

		if( y+1 < cube_width || wrap_around_y )
		{
			// subtract old subfunction for ynext
			variables_copy[0] = parent->variables[ind];
			variables_copy[1] = parent->variables[ind_ynext];
			result -= subfunction( variables_copy, 2 );

			// add new subfunction for ynext
			variables_copy[0] = solution->touched_variables[i];
			if( ind_touched_ynext == -1 )
				variables_copy[1] = parent->variables[ind_ynext];
			else
				variables_copy[1] = solution->touched_variables[ind_touched_ynext];
			result += subfunction( variables_copy, 2 );
			num_subfunctions_evaluated++;
		}

		if( z+1 < cube_width || wrap_around_z )
		{
			// subtract old subfunction for znext
			variables_copy[0] = parent->variables[ind];
			variables_copy[1] = parent->variables[ind_znext];
			result -= subfunction( variables_copy, 2 );

			// add new subfunction for znext
			variables_copy[0] = solution->touched_variables[i];
			if( ind_touched_znext == -1 )
				variables_copy[1] = parent->variables[ind_znext];
			else
				variables_copy[1] = solution->touched_variables[ind_touched_znext];
			result += subfunction( variables_copy, 2 );
			num_subfunctions_evaluated++;
		}

		delete[] variables_copy;
	}
	solution->setObjectiveValue(parent->getObjectiveValue() + result);
	solution->setConstraintValue(parent->getConstraintValue());
	full_number_of_evaluations++;
	number_of_evaluations += num_subfunctions_evaluated / (double) number_of_subfunctions;*/
	std::set<int> subfunction_indices;
	for( int i = 0; i < solution->getNumberOfTouchedVariables(); i++ )
	{
		int ind = solution->touched_indices[i];
		subfunction_indices.insert(ind);
		for( int x : getNeighborsInGrid(ind) )
			subfunction_indices.insert(x);
	}

	int num_subfunctions_evaluated;
   	if( solution->getNumberOfTouchedVariables() == 1 )
	{
		num_subfunctions_evaluated = 1 + getNeighborsInGrid( solution->touched_indices[0] ).size();
		assert( num_subfunctions_evaluated == subfunction_indices.size() );
	}
	else
		num_subfunctions_evaluated = subfunction_indices.size();
	
	assert( num_subfunctions_evaluated <= number_of_subfunctions );
	evaluatePartialSolutionBlackBox( parent, solution );
	number_of_evaluations--;
	number_of_evaluations += num_subfunctions_evaluated / (double) number_of_subfunctions;
}

double sorebCubeFunction_t::subfunction( double *vars, int num_vars )
{
	double *rotated_vars = vars;
	if( rotation_angle != 0.0 )
	{
		if( rotation_matrices.find(num_vars) == rotation_matrices.end() )
			rotation_matrices[num_vars] = initializeObjectiveRotationMatrix( rotation_angle, num_vars );
		double **rotation_matrix = rotation_matrices.find(num_vars)->second;
		rotated_vars = rotateVariables( vars, num_vars, rotation_matrix );
	}
	double result = 0.0;
    for(int i = 0; i < num_vars; i++ )
        result += pow( 10.0, conditioning_number*(((double) (i))/((double) (num_vars-1))) )*rotated_vars[i]*rotated_vars[i];
	if( rotation_angle != 0.0 )
		free( rotated_vars );
	return( result ); 
}

double sorebCubeFunction_t::getLowerRangeBound( int dimension )
{
	return( -1e308 );
}
		
double sorebCubeFunction_t::getUpperRangeBound( int dimension )
{
	return( 1e308 );
}

sorebCubeFunction_t::~sorebCubeFunction_t()
{
	for( auto it : rotation_matrices )
	{
		int n = it.first;
		double **rot_mat = it.second;
		for( int i = 0; i < n; i++ )
			free( rot_mat[i] );
		free( rot_mat );
	}
}
		
}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
