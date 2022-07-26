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
#include "fitness.hpp"
#include "embed.hpp"
#include "RealValuedGOMEA.h"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
namespace gomea{
namespace realvalued{

fitness_t::~fitness_t()
{
 	free( lower_range_bound );
    free( upper_range_bound );
}

fitness_t *fitness_t::getFitnessClass( int problem_index, int number_of_parameters, double vtr )
{
	if( problem_index > 1000 )
	{
		double rotation_angle = (problem_index%10)*5;
		problem_index/=10;
		double conditioning_number = problem_index%10;
		problem_index/=10;
		int overlap_size = problem_index%10;
		problem_index/=10;
		int block_size = problem_index;
		//printf("%d %d %d %d\n",rotation_angle,conditioning_number,overlap_size,block_size);
		fitness_t *func = new sorebFunction_t( number_of_parameters, vtr, conditioning_number, rotation_angle, block_size, overlap_size );
		return( func );
	}

#ifdef CECLSGOFUNC
	if( problem_index > 200 && problem_index < 300 )
		return( new CECLSGOFunctions_t( problem_index-200, number_of_parameters, vtr ) );
#endif

	switch( problem_index )
	{
		case 0 : return( new sphereFunction_t( number_of_parameters, vtr ) );
		case 7 : return( new rosenbrockFunction_t( number_of_parameters, vtr ) );
		case 13: return( new sorebFunction_t( number_of_parameters, vtr, 6, 45, 5, 0 ) );
		case 16: return( new osorebFunction_t( number_of_parameters, vtr ) );
		case 17: return( new BD2FunctionHypervolume_t( number_of_parameters, vtr ) );
		case 10: return( new sorebChainFunction_t( number_of_parameters, vtr, 6, -45, 0 ) );
		case 20: return( new sorebGridFunction_t( number_of_parameters, vtr, 6, -45, 0, 0 ) );
		case 21: return( new sorebGridFunction_t( number_of_parameters, vtr, 6, -45, 1, 0 ) );
		case 22: return( new sorebGridFunction_t( number_of_parameters, vtr, 6, -45, 1, 1 ) );
		case 30: return( new sorebCubeFunction_t( number_of_parameters, vtr, 6, -45, 0, 0, 0 ) );
		case 31: return( new sorebCubeFunction_t( number_of_parameters, vtr, 6, -45, 1, 0, 0 ) );
		case 32: return( new sorebCubeFunction_t( number_of_parameters, vtr, 6, -45, 1, 1, 0 ) );
		case 33: return( new sorebCubeFunction_t( number_of_parameters, vtr, 6, -45, 1, 1, 1 ) );
		default: return NULL;
	}
}

/**
 * Returns 1 if x is better than y, 0 otherwise.
 * x is not better than y unless:
 * - x and y are both infeasible and x has a smaller sum of constraint violations, or
 * - x is feasible and y is not, or
 * - x and y are both feasible and x has a smaller objective value than y
 */
short fitness_t::betterFitness( double objective_value_x, double constraint_value_x, double objective_value_y, double constraint_value_y )
{
    short result = 0;

    if( constraint_value_x > 0 ) /* x is infeasible */
    {
        if( constraint_value_y > 0 ) /* Both are infeasible */
        {
            if( constraint_value_x < constraint_value_y )
                result = 1;
        }
    }
    else /* x is feasible */
    {
        if( constraint_value_y > 0 ) /* x is feasible and y is not */
            result = 1;
        else /* Both are feasible */
        {
            if( objective_value_x < objective_value_y )
                result = 1;
        }
    }

    return( result );
}

short fitness_t::betterFitness( solution_t *sol_x, solution_t *sol_y ) 
{
	return( betterFitness( sol_x->objective_value, sol_x->constraint_value, sol_y->objective_value, sol_y->constraint_value ) );
}

void fitness_t::initializeFitnessFunction( void )
{
	initializeRangeBounds();
}

// distance to a box defined by [-infty, ref_x, -infty, ref_y]
double fitness_t::distance_to_box(double ref_x, double ref_y, double p_x, double p_y)
{
	double dx = max(0.0, p_x - ref_x );
	double dy = max(0.0, p_y - ref_y );
	return sqrt(dx*dx + dy*dy);
}

// Based on the Uncrowded Hypervolume improvement by the Inria group,
// but we extened the definition to points that are outside of the reference frame
// we compute the distance to the non-dominated area, within the reference window (r_x,r_y)
// define the area points a(P^(i)_x, P^(i-1)_y), for i = 0...n (n =  number of points in front, i.e., obj0.size())
// and let P^(-1)_y = r_y,   and    P^(n)_x = r_x
double fitness_t::distance_to_front(double p_x, double p_y, const std::vector<double>& obj_x, const std::vector<double>& obj_y, std::vector<size_t> & sorted_obj, double r_x, double r_y)
{
	// if the front is empty, use the reference point for the distance measure
	if(obj_x.size() == 0) {
		return distance_to_box(r_x, r_y, p_x, p_y);
	}

	size_t n = obj_x.size();

	// if not available, get the sorted front
	if(sorted_obj.size() != n)
	{
		sorted_obj.resize(n);
		for (size_t i = 0; i < n; ++i) {
			sorted_obj[i] = i;
		}

		std::sort(std::begin(sorted_obj), std::end(sorted_obj), [&obj_x](double idx, double idy) { return obj_x[(size_t)idx] < obj_x[(size_t)idy]; });
	}
	
	double dist;

	// distance to the 'end' boxes
	double min_dist = min( distance_to_box(obj_x[sorted_obj[0]], r_y, p_x, p_y), distance_to_box(r_x, obj_y[sorted_obj[n-1]], p_x, p_y) );

	// distance to 'inner' boxes
	for(size_t k = 1; k < n; ++k)
	{
		dist = distance_to_box(obj_x[sorted_obj[k]], obj_y[sorted_obj[k-1]], p_x, p_y);
		if(dist < min_dist) {
			min_dist = dist;
		}
	}
	assert(min_dist >= 0); // can be 0 if its at the front!
	return min_dist;
}

bool fitness_t::paretoDominates2D( double xf0, double xf1, double yf0, double yf1 )
{
	bool res = true;
	if( yf0 < xf0 ) res = false;
	else if( yf1 < xf1 ) res = false;
	else if( yf0 == xf0 && yf1 == xf1 ) res = false;
	return( res );
}

double fitness_t::compute2DUncrowdedHypervolume( double *obj_f0, double *obj_f1, int population_size )
{
	std::vector<double> f0_approx_set;
	std::vector<double> f1_approx_set;
	std::vector<double> f0_dominated_set;
	std::vector<double> f1_dominated_set;
	double rx = 1.1, ry = 1.1;
	for( int i = 0; i < population_size; i++ )
	{
		bool in_approx_set = false;
		if( obj_f0[i] < rx && obj_f1[i] < ry )
		{
			in_approx_set = true;
			for( int j = 0; j < f0_approx_set.size(); j++ )
			{
				if( paretoDominates2D( f0_approx_set[j], f1_approx_set[j], obj_f0[i], obj_f1[i] ) )
				{
					in_approx_set = true;
					break;
				}
			}
		}
		if( in_approx_set )
		{
			f0_approx_set.push_back(obj_f0[i]);
			f1_approx_set.push_back(obj_f1[i]);
		}
		else
		{
			f0_dominated_set.push_back(obj_f0[i]);
			f1_dominated_set.push_back(obj_f1[i]);
		}
	}

	std::vector<size_t> sorted_obj;
	double penalty = 0;
	for( int i = 0; i < f0_dominated_set.size(); i++ )
		penalty += distance_to_front( f0_dominated_set[i], f1_dominated_set[i], f0_approx_set, f1_approx_set, sorted_obj, rx, ry);
	if( f0_dominated_set.size() > 0 )
		penalty /= f0_dominated_set.size();
	double hv = compute2DHyperVolume( f0_approx_set, f1_approx_set, sorted_obj, rx, ry );
	double res = -hv + penalty;
	//printf("P[%d] = %10.3e\tHV[%d] = %10.3e\tR = %10.3e\n",f0_dominated_set.size(),penalty,f0_approx_set.size(),hv,res);
	return( res );
}

double fitness_t::compute2DHyperVolume( const std::vector<double> &obj_f0, const std::vector<double> &obj_f1, std::vector<size_t> &sorted, double rx, double ry )
{
    int n = obj_f0.size();
	if( n == 0 ) return 0.0;
	if(sorted.size() != n)
	{
		sorted.resize(n);
		for (size_t i = 0; i < n; ++i) {
			sorted[i] = i;
		}
		std::sort(std::begin(sorted), std::end(sorted), [&obj_f0](double idx, double idy) { return obj_f0[(size_t)idx] < obj_f0[(size_t)idy]; });
	}

    double area = (rx - fmin(rx, obj_f0[sorted[n-1]])) * (ry - fmin(ry, obj_f1[sorted[n-1]]));
    for( int i = n-2; i >= 0; i-- )
        area += (fmin(rx, obj_f0[sorted[i+1]]) - fmin(rx, obj_f0[sorted[i]])) * (ry-fmin(ry, obj_f1[sorted[i]]));

    return area;
}

double fitness_t::compute2DHyperVolume( double *obj_f0, double *obj_f1, int population_size )
{
    int i, n, *sorted;
    double max_0, max_1, area;

    n = population_size;
    max_0 = 1.1;
    max_1 = 1.1;
    sorted = mergeSort( obj_f0, n );

    area = (max_0 - fmin(max_0, obj_f0[sorted[n-1]])) * (max_1 - fmin(max_1, obj_f1[sorted[n-1]]));
    for( i = n-2; i >= 0; i-- )
        area += (fmin(max_0, obj_f0[sorted[i+1]]) - fmin(max_0, obj_f0[sorted[i]])) * (max_1-fmin(max_1, obj_f1[sorted[i]]));

    free( sorted );

    return area;
}

double *fitness_t::rotateVariables( double *variables, int num_variables, double **rotation_matrix )
{
	double *rotated_variables = matrixVectorMultiplication( rotation_matrix, variables, num_variables, num_variables );
    return( rotated_variables );
}

double *fitness_t::rotateVariablesInBlocks( double *variables, int len, int from, int to, double **rotation_matrix )
{
	assert( len % rotation_block_size == 0 );
    int num_blocks = len / rotation_block_size;
	double *rotated_variables = (double*) Malloc( len*sizeof( double ) );
    for( int i = 0; i < from; i++ ) rotated_variables[i] = variables[i];
	double *cluster = (double*) Malloc( rotation_block_size*sizeof( double ) );
    for( int i = 0; i < num_blocks; i++ )
    {
        for( int j = 0; j < rotation_block_size; j++ )
            cluster[j] = variables[from + i*rotation_block_size + j];
        double *rotated_cluster = matrixVectorMultiplication( rotation_matrix, cluster, rotation_block_size, rotation_block_size );
        for( int j = 0; j < rotation_block_size; j++ )
            rotated_variables[from + i*rotation_block_size + j] = rotated_cluster[j];
        free( rotated_cluster );
    }
	free( cluster );
    for( int i = to+1; i < len; i++ ) rotated_variables[i] = variables[i];
    return( rotated_variables );
}

void fitness_t::evaluate( solution_t *solution )
{
	if( !gomealib::utils::embeddingInitialized() )
	{
		int out = gomealib::utils::initializePythonEmbedding("RealValuedGOMEA",PyInit_RealValuedGOMEA);
		assert( out == 0 );
	}
	
	evaluationFunction( solution );
	
	int out = evaluationEmbedded();
	assert( out == 0 );	
    
	if( use_vtr && !vtr_hit_status && solution->constraint_value == 0 && solution->objective_value <= vtr  )
	{
		vtr_hit_status = 1;
		elitist_objective_value = solution->objective_value;
		elitist_constraint_value = solution->constraint_value;
	}

	if( !vtr_hit_status && betterFitness(solution->objective_value, solution->constraint_value, elitist_objective_value, elitist_constraint_value) )
	{
		elitist_objective_value = solution->objective_value;
		elitist_constraint_value = solution->constraint_value;
	}
}

void fitness_t::evaluatePartialSolutionBlackBox( solution_t *parent, partial_solution_t *solution )
{
	// Make backup of parent
	double *var_backup = (double*) Malloc(solution->num_touched_variables * sizeof(double));
	for( int i = 0; i < solution->num_touched_variables; i++ )
	{
		var_backup[i] = parent->variables[solution->touched_indices[i]];
		parent->variables[solution->touched_indices[i]] = solution->touched_variables[i];
	}
	double obj_val_backup = parent->objective_value;
	double cons_val_backup = parent->constraint_value;

	evaluationFunction( parent );

	// Insert calculated objective and constraint values into partial solution	
	solution->objective_value = parent->objective_value;
	solution->constraint_value = parent->constraint_value;

	// Restore parent to original state
	parent->objective_value = obj_val_backup;
	parent->constraint_value = cons_val_backup;
	for( int i = 0; i < solution->num_touched_variables; i++ )
		parent->variables[solution->touched_indices[i]] = var_backup[i];
	free( var_backup );
}

void fitness_t::evaluatePartialSolution( solution_t *parent, partial_solution_t *solution )
{
	if( black_box_optimization || solution->num_touched_variables == number_of_parameters )
	{
		evaluatePartialSolutionBlackBox( parent, solution );
	}
	else
	{
		partialEvaluationFunction( parent, solution );
#ifdef CHECK_PARTIAL_FITNESS
		double fbefore = solution->objective_value;
		evaluatePartialSolutionBlackBox( parent, solution );	
		double fafter = solution->objective_value;
		if( fabs((fbefore-fafter)/fafter) > 1e-4 )
		{
			printf("fbefore = %10.3e; ",fbefore);
			printf("fafter = %10.3e\n",fafter);
			exit(0);
		}
		solution->objective_value = fbefore;
		number_of_evaluations--;
#endif

		if( use_vtr && !vtr_hit_status && solution->constraint_value == 0 && solution->objective_value <= vtr  )
		{
			evaluatePartialSolutionBlackBox( parent, solution );
			if( solution->constraint_value == 0 && solution->objective_value <= vtr  )
			{
				vtr_hit_status = 1;
				elitist_objective_value = solution->objective_value;
				elitist_constraint_value = solution->constraint_value;
			}
		}
	}
	
	if( !vtr_hit_status && betterFitness(solution->objective_value, solution->constraint_value, elitist_objective_value, elitist_constraint_value) )
	{
		elitist_objective_value = solution->objective_value;
		elitist_constraint_value = solution->constraint_value;
	}
}

void fitness_t::partialEvaluationFunction( solution_t *parent, partial_solution_t *solution )
{
	printf("Partial evaluation function not implemented.\n");
	exit(0);
}

bool fitness_t::hasVariableInteractionGraph()
{
	return( variable_interaction_graph.size() > 0 );
}

void fitness_t::initializeVariableInteractionGraph() 
{
}

void fitness_t::initializeRangeBounds()
{
 	lower_range_bound = (double *) Malloc( number_of_parameters*sizeof( double ) );
    upper_range_bound = (double *) Malloc( number_of_parameters*sizeof( double ) );

    for(int i = 0; i < number_of_parameters; i++ )
    {
        lower_range_bound[i] = getLowerRangeBound( i );
        upper_range_bound[i] = getUpperRangeBound( i );
    }
}

sphereFunction_t::sphereFunction_t( int number_of_parameters, double vtr )
{
	this->name = "Sphere function";
	this->number_of_parameters = number_of_parameters;
	this->number_of_subfunctions = number_of_parameters;
	this->vtr = vtr;
	initializeFitnessFunction();
}
		
void sphereFunction_t::evaluationFunction( solution_t *solution )
{
	double result = 0.0;
	for( int i = 0; i < number_of_subfunctions; i++ )
		result += subfunction( solution->variables[i] );
	//solution->buffer = result;

	solution->objective_value = result;
	solution->constraint_value = 0;
	full_number_of_evaluations++;
	number_of_evaluations++;
}

void sphereFunction_t::partialEvaluationFunction( solution_t *parent, partial_solution_t *solution )
{
	double result = 0.0;
	for( int i = 0; i < solution->num_touched_variables; i++ )
	{
		int ind = solution->touched_indices[i];
		result += subfunction( solution->touched_variables[i] );
		result -= subfunction( parent->variables[ind] );
	}
	//solution->buffer = result;
	
	solution->objective_value = parent->objective_value + result;
	solution->constraint_value = parent->constraint_value;
	full_number_of_evaluations++;
	number_of_evaluations += solution->num_touched_variables / (double) number_of_subfunctions;
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

void rosenbrockFunction_t::evaluationFunction( solution_t *solution )
{
	double result = 0.0;
	for( int i = 0; i < number_of_parameters-1; i++ )
		result += subfunction( solution->variables[i], solution->variables[i+1] );
	//solution->buffer = result;

	solution->objective_value = result;
	solution->constraint_value = 0;
	full_number_of_evaluations++;
	number_of_evaluations++;
}

void rosenbrockFunction_t::univariatePartialEvaluationFunction( solution_t *parent, partial_solution_t *solution )
{
	assert( solution->num_touched_variables == 1 );

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

	//solution->buffer = result;
	solution->objective_value = parent->objective_value + result;
	solution->constraint_value = parent->constraint_value;
	full_number_of_evaluations++;
	number_of_evaluations += num_subfunctions_evaluated / (double) number_of_subfunctions;
}


void rosenbrockFunction_t::partialEvaluationFunction( solution_t *parent, partial_solution_t *solution )
{
	if( solution->num_touched_variables == 1 )
	{
		univariatePartialEvaluationFunction( parent, solution );
	}
	else
	{
		int num_subfunctions_evaluated = 0;
		double result = 0.0;
		int *order = mergeSortInt( solution->touched_indices.data(), solution->num_touched_variables );
		for( int i = 0; i < solution->num_touched_variables; i++ )
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
				if( !(i < solution->num_touched_variables-1 && solution->touched_indices[i+1] == ind+1) )
				{
					result += subfunction( solution->touched_variables[i], y );
					result -= subfunction( parent->variables[ind], parent->variables[ind+1] );
					num_subfunctions_evaluated++;
				}
			}
		}
		free( order );

		//solution->buffer = result;
		solution->objective_value = parent->objective_value + result;
		solution->constraint_value = parent->constraint_value;
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

void sorebFunction_t::evaluationFunction( solution_t *solution )
{
	double result = 0.0;
	for( int i = 0; i < number_of_subfunctions; i++ )
		result += subfunction( &solution->variables[getStartingIndexOfBlock(i)], rotation_block_size );

	solution->objective_value = result;
	solution->constraint_value = 0;
	full_number_of_evaluations++;
	number_of_evaluations++;
}

void sorebFunction_t::partialEvaluationFunction( solution_t *parent, partial_solution_t *solution )
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
	printf("fbefore = %10.3e\n",parent->objective_value);

	printf("MODIFIED VARS:\n");
	for( int i = 0; i < solution->num_touched_variables; i++ )
		printf("%d ",solution->touched_indices[i]);
	printf("\n");
	for( int i = 0; i < solution->num_touched_variables; i++ )
		printf("%.3lf ",solution->touched_variables[i]);
	printf("\n");*/
	for( int i = 0; i < solution->num_touched_variables; i++ )
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
				while( i+j < solution->num_touched_variables && solution->touched_indices[i+j]-block_start < rotation_block_size ) // find other touched variables in this block and put them in local array
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
	
	solution->objective_value = parent->objective_value + result;
	solution->constraint_value = parent->constraint_value;
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

void osorebFunction_t::evaluationFunction( solution_t *solution )
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

	solution->objective_value = result;
	solution->constraint_value = 0;
	full_number_of_evaluations++;
	number_of_evaluations++;
}

void osorebFunction_t::partialEvaluationFunction( solution_t *parent, partial_solution_t *solution )
{
	int num_subfunctions_evaluated = 0;
	double result = 0.0;
	for( int i = 0; i < solution->num_touched_variables; i++ )
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
		while( i+j < solution->num_touched_variables && block_ind == solution->touched_indices[i+j] / rotation_block_size )
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
	//solution->buffer = result;
	
	solution->objective_value = parent->objective_value + result;
	solution->constraint_value = parent->constraint_value;
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

void sorebChainFunction_t::evaluationFunction( solution_t *solution )
{
	double result = 0.0;
	for( int i = 0; i < number_of_subfunctions; i++ )
		result += subfunction( &solution->variables[i], rotation_block_size );

	solution->objective_value = result;
	solution->constraint_value = 0;
	full_number_of_evaluations++;
	number_of_evaluations++;
}

void sorebChainFunction_t::partialEvaluationFunction( solution_t *parent, partial_solution_t *solution )
{
	int num_subfunctions_evaluated = 0;
	double result = 0.0;
	for( int i = 0; i < solution->num_touched_variables; i++ )
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
		if( ind < number_of_parameters-1 && !(i < solution->num_touched_variables-1 && solution->touched_indices[i+1] == ind+1) )
		{
			variables_copy[0] = parent->variables[ind];
			variables_copy[1] = parent->variables[ind+1];
			result -= subfunction( variables_copy, rotation_block_size );
			
			variables_copy[0] = solution->touched_variables[i];
			if( i+1 < solution->num_touched_variables && solution->touched_indices[i+1] == ind+1 )
				variables_copy[1] = solution->touched_variables[i+1];
			result += subfunction( variables_copy, rotation_block_size );
			num_subfunctions_evaluated++;
		}
		delete[] variables_copy;
	}
	
	solution->objective_value = parent->objective_value + result;
	solution->constraint_value = parent->constraint_value;
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

void sorebGridFunction_t::evaluationFunction( solution_t *solution )
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

	solution->objective_value = result;
	//printf("\n");
	//exit(0);
	solution->constraint_value = 0;
	full_number_of_evaluations++;
	number_of_evaluations++;
}

void sorebGridFunction_t::partialEvaluationFunction( solution_t *parent, partial_solution_t *solution )
{
	/*double result = 0.0;
	double *variables_copy = new double[2];
	int num_subfunctions_evaluated = 0;
	for( int i = 0; i < solution->num_touched_variables; i++ )
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
	solution->objective_value = parent->objective_value + result;
	solution->constraint_value = parent->constraint_value;
	full_number_of_evaluations++;
	number_of_evaluations += num_subfunctions_evaluated / (double) number_of_subfunctions;*/
	
	std::set<int> subfunction_indices;
	for( int i = 0; i < solution->num_touched_variables; i++ )
	{
		int ind = solution->touched_indices[i];
		subfunction_indices.insert(ind);
		for( int x : getNeighborsInGrid(ind) )
			subfunction_indices.insert(x);
	}

	int num_subfunctions_evaluated;
   	if( solution->num_touched_variables == 1 )
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

void sorebCubeFunction_t::evaluationFunction( solution_t *solution )
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

	solution->objective_value = result;
	solution->constraint_value = 0;
	full_number_of_evaluations++;
	number_of_evaluations++;
}

void sorebCubeFunction_t::partialEvaluationFunction( solution_t *parent, partial_solution_t *solution )
{
	/*int num_subfunctions_evaluated = 0;
	double result = 0.0;
	double *variables_copy = new double[2];
	for( int i = 0; i < solution->num_touched_variables; i++ )
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
	solution->objective_value = parent->objective_value + result;
	solution->constraint_value = parent->constraint_value;
	full_number_of_evaluations++;
	number_of_evaluations += num_subfunctions_evaluated / (double) number_of_subfunctions;*/
	std::set<int> subfunction_indices;
	for( int i = 0; i < solution->num_touched_variables; i++ )
	{
		int ind = solution->touched_indices[i];
		subfunction_indices.insert(ind);
		for( int x : getNeighborsInGrid(ind) )
			subfunction_indices.insert(x);
	}

	int num_subfunctions_evaluated;
   	if( solution->num_touched_variables == 1 )
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

#ifdef CECLSGOFUNC
CECLSGOFunctions_t::CECLSGOFunctions_t( int id, int number_of_parameters, double vtr )
{
	this->name = "CEC Large-scale Global Optimization Benchmarks";
	this->number_of_parameters = number_of_parameters;
	this->vtr = vtr;
	if (id==1){this->function = new F1();}
	else if (id==2){this->function = new F2();}
	else if (id==3){this->function = new F3();}
	else if (id==4){this->function = new F4();}
	else if (id==5){this->function = new F5();}
	else if (id==6){this->function = new F6();}
	else if (id==7){this->function = new F7();}
	else if (id==8){this->function = new F8();}
	else if (id==9){this->function = new F9();}
	else if (id==10){this->function = new F10();}
	else if (id==11){this->function = new F11();}
	else if (id==12){this->function = new F12();}
	else if (id==13){this->function = new F13();}
	else if (id==14){this->function = new F14();}
	else if (id==15){this->function = new F15();}
	initializeFitnessFunction();
}

double CECLSGOFunctions_t::getLowerRangeBound( int dimension )
{
	return( -1e308 );
}
		
double CECLSGOFunctions_t::getUpperRangeBound( int dimension )
{
	return( 1e308 );
}

void CECLSGOFunctions_t::evaluationFunction( solution_t *solution )
{
	solution->objective_value = function->compute((double*)solution->variables.mem);
	/*for( int i = 0; i < number_of_parameters; i++ )
		printf("%10.3e ",solution->variables[i]);
	printf("[ %10.3e ]\n",solution->objective_value);*/
	
	solution->constraint_value = 0;
	full_number_of_evaluations++;
	number_of_evaluations++;
}
#endif

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Problems -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Returns the name of an installed problem.
 */
char *fitness_t::installedProblemName( int index )
{
	if( index > 100 )
		return( (char*) "CEC Large Scale Optimization Benchmarks" );

	if( index >= 20 && index < 100 )
		return( (char*) "Chain of Sum of Rotated Ellipsoid Blocks with variable conditioning number and rotation angle" );
	
	if( index > 1000 )
	{
		return( (char*) "Chain of Sum of Rotated Ellipsoid Blocks with variable block size, overlap size, conditioning number and rotation angle" );
	}
	
    switch( index )
    {
		case  0: return( (char *) "Sphere" );
		case  7: return( (char *) "Rosenbrock" );
		case 13: return( (char *) "Sum of Rotated Ellipsoid Blocks" );
		case 16: return( (char *) "Overlapping Sum of Rotated Ellipsoid Blocks" );
		case 17: return( (char *) "BD2 function hypervolume" );
		case 10: return( (char *) "Chain of Rotated Ellipsoid Blocks" );
		case 20: return( (char *) "Grid of Rotated Ellipsoid Blocks" );
		case 21: return( (char *) "Band of Rotated Ellipsoid Blocks (wrap around x)" );
		case 22: return( (char *) "Torus of Rotated Ellipsoid Blocks (wrap around x and y)" );
		case 30: return( (char *) "Cube of Rotated Ellipsoid Blocks" );
		case 31: return( (char *) "Cube of Rotated Ellipsoid Blocks (wrap around x)" );
		case 32: return( (char *) "Cube of Rotated Ellipsoid Blocks (wrap around x and y)" );
		case 33: return( (char *) "Cube of Rotated Ellipsoid Blocks (wrap around x, y and z)" );
	}

    return( NULL );
}

/**
 * Returns the number of problems installed.
 */
int fitness_t::numberOfInstalledProblems( void )
{
    static int result = -1;

    if( result == -1 )
    {
        result = 0;
        while( installedProblemName( result ) != NULL )
            result++;
    }

    return( result );
}

/**
 * Writes the names of all installed problems to the standard output.
 */
void fitness_t::printAllInstalledProblems( void )
{
    int i, n;

    n = numberOfInstalledProblems();
    printf("Installed optimization problems:\n");
    for( i = 0; i < n; i++ )
        printf("%3d: %s\n", i, installedProblemName( i ));

    exit( 0 );
}

/**
 * Returns whether a parameter is inside the range bound of
 * every problem.
 */
short fitness_t::isParameterInRangeBounds( double parameter, int dimension )
{
    if( parameter < getLowerRangeBound( dimension ) ||
		parameter > getUpperRangeBound( dimension ) ||
		isnan( parameter ) )
    {
        return( 0 );
    }

    return( 1 );
}

void fitness_t::sphereFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
    int    i;
    double result;

    result = 0.0;
    for( i = 0; i < number_of_parameters; i++ )
        result += parameters[i]*parameters[i];

    *objective_value  = result;
    *constraint_value = 0;
}

void fitness_t::sphereFunctionPartialProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before )
{
    int    i;
    double result;

    result = objective_value_before;
    for( i = 0; i < number_of_touched_parameters; i++ )
    {
        result += touched_parameters[i]*touched_parameters[i];
        result -= parameters_before[i]*parameters_before[i];
    }

    *objective_value  = result;
    *constraint_value = 0;
}

void fitness_t::ellipsoidFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
    int    i;
    double result;

    result = 0.0;
    for( i = 0; i < number_of_parameters; i++ )
        result += pow( 10.0, 6.0*(((double) (i))/((double) (number_of_parameters-1))) )*parameters[i]*parameters[i];

    *objective_value  = result;
    *constraint_value = 0;
}

void fitness_t::ellipsoidFunctionPartialProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before )
{
    int    i;
    double result;

    result = objective_value_before;
    for( i = 0; i < number_of_touched_parameters; i++ )
    {
        result += pow( 10.0, 6.0*(((double) (touched_parameters_indices[i]))/((double) (number_of_parameters-1))) )*touched_parameters[i]*touched_parameters[i];
        result -= pow( 10.0, 6.0*(((double) (touched_parameters_indices[i]))/((double) (number_of_parameters-1))) )*parameters_before[i]*parameters_before[i];
    }

    *objective_value  = result;
    *constraint_value = 0;
}

double ellipsoidFunctionLowerRangeBound( int dimension )
{
    return( -1e+308 );
}

double ellipsoidFunctionUpperRangeBound( int dimension )
{
    return( 1e+308 );
}

void fitness_t::cigarFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
    int    i;
    double result;

    result = parameters[0]*parameters[0];
    for( i = 1; i < number_of_parameters; i++ )
    {
        result += pow( 10.0, 6.0 )*parameters[i]*parameters[i];
    }

    *objective_value  = result;
    *constraint_value = 0;
}

double cigarFunctionLowerRangeBound( int dimension )
{
    return( -1e+308 );
}

double cigarFunctionUpperRangeBound( int dimension )
{
    return( 1e+308 );
}

void fitness_t::tabletFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
    int    i;
    double result;

    result = pow( 10.0, 6.0 )*parameters[0]*parameters[0];
    for( i = 1; i < number_of_parameters; i++ )
        result += parameters[i]*parameters[i];

    *objective_value  = result;
    *constraint_value = 0;
}

double tabletFunctionLowerRangeBound( int dimension )
{
    return( -1e+308 );
}

double tabletFunctionUpperRangeBound( int dimension )
{
    return( 1e+308 );
}

void fitness_t::cigarTabletFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
    int    i;
    double result;

    result = parameters[0]*parameters[0];
    for( i = 1; i < number_of_parameters-1; i++ )
        result += pow( 10.0, 4.0 )*parameters[i]*parameters[i];
    result += pow( 10.0, 8.0 )*parameters[number_of_parameters-1]*parameters[number_of_parameters-1];

    *objective_value  = result;
    *constraint_value = 0;
}

double cigarTabletFunctionLowerRangeBound( int dimension )
{
    return( -1e+308 );
}

double cigarTabletFunctionUpperRangeBound( int dimension )
{
    return( 1e+308 );
}

void fitness_t::twoAxesFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
    int    i;
    double result;

    result = 0.0;
    for( i = 0; i <= (number_of_parameters/2)-1; i++ )
        result += pow( 10.0, 6.0 )*parameters[i]*parameters[i];
    for( i = (number_of_parameters/2); i < number_of_parameters; i++ )
        result += parameters[i]*parameters[i];

    *objective_value  = result;
    *constraint_value = 0;
}

double twoAxesFunctionLowerRangeBound( int dimension )
{
    return( -1e+308 );
}

double twoAxesFunctionUpperRangeBound( int dimension )
{
    return( 1e+308 );
}

void fitness_t::differentPowersFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
    int    i;
    double result;

    result = 0.0;
    for( i = 0; i < number_of_parameters; i++ )
        result += pow( fabs(parameters[i]), 2.0 + 10.0*(((double) (i))/((double) (number_of_parameters-1))) );

    *objective_value  = result;
    *constraint_value = 0;
}

double differentPowersFunctionLowerRangeBound( int dimension )
{
    return( -1e+308 );
}

double differentPowersFunctionUpperRangeBound( int dimension )
{
    return( 1e+308 );
}

void fitness_t::rosenbrockFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
    int    i;
    double result;

    result = 0.0;
    for( i = 0; i < number_of_parameters-1; i++ )
        result += 100*(parameters[i+1]-parameters[i]*parameters[i])*(parameters[i+1]-parameters[i]*parameters[i]) + (1.0-parameters[i])*(1.0-parameters[i]);

    *objective_value  = result;
    *constraint_value = 0;
}

void fitness_t::rosenbrockFunctionPartialProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before )
{
    int    i, j;
    double result;

    result = objective_value_before;
    if( number_of_touched_parameters == 1 )
    {
        i = touched_parameters_indices[0];
        if( i > 0 )
        {
            result += 100*(parameters[i]-parameters[i-1]*parameters[i-1])*(parameters[i]-parameters[i-1]*parameters[i-1]) + (1.0-parameters[i-1])*(1.0-parameters[i-1]);
            result -= 100*(parameters_before[0]-parameters[i-1]*parameters[i-1])*(parameters_before[0]-parameters[i-1]*parameters[i-1]) + (1.0-parameters[i-1])*(1.0-parameters[i-1]);
        }
        if( i < number_of_parameters-1 )
        {
            result += 100*(parameters[i+1]-parameters[i]*parameters[i])*(parameters[i+1]-parameters[i]*parameters[i]) + (1.0-parameters[i])*(1.0-parameters[i]);
            result -= 100*(parameters[i+1]-parameters_before[0]*parameters_before[0])*(parameters[i+1]-parameters_before[0]*parameters_before[0]) + (1.0-parameters_before[0])*(1.0-parameters_before[0]);
        }
    }
    else if( number_of_touched_parameters == 2 && touched_parameters_indices[1]-touched_parameters_indices[0] == 1 )
    {
        i = touched_parameters_indices[0];
        j = touched_parameters_indices[1];

        result += 100*(parameters[j]-parameters[i]*parameters[i])*(parameters[j]-parameters[i]*parameters[i]) + (1.0-parameters[i])*(1.0-parameters[i]);
        result -= 100*(parameters_before[1]-parameters_before[0]*parameters_before[0])*(parameters_before[1]-parameters_before[0]*parameters_before[0]) + (1.0-parameters_before[0])*(1.0-parameters_before[0]);
        if( i > 0 )
        {
            result += 100*(parameters[i]-parameters[i-1]*parameters[i-1])*(parameters[i]-parameters[i-1]*parameters[i-1]) + (1.0-parameters[i-1])*(1.0-parameters[i-1]);
            result -= 100*(parameters_before[0]-parameters[i-1]*parameters[i-1])*(parameters_before[0]-parameters[i-1]*parameters[i-1]) + (1.0-parameters[i-1])*(1.0-parameters[i-1]);
        }
        if( j < number_of_parameters-1 )
        {
            result += 100*(parameters[j+1]-parameters[j]*parameters[j])*(parameters[j+1]-parameters[j]*parameters[j]) + (1.0-parameters[j])*(1.0-parameters[j]);
            result -= 100*(parameters[j+1]-parameters_before[1]*parameters_before[1])*(parameters[j+1]-parameters_before[1]*parameters_before[1]) + (1.0-parameters_before[1])*(1.0-parameters_before[1]);
        }
    }
    else
        rosenbrockFunctionProblemEvaluation( parameters, &result, constraint_value );

    *objective_value = result;
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

double rosenbrockFunctionLowerRangeBound( int dimension )
{
    return( -1e+308 );
}

double rosenbrockFunctionUpperRangeBound( int dimension )
{
    return( 1e+308 );
}

void fitness_t::parabolicRidgeFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
    int    i;
    double sum, result;

    sum = 0;
    for( i = 1; i < number_of_parameters; i++ )
        sum += parameters[i]*parameters[i];

    result = -parameters[0] + 100.0*sum;

    *objective_value  = result;
    *constraint_value = 0;
}

double parabolicRidgeFunctionLowerRangeBound( int dimension )
{
    return( -1e+308 );
}

double parabolicRidgeFunctionUpperRangeBound( int dimension )
{
    return( 1e+308 );
}

void fitness_t::sharpRidgeFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
    int    i;
    double sum, result;

    sum = 0;
    for( i = 1; i < number_of_parameters; i++ )
        sum += parameters[i]*parameters[i];

    result = -parameters[0] + 100.0*sqrt( sum );

    *objective_value  = result;
    *constraint_value = 0;
}

double sharpRidgeFunctionLowerRangeBound( int dimension )
{
    return( -1e+308 );
}

double sharpRidgeFunctionUpperRangeBound( int dimension )
{
    return( 1e+308 );
}

void fitness_t::griewankFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
    int    i;
    double yi, sum, prod, result;

    sum  = 0;
    prod = 1.0;
    for( i = 0; i < number_of_parameters; i++ )
    {
        yi    = parameters[i] - 100.0;
        sum  += yi*yi;

        yi    = (parameters[i] - 100.0)/sqrt((double) (i+1));
        prod *= cos( yi );
    }

    result = sum/4000.0 - prod + 1.0;

    *objective_value  = result;
    *constraint_value = 0;
}

double griewankFunctionLowerRangeBound( int dimension )
{
    return( -1e+308 );
}

double griewankFunctionUpperRangeBound( int dimension )
{
    return( 1e+308 );
}

void fitness_t::michalewiczFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
    int    i;
    double result;

    result  = 0.0;
    for( i = 0; i < number_of_parameters; i++ )
        result += -sin(parameters[i])*pow(sin(((i+1)*parameters[i]*parameters[i])/PI),20.0);

    *objective_value  = result;
    *constraint_value = 0;
}

void fitness_t::michalewiczFunctionPartialProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before )
{
    int    i;
    double result;

    result  = objective_value_before;
    for( i = 0; i < number_of_touched_parameters; i++ )
    {
        result += -sin(touched_parameters[i])*pow(sin(((touched_parameters_indices[i]+1)*touched_parameters[i]*touched_parameters[i])/PI),20.0);
        result -= -sin(parameters_before[i])*pow(sin(((touched_parameters_indices[i]+1)*parameters_before[i]*parameters_before[i])/PI),20.0);
    }

    *objective_value  = result;
    *constraint_value = 0;
}

double michalewiczFunctionLowerRangeBound( int dimension )
{
    return( 0.0 );
}

double michalewiczFunctionUpperRangeBound( int dimension )
{
    return( PI );
}

void fitness_t::rastriginFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
    int    i;
    double result;

    result  = 10*number_of_parameters;
    for( i = 0; i < number_of_parameters; i++ )
        result += parameters[i]*parameters[i] - 10.0*cos(2.0*PI*parameters[i]);

    *objective_value  = result;
    *constraint_value = 0;
}

void fitness_t::rastriginFunctionPartialProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before )
{
    int    i;
    double result;

    result  = objective_value_before;
    for( i = 0; i < number_of_touched_parameters; i++ )
    {
        result -= parameters_before[i]*parameters_before[i] - 10.0*cos(2.0*PI*parameters_before[i]);
        result += touched_parameters[i]*touched_parameters[i] - 10.0*cos(2.0*PI*touched_parameters[i]);
    }

    *objective_value  = result;
    *constraint_value = 0;
}

double rastriginFunctionLowerRangeBound( int dimension )
{
    return( -1e+308 );
}

double rastriginFunctionUpperRangeBound( int dimension )
{
    return( 1e+308 );
}

void fitness_t::sumOfEllipsoidsFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
    int    i, j;
    double result;

    result = 0.0;
    for( i = 0; i < number_of_parameters; i++ )
    {
        j = i % rotation_block_size;
        result += pow( 10.0, 6.0*(((double) (j))/((double) (rotation_block_size-1))) )*parameters[i]*parameters[i];
    }

    *objective_value  = result;
    *constraint_value = 0;
}

void fitness_t::sumOfEllipsoidsFunctionPartialProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before )
{
    int    i, j;
    double result;

    result = objective_value_before;

    for( i = 0; i < number_of_touched_parameters; i++ )
    {
        j = touched_parameters_indices[i] % rotation_block_size;
        result -= pow( 10.0, 6.0*(((double) (j))/((double) (rotation_block_size-1))) )*parameters_before[i]*parameters_before[i];
        result += pow( 10.0, 6.0*(((double) (j))/((double) (rotation_block_size-1))) )*touched_parameters[i]*touched_parameters[i];
    }

    *objective_value  = result;
    *constraint_value = 0;
}

double sumOfEllipsoidsFunctionLowerRangeBound( int dimension )
{
    return( -1e+308 );
}

double sumOfEllipsoidsFunctionUpperRangeBound( int dimension )
{
    return( 1e+308 );
}

void fitness_t::ciasBRFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
    int    i, j, nc;
    double result, xi0, xi1, xj0, xj1, distance;

    nc = number_of_parameters/2;

    for( i = 0; i < number_of_parameters; i++ )
    {
        if( parameters[i] < 0 )
            parameters[i] = 0;

        if( parameters[i] > 1 )
            parameters[i] = 1;
    }

    result = -1.0;
    for( i = 0; i < nc; i++ )
        for( j = i+1; j < nc; j++ )
        {
            xi0      = parameters[2*i];
            xi1      = parameters[2*i+1];
            xj0      = parameters[2*j];
            xj1      = parameters[2*j+1];
            distance = (xi0-xj0)*(xi0-xj0) + (xi1-xj1)*(xi1-xj1);
            if( result < 0 || distance < result )
                result = distance;
        }
    result = sqrt( result );

    *objective_value  = -result;
    *constraint_value = 0;
}

double ciasBRFunctionLowerRangeBound( int dimension )
{
    return( -1e+308 );
}

double ciasBRFunctionUpperRangeBound( int dimension )
{
    return( 1e+308 );
}

void fitness_t::trapSphereFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value )
{
    int    i, j, m, k, u;
    double result;

    *objective_value = 0;
    *constraint_value = 0;

    k = 5;
    m = (number_of_parameters / 2) / k;
    for( i = 0; i < m; i++ )
    {
        u = 0;
        for( j = 0; j < k; j++ )
            u += parameters[i*k+j]<0.5?0:1;

        if( u == k )
            result = 1.0;
        else
            result = ((double) (k-1-u))/((double) k);
        *objective_value += (1.0-result);
    }
    for( i = number_of_parameters/2; i < number_of_parameters; i++ )
    {
        *objective_value += parameters[i]*parameters[i];
    }
}

double trapSphereFunctionLowerRangeBound( int dimension )
{
    return( -1e+308 );
}

double trapSphereFunctionUpperRangeBound( int dimension )
{
    return( 1e+308 );
}

/**
 * Computes the rotation matrix to be applied to any solution
 * before evaluating it (i.e. turns the evaluation functions
 * into rotated evaluation functions).
 */
double **fitness_t::initializeObjectiveRotationMatrix( double rotation_angle, int rotation_block_size )
{
    int      i, j, index0, index1;
    double **matrix, **product, theta, cos_theta, sin_theta;

    if( rotation_angle == 0.0 )
        return NULL;

    matrix = (double **) Malloc( rotation_block_size*sizeof( double * ) );
    for( i = 0; i < rotation_block_size; i++ )
        matrix[i] = (double *) Malloc( rotation_block_size*sizeof( double ) );

    double **rotation_matrix = (double **) Malloc( rotation_block_size*sizeof( double * ) );
    for( i = 0; i < rotation_block_size; i++ )
        rotation_matrix[i] = (double *) Malloc( rotation_block_size*sizeof( double ) );

    /* Initialize the rotation matrix to the identity matrix */
    for( i = 0; i < rotation_block_size; i++ )
    {
        for( j = 0; j < rotation_block_size; j++ )
            rotation_matrix[i][j] = 0.0;
        rotation_matrix[i][i] = 1.0;
    }

    /* Construct all rotation matrices (quadratic number) and multiply */
    theta     = (rotation_angle/180.0)*PI;
    cos_theta = cos( theta );
    sin_theta = sin( theta );
    for( index0 = 0; index0 < rotation_block_size-1; index0++ )
    {
        for( index1 = index0+1; index1 < rotation_block_size; index1++ )
        {
            for( i = 0; i < rotation_block_size; i++ )
            {
                for( j = 0; j < rotation_block_size; j++ )
                    matrix[i][j] = 0.0;
                matrix[i][i] = 1.0;
            }
            matrix[index0][index0] = cos_theta;
            matrix[index0][index1] = -sin_theta;
            matrix[index1][index0] = sin_theta;
            matrix[index1][index1] = cos_theta;
	
            product = matrixMatrixMultiplication( matrix, rotation_matrix, rotation_block_size, rotation_block_size, rotation_block_size );
            for( i = 0; i < rotation_block_size; i++ )
                for( j = 0; j < rotation_block_size; j++ )
                    rotation_matrix[i][j] = product[i][j];
			
			/*printf("R[%d][%d]\n",index0,index1);
			for( i = 0; i < rotation_block_size; i++ )
			{
				for( j = 0; j < rotation_block_size; j++ )
					printf("%10.3e ",rotation_matrix[i][j]);
				printf("\n");
			}*/

            for( i = 0; i < rotation_block_size; i++ )
                free( product[i] );
            free( product );
        }
    }

	/*for( i = 0; i < rotation_block_size; i++ )
	{
		for( j = 0; j < rotation_block_size; j++ )
			printf("%10.3e ",rotation_matrix[i][j]);
		printf("\n");
	}
	exit(0);*/

	/*printf("[%10.3e %10.3e]\n",cos_theta,-sin_theta);
	printf("[%10.3e %10.3e]\n",sin_theta,cos_theta);*/
	/*double *vars = new double[rotation_block_size];
	for( i = 0; i < rotation_block_size; i++ )
		vars[i] = 0.0;
	vars[0] = 1.0;
	double *rot_vars = rotateVariables( vars, rotation_block_size, rotation_matrix );
	for( j = 0; j < rotation_block_size; j++ )
		printf("%10.3e ",rot_vars[j]);
	printf("\n");
	vars[0] = 0.0;
	vars[1] = 1.0;
	rot_vars = rotateVariables( vars, rotation_block_size, rotation_matrix );
	for( j = 0; j < rotation_block_size; j++ )
		printf("%10.3e ",rot_vars[j]);
	printf("\n");
	vars[1] = 0.0;
	vars[2] = 1.0;
	rot_vars = rotateVariables( vars, rotation_block_size, rotation_matrix );
	for( j = 0; j < rotation_block_size; j++ )
		printf("%10.3e ",rot_vars[j]);
	printf("\n");
	exit(0);*/

    for( i = 0; i < rotation_block_size; i++ )
        free( matrix[i] );
    free( matrix );

	return( rotation_matrix );
}

void fitness_t::ezilaitiniObjectiveRotationMatrix( double **rotation_matrix, double rotation_angle, int rotation_block_size )
{
    int i;

    if( rotation_angle == 0.0 )
        return;

    for( i = 0; i < rotation_block_size; i++ )
        free( rotation_matrix[i] );
    free( rotation_matrix );
}

int fitness_t::evaluationEmbedded()
{
	if (fitness_embedded() < 0) {
        PyErr_Print();
        fprintf(stderr, "Error in Python code, exception was printed.\n");
		gomealib::utils::freePythonEmbedding();
		exit(1);
    }
	return 0;
}
		
BD2FunctionHypervolume_t::BD2FunctionHypervolume_t( int number_of_parameters, double vtr )
{
	this->name = "BD2 function hypervolume";
	this->number_of_parameters = number_of_parameters;
	this->vtr = vtr;
	this->front_size = 10;
	assert( number_of_parameters % front_size == 0 );
	this->subfunction_size = number_of_parameters / front_size;
	//printf("%d %d %d\n",number_of_parameters,front_size,subfunction_size);
	initializeFitnessFunction();
}

void BD2FunctionHypervolume_t::evaluationFunction( solution_t *solution )
{
	double result = 0.0;
	double *objective_values_front_f0 = new double[front_size];
	double *objective_values_front_f1 = new double[front_size];
	//printf("%d %d %d\n",number_of_parameters,front_size,subfunction_size);
	//solution->print();
	for( int i = 0; i < front_size; i++ )
	{
		objective_values_front_f0[i] = subfunction_f0( &(solution->variables[i*subfunction_size]) ); 
		objective_values_front_f1[i] = subfunction_f1( &(solution->variables[i*subfunction_size]) ); 
		//printf("%10.3e\n",solution->variables[i*subfunction_size]);
		//printf("[%d] %10.3e %10.3e\n",i*subfunction_size,objective_values_front_f0[i], objective_values_front_f1[i]);
	}
	result = compute2DUncrowdedHypervolume( objective_values_front_f0, objective_values_front_f1, front_size );

	delete[] objective_values_front_f0;
	delete[] objective_values_front_f1;

	solution->objective_value = result;
	//printf("%10.3e\n",result);
	solution->constraint_value = 0;
	full_number_of_evaluations++;
	number_of_evaluations++;
}
		
double BD2FunctionHypervolume_t::getLowerRangeBound( int dimension )
{
	return( -1e308 );
}

double BD2FunctionHypervolume_t::getUpperRangeBound( int dimension )
{
	return( 1e308 );
}

double BD2FunctionHypervolume_t::subfunction_f0( double *x )
{
	double result = 0.0;
	for( int i = 0; i < subfunction_size; i++ )
		result += x[i] * x[i];
	result /= (double) subfunction_size;
	//printf("result f0 %10.3e\n",result);
	return( result );
}

double BD2FunctionHypervolume_t::subfunction_f1( double *x )
{
	double result = 0.0;
	for( int i = 0; i < subfunction_size-1; i++ )
		result += 100.0*(x[i+1]-x[i]*x[i])*(x[i+1]-x[i]*x[i]) + (1.0-x[i])*(1.0-x[i]);
	result /= (double) (subfunction_size-1);
	//printf("result f1 %10.3e\n",result);
	return( result );
}

}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
