//#include "utils/embed.hpp"
#include "gomea/src/fitness/fitness.hpp"
#include "gomea/src/fitness/benchmarks-rv.hpp"

namespace gomea{
namespace fitness{

fitness_t::fitness_t( int number_of_parameters, double vtr )
{
	this->name = "Fitness function (C++)";
	this->number_of_parameters = number_of_parameters;
	this->vtr = vtr;
}

fitness_t *fitness_t::getFitnessClass( int problem_index, int number_of_parameters, double vtr )
{
	switch( problem_index )
	{
		case 0 : return( new sphereFunction_t( number_of_parameters, vtr ) );
		case 7 : return( new rosenbrockFunction_t( number_of_parameters, vtr ) );
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

short fitness_t::betterFitness( solution_t<double> *sol_x, solution_t<double> *sol_y ) 
{
	return( betterFitness( sol_x->getObjectiveValue(), sol_x->getConstraintValue(), sol_y->getObjectiveValue(), sol_y->getConstraintValue() ) );
}

int fitness_t::getNumberOfSubfunctions()
{
	return number_of_parameters;
}

double *fitness_t::rotateVariables( double *variables, int num_variables, double **rotation_matrix )
{
	double *rotated_variables = gomea::utils::matrixVectorMultiplication( rotation_matrix, variables, num_variables, num_variables );
    return( rotated_variables );
}

double *fitness_t::rotateVariablesInBlocks( double *variables, int len, int from, int to, double **rotation_matrix )
{
	assert( len % rotation_block_size == 0 );
    int num_blocks = len / rotation_block_size;
	double *rotated_variables = new double[len];
    for( int i = 0; i < from; i++ )
		rotated_variables[i] = variables[i];
	double *cluster = new double[rotation_block_size];
    for( int i = 0; i < num_blocks; i++ )
    {
        for( int j = 0; j < rotation_block_size; j++ )
            cluster[j] = variables[from + i*rotation_block_size + j];
        double *rotated_cluster = gomea::utils::matrixVectorMultiplication( rotation_matrix, cluster, rotation_block_size, rotation_block_size );
        for( int j = 0; j < rotation_block_size; j++ )
            rotated_variables[from + i*rotation_block_size + j] = rotated_cluster[j];
        delete[] rotated_cluster;
    }
	delete[] cluster;

    for( int i = to+1; i < len; i++ )
		rotated_variables[i] = variables[i];

    return( rotated_variables );
}

void fitness_t::evaluate( solution_t<double> *solution )
{
	
	evaluationFunction( solution );
	
	//int out = evaluationEmbedded();
	//assert( out == 0 );	
    
	if( use_vtr && !vtr_hit_status && solution->getConstraintValue() == 0 && solution->getObjectiveValue() <= vtr  )
	{
		vtr_hit_status = 1;
		elitist_objective_value = solution->getObjectiveValue();
		elitist_constraint_value = solution->getConstraintValue();
	}

	if( !vtr_hit_status && betterFitness(solution->getObjectiveValue(), solution->getConstraintValue(), elitist_objective_value, elitist_constraint_value) )
	{
		elitist_objective_value = solution->getObjectiveValue();
		elitist_constraint_value = solution->getConstraintValue();
	}
}

void fitness_t::evaluatePartialSolutionBlackBox( solution_t<double> *parent, partial_solution_t<double> *solution )
{
	// Make backup of parent
	double *var_backup = new double[solution->getNumberOfTouchedVariables()];
	for( int i = 0; i < solution->getNumberOfTouchedVariables(); i++ )
	{
		var_backup[i] = parent->variables[solution->touched_indices[i]];
		parent->variables[solution->touched_indices[i]] = solution->touched_variables[i];
	}
	double obj_val_backup = parent->getObjectiveValue();
	double cons_val_backup = parent->getConstraintValue();

	evaluationFunction( parent );

	// Insert calculated objective and constraint values into partial solution	
	solution->setObjectiveValue(parent->getObjectiveValue());
	solution->setConstraintValue(parent->getConstraintValue());

	// Restore parent to original state
	parent->setObjectiveValue(obj_val_backup);
	parent->setConstraintValue(cons_val_backup);
	for( int i = 0; i < solution->getNumberOfTouchedVariables(); i++ )
		parent->variables[solution->touched_indices[i]] = var_backup[i];
	delete[] var_backup;
}

void fitness_t::evaluatePartialSolution( solution_t<double> *parent, partial_solution_t<double> *solution )
{
	if( black_box_optimization || solution->getNumberOfTouchedVariables() == number_of_parameters )
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

		if( use_vtr && !vtr_hit_status && solution->getConstraintValue() == 0 && solution->getObjectiveValue() <= vtr  )
		{
			evaluatePartialSolutionBlackBox( parent, solution );
			if( solution->getConstraintValue() == 0 && solution->getObjectiveValue() <= vtr  )
			{
				vtr_hit_status = 1;
				elitist_objective_value = solution->getObjectiveValue();
				elitist_constraint_value = solution->getConstraintValue();
			}
		}
	}
	
	if( !vtr_hit_status && betterFitness(solution->getObjectiveValue(), solution->getConstraintValue(), elitist_objective_value, elitist_constraint_value) )
	{
		elitist_objective_value = solution->getObjectiveValue();
		elitist_constraint_value = solution->getConstraintValue();
	}
}

void fitness_t::partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution )
{
	printf("Partial evaluation function not implemented.\n");
	exit(0);
}

void fitness_t::initialize()
{
	initializeSubfunctionDependencyMap();
	initializeVariableInteractionGraph();
}

void fitness_t::initializeSubfunctionDependencyMap()
{
	for( int i = 0; i < number_of_parameters; i++ )
	{
		subfunction_dependency_map[i] = std::set<int>();
	}
	for( int i = 0; i < getNumberOfSubfunctions(); i++ )
	{
		vec_t<int> dependent_variables = inputsToSubfunction(i);
		for( int j : dependent_variables )
		{
			subfunction_dependency_map[j].insert(i);
		}
	}
}

vec_t<int> fitness_t::inputsToSubfunction( int subfunction_index )
{
	vec_t<int> dependencies;
	int parameter_index = subfunction_index;
	dependencies.push_back(parameter_index);
	return dependencies;
}

bool fitness_t::hasVariableInteractionGraph()
{
	return( variable_interaction_graph.size() > 0 );
}

void fitness_t::initializeVariableInteractionGraph() 
{
	return;
}

double fitness_t::getLowerRangeBound( int dimension )
{
	return( -INFINITY );
}
		
double fitness_t::getUpperRangeBound( int dimension )
{
	return( INFINITY );
}


/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Problems -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
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

/**
 * Computes the rotation matrix to be applied to any solution
 * before evaluating it (i.e. turns the evaluation functions
 * into rotated evaluation functions).
 */
double **fitness_t::initializeObjectiveRotationMatrix( double rotation_angle, int rotation_block_size )
{
    if( rotation_angle == 0.0 )
        return NULL;

    double **matrix = new double*[rotation_block_size];
    for( int i = 0; i < rotation_block_size; i++ )
        matrix[i] = new double[rotation_block_size];

    double **rotation_matrix = new double*[rotation_block_size];
    for( int i = 0; i < rotation_block_size; i++ )
        rotation_matrix[i] = new double[rotation_block_size];

    /* Initialize the rotation matrix to the identity matrix */
    for( int i = 0; i < rotation_block_size; i++ )
    {
        for( int j = 0; j < rotation_block_size; j++ )
            rotation_matrix[i][j] = 0.0;
        rotation_matrix[i][i] = 1.0;
    }

    /* Construct all rotation matrices (quadratic number) and multiply */
    double theta     = (rotation_angle/180.0)*M_PI;
    double cos_theta = cos( theta );
    double sin_theta = sin( theta );
    for( int index0 = 0; index0 < rotation_block_size-1; index0++ )
    {
        for( int index1 = index0+1; index1 < rotation_block_size; index1++ )
        {
            for( int i = 0; i < rotation_block_size; i++ )
            {
                for( int j = 0; j < rotation_block_size; j++ )
                    matrix[i][j] = 0.0;
                matrix[i][i] = 1.0;
            }
            matrix[index0][index0] = cos_theta;
            matrix[index0][index1] = -sin_theta;
            matrix[index1][index0] = sin_theta;
            matrix[index1][index1] = cos_theta;
	
            double **product = gomea::utils::matrixMatrixMultiplication( matrix, rotation_matrix, rotation_block_size, rotation_block_size, rotation_block_size );
            for( int i = 0; i < rotation_block_size; i++ )
                for( int j = 0; j < rotation_block_size; j++ )
                    rotation_matrix[i][j] = product[i][j];

            for( int i = 0; i < rotation_block_size; i++ )
                delete[] product[i];
            delete[] product;
        }
    }

    for( int i = 0; i < rotation_block_size; i++ )
	{
        delete[] matrix[i];
	}
	delete[] matrix;

	return( rotation_matrix );
}

void fitness_t::ezilaitiniObjectiveRotationMatrix( double **rotation_matrix, double rotation_angle, int rotation_block_size )
{
    int i;

    if( rotation_angle == 0.0 )
        return;

    for( i = 0; i < rotation_block_size; i++ )
        delete[] rotation_matrix[i];
    delete[] rotation_matrix;
}

/*int fitness_t::evaluationEmbedded()
{
	if (fitness_embedded() < 0) {
        PyErr_Print();
        fprintf(stderr, "Error in Python code, exception was printed.\n");
		gomea::utils::freePythonEmbedding();
		exit(1);
    }
	return 0;
}*/
		
}}