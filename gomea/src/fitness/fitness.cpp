//#include "utils/embed.hpp"
#include "gomea/src/fitness/fitness.hpp"

namespace gomea{
namespace fitness{

template<>
fitness_t<char>::fitness_t() : fitness_t(-1,MAX) {}
template<>
fitness_t<char>::fitness_t( int number_of_variables ) : fitness_t(number_of_variables,MAX) {}
template<>
fitness_t<char>::fitness_t( int number_of_variables, double vtr ) : fitness_t(number_of_variables,MAX) {
	this->vtr = vtr;
	this->use_vtr = true;
}
template<>
fitness_t<char>::fitness_t( int number_of_variables, int alphabet_size ) : fitness_t(number_of_variables,MAX) {
	this->alphabet_size = alphabet_size;
}
template<>
fitness_t<char>::fitness_t( int number_of_variables, int alphabet_size, double vtr ) : fitness_t(number_of_variables,MAX) {
	this->alphabet_size = alphabet_size;
	this->vtr = vtr;
	this->use_vtr = true;
}

template<>
fitness_t<double>::fitness_t() : fitness_t(-1,MIN) {}
template<>
fitness_t<double>::fitness_t( int number_of_variables ) : fitness_t(number_of_variables,MIN) {}
template<>
fitness_t<double>::fitness_t( int number_of_variables, double vtr ) : fitness_t(number_of_variables,MIN) {
	this->vtr = vtr;
	this->use_vtr = true;
}

template<class T>
fitness_t<T>::fitness_t( int number_of_variables, opt_mode optimization_mode )
	: name("Fitness function"), number_of_variables(number_of_variables), optimization_mode(optimization_mode){}

template<class T>
fitness_t<T>::~fitness_t(){
	variable_interaction_graph.clear();
	subfunction_dependency_map.clear();
}

/**
 * Returns 1 if x is better than y, 0 otherwise.
 * x is not better than y unless:
 * - x and y are both infeasible and x has a smaller sum of constraint violations, or
 * - x is feasible and y is not, or
 * - x and y are both feasible and x has a smaller objective value than y
 */
template<class T>
bool fitness_t<T>::betterFitness( double objective_value_x, double constraint_value_x, double objective_value_y, double constraint_value_y )
{
    bool result = false;

    if( constraint_value_x > 0 ) /* x is infeasible */
    {
        if( constraint_value_y > 0 ) /* Both are infeasible */
        {
            if( constraint_value_x < constraint_value_y )
                result = true;
        }
    }
    else /* x is feasible */
    {
        if( constraint_value_y > 0 ) /* x is feasible and y is not */
            result = true;
        else /* Both are feasible */
        {
            if( optimization_mode == MIN && objective_value_x < objective_value_y )
                result = true;
            else if( optimization_mode == MAX && objective_value_x > objective_value_y )
                result = true;
        }
    }

    return( result );
}

template<class T>
bool fitness_t<T>::betterFitness( solution_t<T> *sol_x, solution_t<T> *sol_y ) 
{
	return( betterFitness( sol_x->getObjectiveValue(), sol_x->getConstraintValue(), sol_y->getObjectiveValue(), sol_y->getConstraintValue() ) );
}

template<class T>
void fitness_t<T>::evaluate( solution_t<T> *solution )
{
	checkTermination();

	solution->initObjectiveValues( number_of_objectives );

	auto t = utils::getTimestamp();
	evaluationFunction( solution );
	utils::addToTimer("eval_time",t);

	if( use_vtr && !vtr_hit_status && betterFitness(solution->getObjectiveValue(), solution->getConstraintValue(), vtr, 0.0)  )
	{
		vtr_hit_status = true;
		elitist_objective_value = solution->getObjectiveValue();
		elitist_constraint_value = solution->getConstraintValue();
		elitist_was_written = false;
	}

	if( !vtr_hit_status && betterFitness(solution->getObjectiveValue(), solution->getConstraintValue(), elitist_objective_value, elitist_constraint_value) )
	{
		elitist_objective_value = solution->getObjectiveValue();
		elitist_constraint_value = solution->getConstraintValue();
		elitist_was_written = false;
	}
}

template<class T>
void fitness_t<T>::evaluatePartialSolutionBlackBox( solution_t<T> *parent, partial_solution_t<T> *solution )
{
	checkTermination();
	
	solution->initObjectiveValues( number_of_objectives );
	
	auto t = utils::getTimestamp();
	// Make backup of parent
	double *var_backup = new double[solution->getNumberOfTouchedVariables()];
	for( int i = 0; i < solution->getNumberOfTouchedVariables(); i++ )
	{
		int ind = solution->touched_indices[i];
		var_backup[i] = parent->variables[ind];
		parent->variables[ind] = solution->touched_variables[i];
	}
	double obj_val_backup = parent->getObjectiveValue();
	double cons_val_backup = parent->getConstraintValue();
	//std::vector<double> buffer_backup = parent->getFitnessBuffers();

	evaluationFunction( parent );

	// Insert calculated objective and constraint values into partial solution	
	solution->setObjectiveValue(parent->getObjectiveValue());
	solution->setConstraintValue(parent->getConstraintValue());
	//solution->setFitnessBuffers(parent->getFitnessBuffers());

	// Restore parent to original state
	parent->setObjectiveValue(obj_val_backup);
	parent->setConstraintValue(cons_val_backup);
	//parent->setFitnessBuffers(buffer_backup);
	for( int i = 0; i < solution->getNumberOfTouchedVariables(); i++ )
		parent->variables[solution->touched_indices[i]] = var_backup[i];
	delete[] var_backup;
	
	utils::addToTimer("eval_time",t);
	if( use_vtr && !vtr_hit_status && betterFitness(solution->getObjectiveValue(), solution->getConstraintValue(), vtr, 0.0)  )
	{
		vtr_hit_status = true;
		elitist_objective_value = solution->getObjectiveValue();
		elitist_constraint_value = solution->getConstraintValue();
		elitist_was_written = false;
	}

	if( !vtr_hit_status && betterFitness(solution->getObjectiveValue(), solution->getConstraintValue(), elitist_objective_value, elitist_constraint_value) )
	{
		elitist_objective_value = solution->getObjectiveValue();
		elitist_constraint_value = solution->getConstraintValue();
		elitist_was_written = false;
	}
}

/*template<class T>
void fitness_t<T>::evaluatePartialSolution( solution_t<T> *parent, partial_solution_t<T> *solution )
{
	std::set<int> dependent_subfunctions;
	for( int ind : solution->touched_indices )
	{
		assert( this->subfunction_dependency_map[ind].size() > 0 );
		dependent_subfunctions.insert(this->subfunction_dependency_map[ind].begin(), this->subfunction_dependency_map[ind].end());
	}
	evaluatePartialSolution( parent, solution, dependent_subfunctions );
}*/

//void fitness_t<T>::evaluatePartialSolution( solution_t<T> *parent, partial_solution_t<T> *solution, const std::set<int> &dependent_subfunctions )
template<class T>
void fitness_t<T>::evaluatePartialSolution( solution_t<T> *parent, partial_solution_t<T> *solution ) 
{
	checkTermination();
	
	solution->initObjectiveValues( number_of_objectives );
	
	if( black_box_optimization )
	{
		evaluatePartialSolutionBlackBox( parent, solution );
	}
	else
	{
		//partialEvaluationFunction( parent, solution, dependent_subfunctions );
		auto t = utils::getTimestamp();
		partialEvaluationFunction( parent, solution );
		utils::addToTimer("eval_time",t);

#ifdef CHECK_PARTIAL_FITNESS
		double fbefore = solution->objective_value;
		evaluatePartialSolutionBlackBox( parent, solution );	
		double fafter = solution->objective_value;
		if( std::abs((fbefore-fafter)/fafter) > 1e-4 )
		{
			printf("fbefore = %10.3e; ",fbefore);
			printf("fafter = %10.3e\n",fafter);
			exit(0);
		}
		solution->objective_value = fbefore;
		number_of_evaluations--;
#endif

		if( use_vtr && !vtr_hit_status && betterFitness(solution->getObjectiveValue(), solution->getConstraintValue(), vtr, 0.0)  )
		{
			evaluatePartialSolutionBlackBox( parent, solution );
			if( betterFitness(solution->getObjectiveValue(), solution->getConstraintValue(), vtr, 0.0)  )
			{
				vtr_hit_status = true;
				elitist_objective_value = solution->getObjectiveValue();
				elitist_constraint_value = solution->getConstraintValue();
				elitist_was_written = false;
			}
		}
	}
	
	if( !vtr_hit_status && betterFitness(solution->getObjectiveValue(), solution->getConstraintValue(), elitist_objective_value, elitist_constraint_value) )
	{
		elitist_objective_value = solution->getObjectiveValue();
		elitist_constraint_value = solution->getConstraintValue();
		elitist_was_written = false;
	}
}

//void fitness_t<T>::partialEvaluationFunction( solution_t<T> *parent, partial_solution_t<T> *solution, const std::set<int> &dependent_subfunctions )
template<class T>
void fitness_t<T>::partialEvaluationFunction( solution_t<T> *parent, partial_solution_t<T> *solution )
{
	printf("Partial evaluation function not implemented.\n");
	exit(0);
}

template<class T>
void fitness_t<T>::initialize()
{
	initializeVariableInteractionGraph();
}

template<class T>
void fitness_t<T>::initializeRun()
{
	number_of_evaluations = 0.0;	// discounted in GBO
	full_number_of_evaluations = 0; // not discounted in GBO
	if( opt_mode::MIN == optimization_mode )
		elitist_objective_value = INFINITY;
	else
		elitist_objective_value = -INFINITY;
	elitist_constraint_value = INFINITY;
	vtr_hit_status = false;
}

template<class T>
void fitness_t<T>::initializeVariableInteractionGraph()
{
	return;
}

template<class T>
bool fitness_t<T>::hasVariableInteractionGraph()
{
	return( this->variable_interaction_graph.size() > 0 );
}

template<class T>
void fitness_t<T>::printVariableInteractionGraph() 
{
	for( int i = 0; i < this->number_of_variables; i++ )
	{
		std::set<int> deps = this->variable_interaction_graph[i];
		printf("[%d]->[",i);
		int c = 0;
		for( int v : deps )
		{
			c++;
			printf("%d",v);
			if( c < deps.size() )
				printf(",");
		}
		printf("]\n");
	}
}


template<class T>
vec_t<vec_t<double>> fitness_t<T>::getSimilarityMatrix( int similarity_measure_index )
{
	if( similarity_matrix.size() == 0 )
	{
		similarity_matrix.resize(number_of_variables);
		for (size_t i = 0; i < number_of_variables; i++)
		{
			similarity_matrix[i].resize(number_of_variables);
			similarity_matrix[i][i] = 1e100;
			for (size_t j = 0; j < i; j++)
			{
				double sim;
				if( similarity_measure_index == 2 )
					sim = getSimilarityMeasure(i, j);
				else
					sim = utils::randomRealUniform01();
				similarity_matrix[i][j] = sim;
				similarity_matrix[j][i] = sim;
			}
		}
	}
	return similarity_matrix;
}

template<class T>
double fitness_t<T>::getSimilarityMeasure( size_t var_a, size_t var_b )
{
	assert(0);
	throw std::runtime_error("Fitness function does not implement getSimilarityMeasure(size_t,size_t).");
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Problems -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
double fitness_generic_t::getLowerRangeBound( int dimension )
{
	return( -INFINITY );
}
		
double fitness_generic_t::getUpperRangeBound( int dimension )
{
	return( INFINITY );
}

/**
 * Returns whether a parameter is inside the range bound of
 * every problem.
 */
bool fitness_generic_t::isParameterInRangeBounds( double parameter, int dimension )
{
    if( parameter < getLowerRangeBound( dimension ) ||
		parameter > getUpperRangeBound( dimension ) ||
		std::isnan( parameter ) )
    {
        return( false );
    }

    return( true );
}

/**
 * Computes the rotation matrix to be applied to any solution
 * before evaluating it (i.e. turns the evaluation functions
 * into rotated evaluation functions).
 */
template<>
double **fitness_t<double>::initializeObjectiveRotationMatrix( double rotation_angle, int rotation_block_size )
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
    double theta     = (rotation_angle/180.0)*MY_PI;
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

template<>
void fitness_t<double>::ezilaitiniObjectiveRotationMatrix( double **rotation_matrix, double rotation_angle, int rotation_block_size )
{
    int i;

    if( rotation_angle == 0.0 )
        return;

    for( i = 0; i < rotation_block_size; i++ )
        delete[] rotation_matrix[i];
    delete[] rotation_matrix;
}

template<>
double *fitness_t<double>::rotateVariables( double *variables, int num_variables, double **rotation_matrix )
{
	double *rotated_variables = gomea::utils::matrixVectorMultiplication( rotation_matrix, variables, num_variables, num_variables );
    return( rotated_variables );
}

template<>
double *fitness_t<double>::rotateVariablesInBlocks( double *variables, int len, int from, int to, double **rotation_matrix )
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

template<class T>
void fitness_t<T>::checkTermination() 
{
	checkEvaluationLimitTerminationCondition();

	checkTimeLimitTerminationCondition();
}

template<class T>
void fitness_t<T>::checkEvaluationLimitTerminationCondition() 
{
	if( maximum_number_of_evaluations > 0 && number_of_evaluations >= maximum_number_of_evaluations )
	{
        throw utils::terminationException("evaluations");
	}
}

template<class T>
void fitness_t<T>::checkTimeLimitTerminationCondition() 
{
	if( maximum_number_of_seconds > 0 && utils::getElapsedTimeSinceStartSeconds(utils::start_time) >= maximum_number_of_seconds )
	{
        throw utils::terminationException("time");
	}
}

template<class T>
int fitness_t<T>::getNumberOfVariables()
{
	return this->number_of_variables;
}

template<class T>
double fitness_t<T>::getVTR()
{
	return this->vtr;
}

template<class T>
int fitness_t<T>::getAlphabetSize()
{
	throw std::runtime_error("Fitness function of this type does not implement getAlphabetSize().");
}

template<>
int fitness_t<char>::getAlphabetSize()
{
	return this->alphabet_size;
}

template class fitness_t<char>;
template class fitness_t<double>;

}}
