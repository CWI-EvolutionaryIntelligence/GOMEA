#pragma once

#include "fitness/fitness_basic.hpp"
#include "common/solution.hpp"
#include "common/partial_solution.hpp"
#include "common/gomea_defs.hpp"
#include "utils/tools.hpp"
//#include "real_valued_gomea/cython/RealValuedGOMEA.h"
//#include "fitness/cython/Fitness.h"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{
namespace fitness{

class sphereFunction_t : public fitness_t 
{
	public:
		sphereFunction_t( int number_of_parameters, double vtr );
		
	private:
		void evaluationFunction( solution_t<double> *solution );
		void partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution );
		double subfunction( double x );
};

class rosenbrockFunction_t : public fitness_t 
{
	public:
		rosenbrockFunction_t( int number_of_parameters, double vtr );

		int getNumberOfSubfunctions();

		void initializeVariableInteractionGraph();
		
	private:
		void evaluationFunction( solution_t<double> *solution );
		void partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution );
		void univariatePartialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution );
		double subfunction( double x, double y );
};

class sorebFunction_t : public fitness_t 
{
	public:
		double **rotation_matrix;
		double conditioning_number, rotation_angle;
		int overlap_size;

		sorebFunction_t( int number_of_parameters, double vtr, double conditioning_number, double rotation_angle, int block_size, int overlap_size );
		~sorebFunction_t();
		
		int getNumberOfSubfunctions();
		
		void initializeVariableInteractionGraph();

	private:
		void evaluationFunction( solution_t<double> *solution );
		void partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution );
		double subfunction( double *vars, int num_vars );

		int getStartingIndexOfBlock( int block_index );
		int getIndexOfFirstBlock( int var );
};

class osorebFunction_t : public fitness_t 
{
	public:
		double **rotation_matrix_big, **rotation_matrix_small;
		int number_of_large_rotated_blocks;
		int number_of_small_rotated_blocks;

		int getNumberOfSubfunctions();
		
		osorebFunction_t( int number_of_parameters, double vtr );
		~osorebFunction_t();
		
	private:
		void evaluationFunction( solution_t<double> *solution );
		void partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution );
		double subfunction( double *vars, int num_vars );
};

class sorebChainFunction_t : public fitness_t 
{
	public:
		double **rotation_matrix;
		bool wrap_around;

		sorebChainFunction_t( int number_of_parameters, double vtr, double conditioning_number, double rotation_angle, bool wrap_around );
		~sorebChainFunction_t();
		
		int getNumberOfSubfunctions();
		
		void initializeVariableInteractionGraph();

	private:
		double conditioning_number;

		void evaluationFunction( solution_t<double> *solution );
		void partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution );
		double subfunction( double *vars, int num_vars );
};

class sorebGridFunction_t : public fitness_t 
{
	public:
		bool wrap_around_x;
	   	bool wrap_around_y;
		int grid_width;
		std::map<int,double**> rotation_matrices;

		sorebGridFunction_t( int number_of_parameters, double vtr, double conditioning_number, double rotation_angle, bool wrap_around_x, bool wrap_around_y );
		~sorebGridFunction_t();
		
		void initializeVariableInteractionGraph();

	private:
		double conditioning_number;

		std::set<int> getNeighborsInGrid( int ind );

		void evaluationFunction( solution_t<double> *solution );
		void partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution );
		double subfunction( double *vars, int num_vars );
};

class sorebCubeFunction_t : public fitness_t 
{
	public:
		double **rotation_matrix;
		bool wrap_around_x;
	   	bool wrap_around_y;
	   	bool wrap_around_z;
		int cube_width;
		std::map<int,double**> rotation_matrices;

		sorebCubeFunction_t( int number_of_parameters, double vtr, double conditioning_number, double rotation_angle, bool wrap_around_x, bool wrap_around_y, bool wrap_around_z );
		~sorebCubeFunction_t();
		
		void initializeVariableInteractionGraph();

	private:
		double conditioning_number;

		std::set<int> getNeighborsInGrid( int ind );

		void evaluationFunction( solution_t<double> *solution );
		void partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution );
		double subfunction( double *vars, int num_vars );
};

class BD2FunctionHypervolume_t: public fitness_t 
{
	public:
		BD2FunctionHypervolume_t( int number_of_parameters, double vtr );
		
		int front_size;
		int subfunction_size;
	private:
		void evaluationFunction( solution_t<double> *solution );
		//void partialEvaluationFunction( solution_t<double> *parent, partial_solution_t<double> *solution );
		double subfunction_f0( double *x );
		double subfunction_f1( double *x );
};


}}
