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

#pragma once

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "tools.hpp"
#include "fitness_buffer.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{
namespace realvalued{

template<class T>
class partial_solution_t
{
	public:
		vec_t<int> touched_indices;
		vec_t<T> touched_variables;
		vec_t<T> sample_zs; // Samples z~N(0,I), later transformed to N(mu,C)
		vec_t<T> sample_means;

		double buffer;
		short is_accepted = 0;
		short improves_elitist = 0;

		partial_solution_t( int num_touched_variables );
		partial_solution_t( vec_t<T> &touched_variables, vec_t<int> &touched_indices );
		partial_solution_t( vec_t<T> &touched_variables, vec_t<T> &sample_zs, vec_t<int> &touched_indices );
		partial_solution_t( partial_solution_t<T> &other);

		int getNumberOfTouchedVariables();
		double getObjectiveValue();
		double getObjectiveValue( int objective_value_index );
		double getConstraintValue();

		void setObjectiveValue( double v );
		void setObjectiveValue( double v, int objective_value_index );
		void setConstraintValue( double v );

		int getTouchedIndex( int ind );
		void setSampleMean( vec_t<T> &means );

		void print();
	private:
		double objective_value = 1e308;
		double constraint_value = 1e308;
		
		std::map<int,int> touched_index_map;
};

}}
