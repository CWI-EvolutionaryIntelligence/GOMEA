#include "gomea/src/fitness/benchmarks-rv.hpp"

namespace gomea{
namespace fitness{

circlesInASquareBBO_t::circlesInASquareBBO_t( int number_of_variables, double vtr ) : BBOFitnessFunction_t(number_of_variables)
{
	assert( number_of_variables % 2 == 0 );
	assert( number_of_variables > 2 );
	this->name = "Circles in a square";
	this->initialize();
	this->initializeOptima();
	if( number_of_variables/2 < optima.size() )
	{
		this->use_vtr = true;
		this->vtr = -1*(optima[number_of_variables/2]-vtr);
	}
}

void circlesInASquareBBO_t::initializeOptima()
{
	optima.resize(31);
	for( int i = 0; i < optima.size(); i++ )
		optima[i] = 1e308; 
	optima[2] = 1.414213562373;
	optima[3] = 1.035276180410;
	optima[4] = 1.000000000000;
	optima[5] = 0.707106781186;
	optima[6] = 0.600925212577;
	optima[7] = 0.535898384862;
	optima[8] = 0.517638090205;
	optima[9] = 0.500000000000;
	optima[10]= 0.421279543983;
	optima[11]= 0.398207310236;
	optima[12]= 0.388730126323;
	optima[13]= 0.366096007696;
	optima[14]= 0.348915260374;
	optima[15]= 0.341081377402;
	optima[16]= 0.333333333333;
	optima[17]= 0.306153985300;
	optima[18]= 0.300462606288;
	optima[19]= 0.289541991994;
	optima[20]= 0.286611652351;
	optima[21]= 0.271812255359;
	optima[22]= 0.267958401550;
	optima[23]= 0.258819045102;
	optima[24]= 0.254333095030;
	optima[25]= 0.250000000000;
	optima[26]= 0.238734757241;
	optima[27]= 0.235849528301;
	optima[28]= 0.230535493642;
	optima[29]= 0.226882900744;
	optima[30]= 0.224502964531;
}

double circlesInASquareBBO_t::objectiveFunction( int objective_index, vec_t<double> &variables )
{
	double min_dist = 1e308;
	for( int i = 0; i < number_of_variables; i+=2 )
	{
		for( int j = i+2; j < number_of_variables; j+=2 )
		{
			double dist = utils::distanceEuclidean2D(variables[i],variables[i+1],variables[j],variables[j+1]);
			min_dist = fmin(min_dist,dist);
		}
	}
	return -min_dist;
}
		
double circlesInASquareBBO_t::getLowerRangeBound( int dimension )
{
	return 0.0;
}

double circlesInASquareBBO_t::getUpperRangeBound( int dimension )
{
	return 1.0;
}

}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
