/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/fitness/benchmarks-discrete.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include <fstream>

namespace gomea{
namespace fitness{

using namespace gomea;

maxCut_t::maxCut_t( std::string input_file, std::string vtr_file ) : GBOFitnessFunction_t<char>(number_of_variables)
{
	this->name = "MaxCut function";
	readInputFile(input_file);
	readVTRFile(vtr_file);
	this->initialize();
}

void maxCut_t::readInputFile( std::string input_file )
{
	std::ifstream inFile(input_file, std::ifstream::in);
	if (inFile.fail())
	{
		throw std::runtime_error("Problem instance file does not exist!");
	}
	int N, numEdges;
	inFile >> N >> numEdges;

	this->number_of_variables = N;
	for (int i = 0; i < numEdges; ++i)
	{
		int v1, v2;
		double w;
		inFile >> v1 >> v2 >> w;
		vec_t<int> edge;
		edge.push_back(v1-1);
		edge.push_back(v2-1);
		edges.push_back(edge);
		edge_weights.push_back(w);
	}

	inFile.close();
}

void maxCut_t::readVTRFile( std::string vtr_file )
{
	if( vtr_file.size() > 0 )
	{
		this->use_vtr = true;
		std::ifstream inFile(vtr_file, std::ifstream::in);
		if (inFile.fail())
		{
			throw std::runtime_error("VTR file does not exist!");
		}
		double vtr;
		inFile >> vtr;
		this->vtr = vtr;
		inFile.close();
	}
}

int maxCut_t::getNumberOfSubfunctions() 
{
	return edges.size();
}
		
double maxCut_t::subfunction( int subfunction_index, vec_t<char> &variables )
{
	int edge_index = subfunction_index; 
	vec_t<int> inputs = inputsToSubfunction(subfunction_index);
	if( variables[inputs[0]] == variables[inputs[1]] )
		return 0.0;
	else
		return edge_weights[subfunction_index];
}

vec_t<int> maxCut_t::inputsToSubfunction( int subfunction_index )
{
	int edge_index = subfunction_index;
	return( edges[edge_index] );
}


}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
