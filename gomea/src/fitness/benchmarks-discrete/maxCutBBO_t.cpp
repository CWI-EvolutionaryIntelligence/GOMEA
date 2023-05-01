/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/fitness/benchmarks-discrete.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include <fstream>

namespace gomea{
namespace fitness{

maxCutBBO_t::maxCutBBO_t( std::string input_file, std::string vtr_file ) : BBOFitnessFunction_t<char>(number_of_variables)
{
	this->name = "MaxCut function";
	readInputFile(input_file);
	readVTRFile(vtr_file);
	this->initialize();
}

void maxCutBBO_t::readInputFile( std::string input_file )
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

void maxCutBBO_t::readVTRFile( std::string vtr_file )
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
		
double maxCutBBO_t::objectiveFunction( int objective_index, vec_t<char> &variables )
{
	double cut = 0.0;
	for( int i = 0; i < edges.size(); i++ )
	{
		vec_t<int> edge = edges[i];
		if (variables[edge[0]] != variables[edge[1]])
			cut += edge_weights[i];
	}
	return cut;
}

}}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
