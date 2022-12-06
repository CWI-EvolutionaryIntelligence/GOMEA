#include "gomea/src/common/linkage_model.hpp"

namespace gomea{
        
std::string linkage_model_t::getTypeName( linkage::linkage_model_type type )
{
    switch (type)
    {
        case linkage::UNIVARIATE: return "Univariate";
        case linkage::MPM: return "Marginal Product Model";
		case linkage::LINKAGE_TREE: return "Linkage Tree";
		case linkage::CUSTOM_LM: return "Custom Linkage Model";
		case linkage::FROM_FILE: return "Linkage Model read from file";
		case linkage::CONDITIONAL: return "Conditional";
    }
    return "Unknown type";
}

linkage_model_pt linkage_model_t::createLinkageTreeFOSInstance(size_t FOSIndex, size_t numberOfVariables, int similarityMeasure, int maximumFOSSetSize, bool is_static )
{
    switch (FOSIndex)
    {
        case 0: return linkage_model_t::linkage_tree(numberOfVariables,similarityMeasure,false,maximumFOSSetSize,is_static);
        case 1: return linkage_model_t::linkage_tree(numberOfVariables,similarityMeasure,true,maximumFOSSetSize,is_static);
        default: return( NULL );
    }
	return NULL;
}

linkage_model_pt linkage_model_t::createFOSInstance( const linkage_config_t &config, size_t numberOfVariables )
{
	if( config.type != linkage::FROM_FILE )
		assert( numberOfVariables > 0 );
	if( config.type == linkage::MPM )
		assert( config.mpm_block_size < numberOfVariables );
	switch( config.type )
	{
		case linkage::UNIVARIATE: return univariate(numberOfVariables);
		case linkage::MPM: return marginal_product_model(numberOfVariables, config.mpm_block_size);
		case linkage::LINKAGE_TREE: return linkage_tree(numberOfVariables, config.lt_similarity_measure, config.lt_filtered, config.lt_maximum_set_size, config.lt_is_static );
		case linkage::CUSTOM_LM: return custom_fos(numberOfVariables,config.FOS);
		case linkage::FROM_FILE: return from_file(config.filename);
	}
	throw std::runtime_error("Unknown or unsuitable linkage model.\n");
}

linkage_model_t::linkage_model_t( size_t numberOfVariables_, size_t block_size ) : linkage_model_t(numberOfVariables_)
{
	if( block_size == 0 )
		block_size = numberOfVariables_;
	if( block_size == 1 )
	{
		for (size_t i = 0; i < numberOfVariables_; i++)
			addGroup(i);
		type = linkage::UNIVARIATE;
	}
	else
	{
		assert(numberOfVariables_ % block_size == 0);
		for (int i = 0; i < numberOfVariables_ / block_size; i++)
		{
			std::vector<int> group;
			for (size_t j = 0; j < block_size; j++)
				group.push_back(i * block_size + j);
			addGroup(group);
		}
		if( numberOfVariables == block_size )
			type = linkage::FULL;
		else
			type = linkage::MPM;
	}
	is_static = true;
	shuffleFOS();
}

linkage_model_t::linkage_model_t( size_t numberOfVariables_, const vec_t<vec_t<int>> &FOS )
{
	size_t tot_size = 0;
	for( vec_t<int> group : FOS )
	{
		addGroup(group);
		tot_size += group.size();
	}
	type = linkage::CUSTOM_LM;
	is_static = true;
	shuffleFOS();
}

linkage_model_t::linkage_model_t(size_t numberOfVariables_, int similarityMeasure_, bool filtered_, int maximumSetSize_, bool is_static_ ) : linkage_model_t(numberOfVariables_)
{
	similarityMeasure = similarityMeasure_;
	filtered = filtered_;
	if (maximumSetSize_ > 0)
		maximumSetSize = maximumSetSize_;
	else
		maximumSetSize = numberOfVariables;
	is_static = is_static_;

	S_Matrix.resize(numberOfVariables);
	for (size_t i = 0; i < numberOfVariables; ++i)
		S_Matrix[i].resize(numberOfVariables);
	type = linkage::LINKAGE_TREE;
}

// Read a FOS from a file
linkage_model_t::linkage_model_t( std::string filename )
{
	char    c, string[1000];
	int     i, j, k;
	FILE *file = fopen( filename.c_str(), "r" );
	if( file != NULL )
	{
		/* Length */
		k = 0;
		int length = 0;
		c = fgetc(file);
		while ((c != EOF))
		{
			while (c != '\n')
				c = fgetc(file);
			length++;
			c = fgetc(file);
		}

		fclose(file);
		fflush(stdout);
		file = fopen(filename.c_str(), "r");

		for (i = 0; i < length; i++)
		{
			std::vector<int> vec;
			c = fgetc(file);
			j = 0;
			while ((c != '\n') && (c != EOF))
			{
				k = 0;
				while ((c == ' ') || (c == '\n') || (c == '\t'))
					c = fgetc(file);
				while ((c != ' ') && (c != '\n') && (c != '\t'))
				{
					string[k] = (char)c;
					c = fgetc(file);
					k++;
				}
				string[k] = '\0';
				// printf("FOS[%d][%d] = %d\n",i,j,(int) atoi( string ));
				int e = ((int)atoi(string));
				this->numberOfVariables = fmax(this->numberOfVariables, e);
				vec.push_back(e);
				j++;
			}
			addGroup(vec);
		}
		fclose(file);
		shuffleFOS();
		type = linkage::CUSTOM_LM;
	}
	else
	{
		sprintf( string, "Error reading file %s.\n", filename.c_str() );
		throw std::runtime_error(string);
	}
	is_static = true;
}
    
linkage_model_pt linkage_model_t::univariate(size_t numberOfVariables_)
{
	linkage_model_pt new_fos = std::shared_ptr<linkage_model_t>(new linkage_model_t(numberOfVariables_,0));
	return( new_fos );
}

linkage_model_pt linkage_model_t::marginal_product_model( size_t numberOfVariables_, size_t block_size )
{
	linkage_model_pt new_fos = std::shared_ptr<linkage_model_t>(new linkage_model_t(numberOfVariables_,block_size));
	return( new_fos );
}
        
linkage_model_pt linkage_model_t::linkage_tree(size_t numberOfVariables_, int similarityMeasure_, bool filtered_, int maximumSetSize_, bool is_static_ )
{
	linkage_model_pt new_fos = std::shared_ptr<linkage_model_t>(new linkage_model_t(numberOfVariables_, similarityMeasure_, filtered_, maximumSetSize_, is_static_));
	return( new_fos );
}
    
linkage_model_pt linkage_model_t::from_file( std::string filename )
{
	linkage_model_pt new_fos = std::shared_ptr<linkage_model_t>(new linkage_model_t(filename));
	return( new_fos );
}

linkage_model_pt linkage_model_t::custom_fos( size_t numberOfVariables_, const vec_t<vec_t<int>> &FOS )
{
	linkage_model_pt new_fos = std::shared_ptr<linkage_model_t>(new linkage_model_t(numberOfVariables_,FOS));
	return( new_fos );
}

void linkage_model_t::addGroup( int var_index )
{
	std::vector<int> vec;
	vec.push_back(var_index);
	addGroup(vec);
}

void linkage_model_t::addGroup( const std::set<int> &group ) 
{
	std::vector<int> vec;
	for( int x : group )
		vec.push_back(x);
	addGroup(vec);
}

void linkage_model_t::addGroup( vec_t<int> group ) 
{
	std::sort(group.begin(),group.end());
	FOSStructure.push_back(group);
}

void linkage_model_t::writeToFileFOS(std::string folder, int populationIndex, int generation)
{
    std::ofstream outFile;
    outFile.open(folder + "/fos/" + std::to_string(populationIndex) + "_" + std::to_string(generation) + ".txt");
    outFile << "FOS size = " << FOSStructure.size() << std::endl;
	for (size_t i = 0; i < FOSStructure.size(); ++i)
    {
    	outFile << "[" << i << " : " << FOSStructure[i].size() << "]";
        outFile << "[ ";
        for (size_t j = 0; j < FOSStructure[i].size(); ++j)
            outFile << FOSStructure[i][j] << " "; 
        outFile << "]" << std::endl;
    }
    outFile.close();
}

void linkage_model_t::setCountersToZero()
{
    improvementCounters.resize(FOSStructure.size());
    usageCounters.resize(FOSStructure.size());
    for (size_t i = 0; i < FOSStructure.size(); ++i)
    {
        improvementCounters[i] = 0;
        usageCounters[i] = 0;
    }
}

void linkage_model_t::writeFOSStatistics(std::string folder, int populationIndex, int generation)
{
    std::ofstream outFile;
    outFile.open(folder + "/fos/statistics_" + std::to_string(populationIndex) + "_" + std::to_string(generation) + ".txt");
    outFile << "FOS_element improvement_counter usage_counter" << std::endl;
    for (size_t i = 0; i < FOSStructure.size(); ++i)
    {
        for (size_t j = 0; j < FOSStructure[i].size(); ++j)
            outFile << FOSStructure[i][j] << "_"; 
        outFile << " " << improvementCounters[i] << " " << usageCounters[i] << std::endl;
    }
    outFile.close();
}

void linkage_model_t::shuffleFOS()
{
    FOSorder.resize(size());
    std::iota(FOSorder.begin(), FOSorder.end(), 0);   
    std::shuffle(FOSorder.begin(), FOSorder.end(), utils::rng);
}

void linkage_model_t::determineParallelFOSOrder(std::map<int,std::set<int>> VIG )
{
    FOSorder.clear();
	
	parallelFOSGroups.clear();
	if( colors.size() == 0 )
		colors = graphColoring(VIG);

	int numColors = 0;
	for( size_t i = 0; i < colors.size(); i++ )
		numColors = std::max(numColors,colors[i]+1);
	
	vec_t<vec_t<int>> parallelFOSGroups;
	parallelFOSGroups.resize(numColors);
	for( int i = 0; i < numColors; i++ )
		parallelFOSGroups[i].resize(0);
	for( size_t i = 0; i < colors.size(); i++ )
	{
		int color = colors[i];
		parallelFOSGroups[color].push_back(i);
	}
	
	vec_t<int> groupIndices(numColors);
    std::iota(groupIndices.begin(), groupIndices.end(), 0);   
    std::shuffle(groupIndices.begin(), groupIndices.end(), utils::rng);

	for( size_t i = 0; i < groupIndices.size(); i++ )
	{
		int groupInd = groupIndices[i];
		for( int fos_ind : parallelFOSGroups[groupInd] )
			FOSorder.push_back(fos_ind);
	}

	assert( FOSorder.size() == FOSStructure.size() ); 
}

// VIG_i is the list of all variables dependent on i
vec_t<int> linkage_model_t::graphColoring( std::map<int,std::set<int>> &VIG )
{
	// FOSReverseMap[i] is a list containing the indices of FOS elements that contain variable x_i
	vec_t<vec_t<int>> FOSReverseMap;
	FOSReverseMap.resize(numberOfVariables);
	for( size_t i = 0; i < FOSStructure.size(); i++ )
		for( int v : FOSStructure[i] )
			FOSReverseMap[v].push_back(i);

	// FOSVIG is a VIG of the dependencies between FOS elements
    vec_t<std::set<int>> FOSVIG;
	FOSVIG.resize(FOSStructure.size());
	for( size_t i = 0; i < FOSStructure.size(); i++ )
	{
		for( int v : FOSStructure[i] ) // loop over variables in FOS element
		{
			for( int x : FOSReverseMap[v] ) // add FOS elements that include dependent variables
				if( x != (int) i ) FOSVIG[i].insert(x);
			for( int neighbor : VIG[v] ) // find dependent variables
				for( int x : FOSReverseMap[neighbor] ) // add FOS elements that include dependent variables
					if( x != (int) i ) FOSVIG[i].insert(x);
		}
	}

	// Welsh-Powell algorithm for graph coloring
	vec_t<int> colors;
	colors.resize(FOSVIG.size());
	std::fill(colors.begin(), colors.end(), -1);

	// Sort vertices from high to low degree
    vec_t<int> order(FOSVIG.size());
    std::iota(order.begin(), order.end(), 0);
	std::sort(order.begin(), order.end(), [&FOSVIG](int i, int j){return FOSVIG[i].size() > FOSVIG[j].size();});
	/*for( size_t i = 0; i < order.size(); i++ )
		printf("%d ",FOSVIG[order[i]].size());
	printf("\n");*/

	// Perform coloring
	int numColors = 0;
	for( size_t i = 0; i < order.size(); i++ )
	{
		int ind = order[i];
		if(i>0)assert(FOSVIG[order[i-1]].size()>=FOSVIG[ind].size());
	
		if( numColors == 0 )
		{
			colors[ind] = 0;
			numColors++;
			continue;
		}

		int numColorsSeen = 0;
		bool colorsSeen[numColors] = {};
		int availableColor = -1;
		for( int v : FOSVIG[ind] )
		{
			int color = colors[v];
			if( color != -1 && !colorsSeen[color] )
			{
				colorsSeen[color] = true;
				numColorsSeen++;
			}
			// Check if all colors have been exhausted
			if( numColorsSeen == numColors )
			{
				availableColor = numColors;
				numColors++;
				break;
			}
		}
		// No new color has to be introduced; use first available color
		if( availableColor == -1 )
		{
			assert( numColorsSeen != numColors );
			for( int j = 0; j < numColors; j++ )
			{
				if( !colorsSeen[j] )
				{
					availableColor = j;
					break;
				}
			}
		}
		assert( availableColor != -1 );
		colors[ind] = availableColor;
	}
	

	return( colors );
}
    
void linkage_model_t::printFOS()
{
	printf("Linkage model: (sim:%d,static:%d,)\n",similarityMeasure,is_static);
	for( size_t i = 0; i < FOSStructure.size(); i++ )
	{
		printf("[%d]{",i);
		int c = 0;
		for( int v : FOSStructure[i] )
		{
			if( c == FOSStructure.size()-1 )
				printf("%d",v);
			else
				printf("%d,",v);
			c++;
		}
		printf("}\n");
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
    
void linkage_model_t::learnLinkageTreeFOS(vec_t<solution_t<char>*> &population, size_t alphabetSize )
{
	vec_t<vec_t<double>> MI_matrix;

    /* Compute Mutual Information matrix */
	if (similarityMeasure == 0) // MI
	{
		MI_matrix = computeMIMatrix(population, alphabetSize);
	}
	else if (similarityMeasure == 1) // normalized MI
	{
		MI_matrix = computeNMIMatrix(population, alphabetSize);
	}
    
	learnLinkageTreeFOS(MI_matrix,false);
}

void linkage_model_t::learnLinkageTreeFOS( vec_t<vec_t<double>> similarity_matrix, bool include_full_fos_element )
{
	assert( type == linkage::LINKAGE_TREE );

    FOSStructure.clear();
    vec_t<int> mpmFOSMap;
    vec_t<int> mpmFOSMapNew;
    
    /* Initialize MPM to the univariate factorization */
    vec_t <int> order(numberOfVariables);
    std::iota(order.begin(), order.end(), 0);
    std::shuffle(order.begin(), order.end(), utils::rng );

    vec_t< vec_t<int> > mpm(numberOfVariables);
    vec_t< vec_t<int> > mpmNew(numberOfVariables);
    
    for (size_t i = 0; i < numberOfVariables; i++)
    {
        mpm[i].push_back(order[i]);  
    }

    /* Initialize LT to the initial MPM */
    FOSStructure.resize(numberOfVariables);

    //vec_t<int> useFOSElement(FOSStructure.size(), true);

    int FOSsIndex = 0;
    for (size_t i = 0; i < numberOfVariables; i++)
    {
        FOSStructure[i] = mpm[i];
        mpmFOSMap.push_back(i);
        FOSsIndex++;
    }

    for (size_t i = 0; i < numberOfVariables; ++i)
    {
        for(size_t j = 0; j < numberOfVariables; j++ )
            S_Matrix[i][j] = similarity_matrix[mpm[i][0]][mpm[j][0]];

        S_Matrix[i][i] = 0;
    }
	//printf("Initialized similarity matrix. (%.3fs)\n",getTime(t)/1000.0);

    vec_t<size_t> NN_chain;
    NN_chain.resize(numberOfVariables+2);
    size_t NN_chain_length = 0;
    bool done = false;
    while (!done)
	{
		if (NN_chain_length == 0)
		{
			NN_chain[NN_chain_length] = utils::rng() % mpm.size();

			NN_chain_length++;
		}

		while (NN_chain_length < 3)
		{
			NN_chain[NN_chain_length] = determineNearestNeighbour(NN_chain[NN_chain_length-1], mpm);
			NN_chain_length++;
		}

		while (NN_chain[NN_chain_length-3] != NN_chain[NN_chain_length-1])
		{
			NN_chain[NN_chain_length] = determineNearestNeighbour(NN_chain[NN_chain_length-1], mpm );
			if( ((S_Matrix[NN_chain[NN_chain_length-1]][NN_chain[NN_chain_length]] == S_Matrix[NN_chain[NN_chain_length-1]][NN_chain[NN_chain_length-2]])) && (NN_chain[NN_chain_length] != NN_chain[NN_chain_length-2]) )
				NN_chain[NN_chain_length] = NN_chain[NN_chain_length-2];

			NN_chain_length++;
			if (NN_chain_length > numberOfVariables)
				break;
		}

		size_t r0 = NN_chain[NN_chain_length-2];
		size_t r1 = NN_chain[NN_chain_length-1];
		bool skipFOSElement = false;
		if( filtered )
		{
			if ( (similarityMeasure == 0 || similarityMeasure == 1 ) )
			{
				if( S_Matrix[r1][r0] >= 1-(1e-6))
					skipFOSElement = true;
			}
			else
			{
				if( S_Matrix[r1][r0] == 0 )
					skipFOSElement = true;
			}
		}
		
		if( r1 >= mpm.size() || r0 >= mpm.size() || mpm[r0].size()+mpm[r1].size() > maximumSetSize )
		{
			NN_chain_length = 1;
			NN_chain[0] = 0;
			if( maximumSetSize < numberOfVariables )
			{
				done = 1;
				for(int i = 1; i < mpm.size(); i++ )
				{
					if( mpm[i].size() + mpm[NN_chain[0]].size() <= maximumSetSize ) done = 0;
					if( mpm[i].size() < mpm[NN_chain[0]].size() ) NN_chain[0] = i;
				}
				if( done ) break;
			}
			continue;
		}

		if (r0 > r1)
		{
			int rswap = r0;
			r0 = r1;
			r1 = rswap;
		}
		NN_chain_length -= 3;

		if (r1 < mpm.size()) 
		{
			if( (int) (mpm[r0].size() + mpm[r1].size()) > maximumSetSize )
			{
				done = true;
				break;
			}
			vec_t<int> indices(mpm[r0].size() + mpm[r1].size());

			size_t i = 0;
			for (size_t j = 0; j < mpm[r0].size(); j++)
			{
				indices[i] = mpm[r0][j];
				i++;
			}

			for (size_t j = 0; j < mpm[r1].size(); j++)
			{
				indices[i] = mpm[r1][j];
				i++;
			}

			if( !skipFOSElement )
			{
				FOSStructure.push_back(indices);
				FOSsIndex++;
				assert(FOSStructure.size() == FOSsIndex);
			}

			double mul0 = (double)mpm[r0].size() / (double)(mpm[r0].size() + mpm[r1].size());
			double mul1 = (double)mpm[r1].size() / (double)(mpm[r0].size() + mpm[r1].size());
			for (size_t i = 0; i < mpm.size(); i++)
			{
				if ((i != r0) && (i != r1))
				{
					S_Matrix[i][r0] = mul0 * S_Matrix[i][r0] + mul1 * S_Matrix[i][r1];
					S_Matrix[r0][i] = S_Matrix[i][r0];
				}
			}

			mpmNew.resize(mpm.size() - 1);
			mpmFOSMapNew.resize(mpmFOSMap.size()-1);
			for (size_t i = 0; i < mpmNew.size(); i++)
			{
				mpmNew[i] = mpm[i];
				mpmFOSMapNew[i] = mpmFOSMap[i];
			}

			mpmNew[r0] = indices;
			mpmFOSMapNew[r0] = FOSsIndex-1;

			if (r1 < mpm.size() - 1)
			{
				mpmNew[r1] = mpm[mpm.size() - 1];
				mpmFOSMapNew[r1] = mpmFOSMap[mpm.size() - 1];

				for (size_t i = 0; i < r1; i++)
				{
					S_Matrix[i][r1] = S_Matrix[i][mpm.size() - 1];
					S_Matrix[r1][i] = S_Matrix[i][r1];
				}

				for (size_t j = r1 + 1; j < mpmNew.size(); j++)
				{
					S_Matrix[r1][j] = S_Matrix[j][mpm.size() - 1];
					S_Matrix[j][r1] = S_Matrix[r1][j];
				}
			}

			for (i = 0; i < NN_chain_length; i++)
			{
				if (NN_chain[i] == mpm.size() - 1)
				{
					NN_chain[i] = r1;
					break;
				}
			}

			mpm = mpmNew;
			mpmFOSMap = mpmFOSMapNew;

			if( include_full_fos_element )
			{
				if( mpm.size() == 1 )
					done = true;
			}
			else
			{
				if (mpm.size() == 2)
					done = true;
			}
		}
	}

  /*if (filtered)
  {
	  printf("Erasing elements from Linkage Tree in %.3fs.\n",getTime(t));
    for (size_t i = 0; i < useFOSElement.size(); ++i)
    {
        if (!useFOSElement[i])
        {
            // cout << "filtered out:\n";
            // for (int j = 0 ; j < FOSStructure[i].size(); ++j)
            //     cout << FOSStructure[i][j] << " ";
            // cout << std::endl;
            FOSStructure[i].clear();
        }
    }
    FOSStructure.erase( std::remove_if(FOSStructure.begin(), FOSStructure.end(), [](std::vec_t<int> e){return e.empty();} ), FOSStructure.end());
  }*/
}

/**
 * Determines nearest neighbour according to similarity values.
 */
int linkage_model_t::determineNearestNeighbour(size_t index, vec_t<vec_t< int> > &mpm)
{
    size_t result = 0;

    if (result == index)
        result++;

    for (size_t i = 1; i < mpm.size(); i++)
    {
        if (i != index)
        {
			if( (int) mpm[i].size() > maximumSetSize)
			{
				if( mpm[i].size() < mpm[result].size() )
					result = i;
			}
			else
			{
				if( (int) mpm[result].size() > maximumSetSize)
					result = i;
				else if((S_Matrix[index][i] > S_Matrix[index][result]) || ((S_Matrix[index][i] == S_Matrix[index][result]) && (mpm[i].size() < mpm[result].size())))
					result = i;
			}
        }
    }
    return result;
}

vec_t<vec_t<double>> linkage_model_t::computeMIMatrix( vec_t<solution_t<char>*> &population, size_t alphabetSize )
{
    vec_t<vec_t<double>> MI_Matrix;
	MI_Matrix.resize(numberOfVariables);        
    for (size_t i = 0; i < numberOfVariables; ++i)
        MI_Matrix[i].resize(numberOfVariables);         

    size_t factorSize;
    double p;
    
    /* Compute joint entropy matrix */
    for (size_t i = 0; i < numberOfVariables; i++)
    {
        for (size_t j = i + 1; j < numberOfVariables; j++)
        {
            vec_t<size_t> indices{i, j};
            vec_t<double> factorProbabilities;
            estimateParametersForSingleBinaryMarginal(population, alphabetSize, indices, factorSize, factorProbabilities);

            MI_Matrix[i][j] = 0.0;
            for(size_t k = 0; k < factorSize; k++)
            {
                p = factorProbabilities[k];
                if (p > 0)
                    MI_Matrix[i][j] += -p * log2(p);
            }
            MI_Matrix[j][i] = MI_Matrix[i][j];
        }

        vec_t<size_t> indices{i};
        vec_t<double> factorProbabilities;
        estimateParametersForSingleBinaryMarginal(population, alphabetSize, indices, factorSize, factorProbabilities);

        MI_Matrix[i][i] = 0.0;
        for (size_t k = 0; k < factorSize; k++)
        {
            p = factorProbabilities[k];
            if (p > 0)
                MI_Matrix[i][i] += -p * log2(p);
        }

    }

    /* Then transform into mutual information matrix MI(X,Y)=H(X)+H(Y)-H(X,Y) */
    for (size_t i = 0; i < numberOfVariables; i++)
    {
        for (size_t j = i + 1; j < numberOfVariables; j++)
        {
            MI_Matrix[i][j] = MI_Matrix[i][i] + MI_Matrix[j][j] - MI_Matrix[i][j];
            MI_Matrix[j][i] = MI_Matrix[i][j];
        }
    }

	return( MI_Matrix );
}

vec_t<vec_t<double>> linkage_model_t::computeNMIMatrix( vec_t<solution_t<char>*> &population, size_t alphabetSize )
{
    vec_t<vec_t<double>> MI_Matrix;
	MI_Matrix.resize(numberOfVariables);        
    for (size_t i = 0; i < numberOfVariables; ++i)
    {
		MI_Matrix[i].resize(numberOfVariables);
	}
    
	double p;
    
    /* Compute joint entropy matrix */
    for (size_t i = 0; i < numberOfVariables; i++)
    {
        for (size_t j = i + 1; j < numberOfVariables; j++)
        {
            vec_t<double> factorProbabilities_joint;
            vec_t<double> factorProbabilities_i;
            vec_t<double> factorProbabilities_j;
            size_t factorSize_joint, factorSize_i, factorSize_j;

            vec_t<size_t> indices_joint{i, j};
            estimateParametersForSingleBinaryMarginal(population, alphabetSize, indices_joint, factorSize_joint, factorProbabilities_joint);
            
            vec_t<size_t> indices_i{i};
            estimateParametersForSingleBinaryMarginal(population, alphabetSize, indices_i, factorSize_i, factorProbabilities_i);
            
            vec_t<size_t> indices_j{j};
            estimateParametersForSingleBinaryMarginal(population, alphabetSize, indices_j, factorSize_j, factorProbabilities_j);

            MI_Matrix[i][j] = 0.0;
            
            double separate = 0.0, joint = 0.0;

            for(size_t k = 0; k < factorSize_joint; k++)
            {
                p = factorProbabilities_joint[k];
                //cout << i << " " << j << " " << p << std::endl;
                if (p > 0)
                    joint += (-p * log2(p));
            }

            for(size_t k = 0; k < factorSize_i; k++)
            {
                p = factorProbabilities_i[k];
                if (p > 0)
                    separate += (-p * log2(p));
            }

            for(size_t k = 0; k < factorSize_j; k++)
            {
                p = factorProbabilities_j[k];
                if (p > 0)
                    separate += (-p * log2(p));
            }

            MI_Matrix[i][j] = 0.0;
            if (joint)
                MI_Matrix[i][j] = separate / joint - 1;
            MI_Matrix[j][i] = MI_Matrix[i][j];

        }

    }

	return( MI_Matrix );
}

void linkage_model_t::writeMIMatrixToFile(vec_t<vec_t<double>> MI_Matrix, std::string folder, int populationIndex, int generation)
{
    std::ofstream outFile;
    outFile.open(folder + "/fos/MI_" + std::to_string(populationIndex) + "_" + std::to_string(generation) + ".txt");
    for (size_t i = 0; i < numberOfVariables; ++i)
    {
        for (size_t j = 0; j < numberOfVariables; ++j)
        {       
            outFile << MI_Matrix[i][j] << " "; 
        }
        outFile << std::endl;
    }
    outFile.close();
}

/**
 * Estimates the cumulative probability distribution of a
 * single binary marginal.
 */
void linkage_model_t::estimateParametersForSingleBinaryMarginal(vec_t<solution_t<char>*> &population, size_t alphabetSize, vec_t<size_t> &indices, size_t &factorSize, vec_t<double> &result)
{
    size_t numberOfIndices = indices.size();
    factorSize = (int)pow(alphabetSize, numberOfIndices);

    result.resize(factorSize);
    fill(result.begin(), result.end(), 0.0);

    for (size_t i = 0; i < population.size(); i++)
    {
        int index = 0;
        int power = 1;
        for (int j = numberOfIndices-1; j >= 0; j--)
        {
            index += (int)population[i]->variables[indices[j]] * power;
            power *= alphabetSize;
        }

        result[index] += 1.0;
    }

    for (size_t i = 0; i < factorSize; i++)
        result[i] /= (double)population.size();
}

}