#include "FOS.hpp"

bool FOSNameByIndex(size_t FOSIndex, string &FOSName)
{
    switch (FOSIndex)
    {
        case 0: FOSName = "Linkage Tree"; break;
        case 1: FOSName = "Filtered Linkage Tree"; break;
            
        default: return false; break;
    }
    return true;
}

FOS_t createFOSInstance(size_t FOSIndex, size_t numberOfVariables, size_t alphabetSize, int similarityMeasure, int maximumFOSSetSize )
{
    switch (FOSIndex)
    {
        case 0: return make_shared<LTFOS>(numberOfVariables, alphabetSize, similarityMeasure, false, maximumFOSSetSize);
        case 1: return make_shared<LTFOS>(numberOfVariables, alphabetSize, similarityMeasure, true, maximumFOSSetSize);
        default: break;
    }
	return NULL;
}

void FOS::writeToFileFOS(string folder, int populationIndex, int generation)
{
    ofstream outFile;
    outFile.open(folder + "/fos/" + to_string(populationIndex) + "_" + to_string(generation) + ".txt");
    outFile << "FOS size = " << FOSStructure.size() << endl;
	for (size_t i = 0; i < FOSStructure.size(); ++i)
    {
    	outFile << "[" << i << " : " << FOSStructure[i].size() << "]";
        outFile << "[ ";
        for (size_t j = 0; j < FOSStructure[i].size(); ++j)
            outFile << FOSStructure[i][j] << " "; 
        outFile << "]" << endl;
    }
    outFile.close();
}

void FOS::setCountersToZero()
{
    improvementCounters.resize(FOSStructure.size());
    usageCounters.resize(FOSStructure.size());
    for (size_t i = 0; i < FOSStructure.size(); ++i)
    {
        improvementCounters[i] = 0;
        usageCounters[i] = 0;
    }
}

void FOS::writeFOSStatistics(string folder, int populationIndex, int generation)
{
    ofstream outFile;
    outFile.open(folder + "/fos/statistics_" + to_string(populationIndex) + "_" + to_string(generation) + ".txt");
    outFile << "FOS_element improvement_counter usage_counter" << endl;
    for (size_t i = 0; i < FOSStructure.size(); ++i)
    {
        for (size_t j = 0; j < FOSStructure[i].size(); ++j)
            outFile << FOSStructure[i][j] << "_"; 
        outFile << " " << improvementCounters[i] << " " << usageCounters[i] << endl;
    }
    outFile.close();
}

void FOS::shuffleFOS(vector<int> &indices, mt19937 *rng)
{
    indices.resize(FOSSize());
    iota(indices.begin(), indices.end(), 0);   
    shuffle(indices.begin(), indices.end(), *rng);
}

void FOS::determineParallelFOSOrder(vector<int> &indices, vector<vector<int>> VIG, mt19937 *rng)
{
    indices.clear();
	
	parallelFOSGroups.clear();
	if( colors.size() == 0 )
		colors = graphColoring(VIG);

	int numColors = 0;
	for( size_t i = 0; i < colors.size(); i++ )
		numColors = max(numColors,colors[i]+1);
	
	vector<vector<int>> parallelFOSGroups;
	parallelFOSGroups.resize(numColors);
	for( size_t i = 0; i < numColors; i++ )
		parallelFOSGroups[i].resize(0);
	for( size_t i = 0; i < colors.size(); i++ )
	{
		int color = colors[i];
		parallelFOSGroups[color].push_back(i);
	}
	
	vector<int> groupIndices(numColors);
    iota(groupIndices.begin(), groupIndices.end(), 0);   
    shuffle(groupIndices.begin(), groupIndices.end(), *rng);

	for( size_t i = 0; i < groupIndices.size(); i++ )
	{
		int groupInd = groupIndices[i];
		for( int fos_ind : parallelFOSGroups[groupInd] )
			indices.push_back(fos_ind);
	}

	assert( indices.size() == FOSStructure.size() ); 
}

// VIG_i is the list of all variables dependent on i
vector<int> FOS::graphColoring( vector<vector<int>> &VIG )
{
	// FOSReverseMap[i] is a list containing the indices of FOS elements that contain variable x_i
	vector<vector<int>> FOSReverseMap;
	FOSReverseMap.resize(numberOfVariables);
	for( size_t i = 0; i < FOSStructure.size(); i++ )
		for( int v : FOSStructure[i] )
			FOSReverseMap[v].push_back(i);

	// FOSVIG is a VIG of the dependencies between FOS elements
    vector<set<int>> FOSVIG;
	FOSVIG.resize(FOSStructure.size());
	for( size_t i = 0; i < FOSStructure.size(); i++ )
	{
		for( int v : FOSStructure[i] ) // loop over variables in FOS element
		{
			for( int x : FOSReverseMap[v] ) // add FOS elements that include dependent variables
				if( x != i ) FOSVIG[i].insert(x);
			for( int neighbor : VIG[v] ) // find dependent variables
				for( int x : FOSReverseMap[neighbor] ) // add FOS elements that include dependent variables
					if( x != i ) FOSVIG[i].insert(x);
		}
	}

	// Welsh-Powell algorithm for graph coloring
	vector<int> colors;
	colors.resize(FOSVIG.size());
	fill(colors.begin(), colors.end(), -1);

	// Sort vertices from high to low degree
    vector<int> order(FOSVIG.size());
    iota(order.begin(), order.end(), 0);
	sort(order.begin(), order.end(), [&FOSVIG](int i, int j){return FOSVIG[i].size() > FOSVIG[j].size();});
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
		//for( size_t j = 0; j < numColors; j++ )
			//colorsSeen[j] = false;
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
			for( size_t j = 0; j < numColors; j++ )
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
	
	/*for( size_t i = 0; i < FOSStructure.size(); i++ )
	{
		printf("[%d][%d][",i,colors[i]);
		for( int v : FOSStructure[i] )
			printf("%d,",v);
		printf("]: {");
		for( int x : FOSVIG[i] )
			printf("%d,",x);
		printf("}\n");
	}*/

	return( colors );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

LTFOS::LTFOS(size_t numberOfVariables_, size_t alphabetSize_, int similarityMeasure_, bool filtered_, int maximumSetSize_): FOS(numberOfVariables_, alphabetSize_)
{
    similarityMeasure = similarityMeasure_;
    filtered = filtered_;
	if( maximumSetSize_ > 0 )
		maximumSetSize = maximumSetSize_;
	else	
		maximumSetSize = numberOfVariables;
    
	S_Matrix.resize(numberOfVariables);
    for (size_t i = 0; i < numberOfVariables; ++i)
        S_Matrix[i].resize(numberOfVariables);
}

void LTFOS::learnFOS(vector<Individual*> &population, mt19937 *rng )
{
	//printf("Initializing MI Matrix.\n");
	long long t = getTimestamp();
	vector<vector<double>> MI_matrix;

    /* Compute Mutual Information matrix */
    if (similarityMeasure == 0) // MI
        MI_matrix = computeMIMatrix(population);
    else if (similarityMeasure == 1) // normalized MI
        MI_matrix = computeNMIMatrix(population);
	//printf("Initialized MI matrix (%.3fs).\n",getTime(t)/1000.0);
    
	learnFOS(MI_matrix, rng);
}

void LTFOS::learnFOS( vector<vector<double>> MI_Matrix, mt19937 *rng )
{
	//printf("Learning Linkage Tree.\n");
	long long t = getTimestamp();
    FOSStructure.clear();
    vector<int> mpmFOSMap;
    vector<int> mpmFOSMapNew;
    
    /* Initialize MPM to the univariate factorization */
    vector <int> order(numberOfVariables);
    iota(order.begin(), order.end(), 0);
    shuffle(order.begin(), order.end(), *rng);

    vector< vector<int> > mpm(numberOfVariables);
    vector< vector<int> > mpmNew(numberOfVariables);
    
    for (size_t i = 0; i < numberOfVariables; i++)
    {
        mpm[i].push_back(order[i]);  
    }

    /* Initialize LT to the initial MPM */
    //int FOSLength = 2 * numberOfVariables - 2;
    FOSStructure.resize(numberOfVariables);

    //vector<int> useFOSElement(FOSStructure.size(), true);

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
            S_Matrix[i][j] = MI_Matrix[mpm[i][0]][mpm[j][0]];

        S_Matrix[i][i] = 0;
    }
	//printf("Initialized similarity matrix. (%.3fs)\n",getTime(t)/1000.0);

    vector<size_t> NN_chain;
    NN_chain.resize(numberOfVariables+2);
    size_t NN_chain_length = 0;
    bool done = false;
    while (!done)
	{
		if (NN_chain_length == 0)
		{
			NN_chain[NN_chain_length] = (*rng)() % mpm.size();
			//std::cout << NN_chain[NN_chain_length] << " | " << mpm.size() << std::endl;

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
		if (filtered && S_Matrix[NN_chain[NN_chain_length-1]][NN_chain[NN_chain_length-2]] >= 1-(1e-6))
			skipFOSElement = true;

		if (r0 > r1)
		{
			int rswap = r0;
			r0 = r1;
			r1 = rswap;
		}
		NN_chain_length -= 3;

		if (r1 < mpm.size()) 
		{
			if( mpm[r0].size() + mpm[r1].size() > maximumSetSize )
			{
				done = true;
				break;
			}
			vector<int> indices(mpm[r0].size() + mpm[r1].size());

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
				//printf("Added FOS element [%d][",FOSsIndex);
				//for(int j : indices) printf("%d ",j);
				//printf("]\n");
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

			if (mpm.size() == 2)
			{
				done = true;
				break;
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
            // cout << endl;
            FOSStructure[i].clear();
        }
    }
    FOSStructure.erase( std::remove_if(FOSStructure.begin(), FOSStructure.end(), [](std::vector<int> e){return e.empty();} ), FOSStructure.end());
  }*/
  //printf("Finished learning Linkage Tree in %.3fs.\n",getTime(t)/1000.0);
}

/**
 * Determines nearest neighbour according to similarity values.
 */
int LTFOS::determineNearestNeighbour(size_t index, vector<vector< int> > &mpm)
{
    size_t result = 0;

    if (result == index)
        result++;

    for (size_t i = 1; i < mpm.size(); i++)
    {
        if (i != index)
        {
			if( mpm[i].size() > maximumSetSize)
			{
				if( mpm[i].size() < mpm[result].size() )
					result = i;
			}
			else
			{
				if( mpm[result].size() > maximumSetSize)
					result = i;
				else if((S_Matrix[index][i] > S_Matrix[index][result]) || ((S_Matrix[index][i] == S_Matrix[index][result]) && (mpm[i].size() < mpm[result].size())))
					result = i;
			}
        }
    }
    return result;
}

vector<vector<double>> LTFOS::computeMIMatrix(vector<Individual*> &population)
{
    vector<vector<double>> MI_Matrix;
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
            vector<size_t> indices{i, j};
            vector<double> factorProbabilities;
            estimateParametersForSingleBinaryMarginal(population, indices, factorSize, factorProbabilities);

            MI_Matrix[i][j] = 0.0;
            for(size_t k = 0; k < factorSize; k++)
            {
                p = factorProbabilities[k];
                if (p > 0)
                    MI_Matrix[i][j] += -p * log2(p);
            }
            MI_Matrix[j][i] = MI_Matrix[i][j];
        }

        vector<size_t> indices{i};
        vector<double> factorProbabilities;
        estimateParametersForSingleBinaryMarginal(population, indices, factorSize, factorProbabilities);

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

vector<vector<double>> LTFOS::computeNMIMatrix(vector<Individual*> &population)
{
    vector<vector<double>> MI_Matrix;
	MI_Matrix.resize(numberOfVariables);        
    for (size_t i = 0; i < numberOfVariables; ++i)
        MI_Matrix[i].resize(numberOfVariables);         
    
	double p;
    
    /* Compute joint entropy matrix */
    for (size_t i = 0; i < numberOfVariables; i++)
    {
        for (size_t j = i + 1; j < numberOfVariables; j++)
        {
            vector<double> factorProbabilities_joint;
            vector<double> factorProbabilities_i;
            vector<double> factorProbabilities_j;
            size_t factorSize_joint, factorSize_i, factorSize_j;

            vector<size_t> indices_joint{i, j};
            estimateParametersForSingleBinaryMarginal(population, indices_joint, factorSize_joint, factorProbabilities_joint);
            
            vector<size_t> indices_i{i};
            estimateParametersForSingleBinaryMarginal(population, indices_i, factorSize_i, factorProbabilities_i);
            
            vector<size_t> indices_j{j};
            estimateParametersForSingleBinaryMarginal(population, indices_j, factorSize_j, factorProbabilities_j);

            MI_Matrix[i][j] = 0.0;
            
            double separate = 0.0, joint = 0.0;

            for(size_t k = 0; k < factorSize_joint; k++)
            {
                p = factorProbabilities_joint[k];
                //cout << i << " " << j << " " << p << endl;
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

void LTFOS::writeMIMatrixToFile(vector<vector<double>> MI_Matrix, string folder, int populationIndex, int generation)
{
    ofstream outFile;
    outFile.open(folder + "/fos/MI_" + to_string(populationIndex) + "_" + to_string(generation) + ".txt");
    for (size_t i = 0; i < numberOfVariables; ++i)
    {
        for (size_t j = 0; j < numberOfVariables; ++j)
        {       
            outFile << MI_Matrix[i][j] << " "; 
        }
        outFile << endl;
    }
    outFile.close();
}

/**
 * Estimates the cumulative probability distribution of a
 * single binary marginal.
 */
void LTFOS::estimateParametersForSingleBinaryMarginal(vector<Individual*> &population, vector<size_t> &indices, size_t &factorSize, vector<double> &result)
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
            index += (int)population[i]->genotype[indices[j]] * power;
            power *= alphabetSize;
        }

        result[index] += 1.0;
    }

    for (size_t i = 0; i < factorSize; i++)
        result[i] /= (double)population.size();
}

