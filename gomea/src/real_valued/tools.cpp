/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "gomea/src/real_valued/tools.hpp"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

namespace gomea{
namespace realvalued{

int *getRanks( double *array, int array_size )
{
    int i, *sorted, *ranks;

    sorted = utils::mergeSort( array, array_size );
    ranks = (int *) utils::Malloc( array_size * sizeof( int ) );
    for( i = 0; i < array_size; i++ ) ranks[sorted[i]] = i;

    free( sorted );
    return( ranks );
}

int *getRanksFromSorted( int *sorted, int array_size )
{
    int i, *ranks;

    ranks = (int *) utils::Malloc( array_size * sizeof( int ) );
    for( i = 0; i < array_size; i++ ) ranks[sorted[i]] = i;

    return( ranks );
}

vecE random1DNormalUnitVector( int length )
{
    vecE result = vecE(length);
    std::normal_distribution<double> distribution(0.0, 1.0);
    for( int i = 0; i < length; i++ )
        result(i) = distribution(utils::rng);
    return result;
}

}}
