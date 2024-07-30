#include "gomea/src/utils/rng.hpp"

namespace gomea{
namespace utils{

long long random_seed = 0;
std::mt19937 rng;

vec_t<int> randomPermutation( int size )
{
    vec_t<int> perm(size);
    iota(perm.begin(), perm.end(), 0);
    std::shuffle( perm.begin(), perm.end(), rng );
    return( perm );
}

vecE random1DNormalUnitVector( int length )
{
    static std::normal_distribution<double> distribution(0.0, 1.0);
	vecE result = vecE(length);
    for(int i = 0; i < length; i++)
        result[i] = distribution(rng);
    return result;
}

double randomRealUniform01()
{
    static std::uniform_real_distribution<double> distribution(0.0,1.0);
	return distribution(rng);
}

int randomInt( int max )
{
    std::uniform_int_distribution<int> distribution(0,max);
	return distribution(rng);
}

void initializeRandomNumberGenerator()
{
	utils::random_seed = static_cast<long long>(std::chrono::system_clock::now().time_since_epoch().count());
	rng.seed(utils::random_seed);
}

void initializeRandomNumberGenerator( long long seed )
{
    utils::random_seed = seed;
	rng.seed(utils::random_seed);
}

}}
