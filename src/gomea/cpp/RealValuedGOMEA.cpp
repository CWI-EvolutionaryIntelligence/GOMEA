#include <Python.h>
#include "real_valued_gomea/rv-gomea.hpp"

/**
 * The main function:
 * - interpret parameters on the command line
 * - run the algorithm with the interpreted parameters
 */
int main( int argc, char **argv )
{
	gomea::realvalued::rvg_t rvgomea = gomea::realvalued::rvg_t(argc, argv);

	rvgomea.run();
	
	return( 0 );
}

