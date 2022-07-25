#include "Python.h"

namespace gomealib{
	namespace utils{

		bool embeddingInitialized();
		int initializePythonEmbedding( const char *name, PyObject *(*initfunc)(void) );
		//int initializePythonEmbedding( const char *name );
		int freePythonEmbedding();
	}
}

