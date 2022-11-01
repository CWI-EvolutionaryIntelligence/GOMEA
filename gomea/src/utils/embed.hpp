#include "Python.h"

namespace gomea{
	namespace utils{

		bool embeddingInitialized();
		int initializePythonEmbedding( const char *name, PyObject *(*initfunc)(void) );
		int freePythonEmbedding();
	}
}

