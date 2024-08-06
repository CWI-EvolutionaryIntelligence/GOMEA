#include "gomea/src/utils/embed.hpp"

#include <cstdlib>
#include <cstddef>
#include <stdio.h>
#include <cassert>

namespace gomea{
	namespace utils{
#ifndef CPP_STANDALONE
		namespace{
			bool embedding_initialized = false;
			wchar_t *embedded_program = NULL;
		}

		bool embeddingInitialized()
		{
			return embedding_initialized;
		}

		int initializePythonEmbedding( const char *name, PyObject *(*initfunc)(void) )
		{
			if( embedding_initialized )
				return 0;

			PyObject *pmodule;
			embedded_program = Py_DecodeLocale(name, NULL);
			if (embedded_program == NULL) {
				fprintf(stderr, "Fatal error: cannot decode program_name %s\n",name);
				exit(1);
			}

			// Pass program name to the Python interpreter
			Py_SetProgramName(embedded_program);

			// Initialize the Python interpreter.  Required.
			// If this step fails, it will be a fatal error.
			Py_Initialize();

			embedding_initialized = true;

			// Optionally import the module; alternatively,
			// import can be deferred until the embedded script
			// imports it.
			pmodule = PyImport_ImportModule(name);
			if (!pmodule) {
				PyErr_Print();
				fprintf(stderr, "Error: could not import module '%s'\n",name);
				freePythonEmbedding();
				return 1;
			}

			return 0;
		}

		int freePythonEmbedding()
		{
			assert( embedding_initialized );

			PyMem_RawFree(embedded_program);
			Py_Finalize();

			embedding_initialized = false;

			return 0;
		}
#endif
	}
}
