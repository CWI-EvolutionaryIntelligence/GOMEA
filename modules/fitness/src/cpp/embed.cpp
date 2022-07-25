#include "embed.hpp"
#include "RealValuedGOMEA.h"
//#include "pyFitness.h"

#include <cstdlib>
#include <cstddef>
#include <stdio.h>
#include <cassert>

namespace gomealib{
	namespace utils{
		namespace{
			bool embedding_initialized = false;
			wchar_t *embedded_program = NULL;
		}

		bool embeddingInitialized()
		{
			return embedding_initialized;
		}

		int initializePythonEmbedding()
		{
			const char* name = "RealValuedGOMEA";
			//const char* name = "pyFitness";
			PyObject *pmodule;
			embedded_program = Py_DecodeLocale(name, NULL);
			if (embedded_program == NULL) {
				fprintf(stderr, "Fatal error: cannot decode program_name %s\n",name);
				exit(1);
			}

			// Add a built-in module, before Py_Initialize
			if (PyImport_AppendInittab(name, PyInit_RealValuedGOMEA) == -1) {
				fprintf(stderr, "Error: could not extend in-built modules table\n");
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
	}
}
