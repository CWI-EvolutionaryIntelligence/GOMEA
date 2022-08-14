from libcpp.vector cimport vector
import numpy as np
cimport numpy as np

cdef extern from "fitness/fitness_basic.hpp" namespace "gomea::fitness":
    cdef cppclass fitness_t:
        fitness_t() except +

cdef extern from "fitness/benchmarks-rv.hpp" namespace "gomea::fitness":
    cdef cppclass sphereFunction_t(fitness_t):
        sphereFunction_t() except +
        sphereFunction_t(int,double) except +

cdef extern from "fitness/benchmarks-rv.hpp" namespace "gomea::fitness":
    cdef cppclass rosenbrockFunction_t(fitness_t):
        rosenbrockFunction_t() except +
        rosenbrockFunction_t(int,double) except +

cdef extern from "fitness/py_fitness.hpp" namespace "gomea::fitness":
    cdef cppclass pyFitnessFunction_t(fitness_t):
        pyFitnessFunction_t() except +
        pyFitnessFunction_t(int,double,PyObject) except +

cdef extern from "fitness/your_fitness.hpp" namespace "gomea::fitness":
    cdef cppclass yourFitnessFunction_t(fitness_t):
        yourFitnessFunction_t() except +
        yourFitnessFunction_t(int,double) except +

cdef class PythonFitnessFunction(FitnessFunction):
    cpdef subfunction( self, int subfunction_index, np.ndarray variables )

cdef class FitnessFunction:
    cdef fitness_t *c_inst

