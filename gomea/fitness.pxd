from cpython cimport PyObject
import numpy as np
cimport numpy as np

cdef extern from "gomea/src/fitness/fitness.hpp" namespace "gomea::fitness":
    cdef cppclass fitness_t[T]:
        fitness_t() except +

cdef extern from "gomea/src/fitness/py_fitness.hpp" namespace "gomea::fitness":
    cdef cppclass pyFitnessFunction_t[T](fitness_t[T]):
        pyFitnessFunction_t() except +
        pyFitnessFunction_t(int,PyObject*) except +
        pyFitnessFunction_t(int,double,PyObject*) except +

cdef extern from "gomea/src/fitness/your_fitness_discrete.hpp" namespace "gomea::fitness":
    cdef cppclass yourFitnessFunctionDiscrete(fitness_t[char]):
        yourFitnessFunctionDiscrete() except +
        yourFitnessFunctionDiscrete(int,double) except +

cdef extern from "gomea/src/fitness/your_fitness_realvalued.hpp" namespace "gomea::fitness":
    cdef cppclass yourFitnessFunctionRealValued(fitness_t[double]):
        yourFitnessFunctionRealValued() except +
        yourFitnessFunctionRealValued(int,double) except +

cdef extern from "gomea/src/fitness/benchmarks-rv.hpp" namespace "gomea::fitness":
    cdef cppclass sphereFunction_t(fitness_t[double]):
        sphereFunction_t() except +
        sphereFunction_t(int,double) except +

cdef extern from "gomea/src/fitness/benchmarks-rv.hpp" namespace "gomea::fitness":
    cdef cppclass rosenbrockFunction_t(fitness_t[double]):
        rosenbrockFunction_t() except +
        rosenbrockFunction_t(int,double) except +

cdef extern from "gomea/src/fitness/benchmarks-discrete.hpp" namespace "gomea::fitness":
    cdef cppclass oneMax_t(fitness_t[char]):
        oneMax_t() except +
        oneMax_t(int) except +

cdef class FitnessFunction:
    cdef fitness_t[char] *c_inst_discrete
    cdef fitness_t[double] *c_inst_realvalued

cdef class GBOFitnessFunction(FitnessFunction):
    cpdef double subfunction( self, int subfunction_index, np.ndarray variables ) except +
    cpdef double objective_function( self, int objective_index, np.ndarray fitness_buffers ) except +
    cpdef double constraint_function( self, np.ndarray fitness_buffers ) except +
    cpdef np.ndarray inputs_to_subfunction( self, int ) except +
    cpdef int number_of_subfunctions( self ) except +

    cpdef int number_of_fitness_buffers( self ) except +
    cpdef int fitness_buffer_index_for_subfunction( self, int ) except +
    
    cpdef double similarity_measure( self, size_t, size_t ) except +

cdef class GBOFitnessFunctionDiscrete(GBOFitnessFunction):
    pass

cdef class GBOFitnessFunctionRealValued(GBOFitnessFunction):
    pass

