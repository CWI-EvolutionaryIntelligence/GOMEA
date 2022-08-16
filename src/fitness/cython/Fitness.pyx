# distutils: language = c++

from Fitness cimport INFINITY, FitnessFunction, fitness_t, sphereFunction_t, rosenbrockFunction_t, pyFitnessFunction_t, yourFitnessFunction_t
from cpython cimport PyObject
from numpy.math cimport INFINITY

include "EmbeddedFitness.pxi"

cdef class SphereFunction(FitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : double = 0.0
    ):
        self.c_inst = new sphereFunction_t(number_of_variables,value_to_reach)

cdef class RosenbrockFunction(FitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : double = 0.0
    ):
        self.c_inst = new rosenbrockFunction_t(number_of_variables,value_to_reach)

cdef class YourFitnessFunction(FitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : double = 0.0
    ):
        self.c_inst = new yourFitnessFunction_t(number_of_variables,value_to_reach)

cdef class PythonFitnessFunction(FitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : double = 0.0
    ):
        self.number_of_variables = number_of_variables
        self.c_inst = new pyFitnessFunction_t(number_of_variables,value_to_reach,self)

    cpdef number_of_subfunctions( self ):
        return -1
    
    cpdef inputs_to_subfunction( self, int subfunction_index ):
        return np.ndarray()
    
    cpdef subfunction( self, int subfunction_index, np.ndarray variables ):
        return 1e308

