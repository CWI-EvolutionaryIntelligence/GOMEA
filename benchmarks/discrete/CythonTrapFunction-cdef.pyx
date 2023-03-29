# distutils: language=c++
# cython: language_level=3, boundscheck=False, wraparound=False

import cython
from libcpp.vector cimport vector
import numpy as np
cimport numpy as np
from gomea.fitness cimport GBOFitnessFunction
from gomea.fitness cimport pyGBOFitnessFunction_t
from cpython cimport PyObject

cdef class DeceptiveTrapFunction(GBOFitnessFunction):
    # Any members must be assigned in __new__ to make them accessible during instantiation of superclass
    cdef int k
    def __cinit__(self, number_of_variables : cython.int, k : cython.int ):
        assert( number_of_variables % k == 0 )
        self.k = k
        self.number_of_variables = number_of_variables
        self.value_to_reach = number_of_variables
        self.c_inst_discrete = new pyGBOFitnessFunction_t[char](number_of_variables,self.value_to_reach,<PyObject*>self)

    cpdef int number_of_subfunctions( self ):
        return int(self.number_of_variables / self.k)
    
    cpdef vector[int] inputs_to_subfunction( self, int subfunction_index ):
        cdef vector[int] indices = vector[int](self.k)
        cdef int i
        for i in range(0,self.k,1):
            indices[i] = self.k*subfunction_index + i
        return indices

    cpdef double subfunction(self, int subfunction_index, np.ndarray variables ):
        cdef vector[int] trap_variables = self.inputs_to_subfunction(subfunction_index)
        cdef int unitation = 0
        for v in trap_variables:
            unitation = unitation + variables[v]
        if unitation == self.k:
            return unitation
        else:
            return self.k - unitation - 1

