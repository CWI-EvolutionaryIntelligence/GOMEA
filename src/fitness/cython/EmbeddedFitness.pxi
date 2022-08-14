# distutils: language = c++

import ctypes
import numpy as np
cimport numpy as np
cimport libc.stdlib
from libcpp.vector cimport vector
from cpython cimport PyObject

np.import_array()

cdef public double gomea_pyfitness_subfunction(obj, int subfunction_index, vector[double] &variables ) except -1:
    fitness_obj = <PythonFitnessFunction?>obj
    
    cdef void *vec_ptr = &variables[0]
    cdef np.npy_intp shape[1]
    shape[0] = <np.npy_intp> variables.size()
    npvars = np.PyArray_SimpleNewFromData(1, shape, np.NPY_FLOAT64, vec_ptr)
    
    #cdef double subf = fitness_obj.subfunction(subfunction_index,variables)
    cdef double subf = fitness_obj.subfunction(subfunction_index,npvars)
    
    return subf