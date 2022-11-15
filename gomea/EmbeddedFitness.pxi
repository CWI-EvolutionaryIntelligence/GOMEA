# distutils: language = c++

import ctypes
import numpy as np
cimport numpy as np
from libcpp.vector cimport vector

np.import_array()

cdef public double gomea_pyfitness_subfunction_realvalued(obj, int subfunction_index, vector[double] &variables ):
    fitness_obj = <PythonFitnessFunction?>obj
    
    cdef void *vec_ptr = &variables[0]
    cdef np.npy_intp shape[1]
    shape[0] = <np.npy_intp> variables.size()
    npvars = np.PyArray_SimpleNewFromData(1, shape, np.NPY_FLOAT64, vec_ptr)
    
    cdef double subf = fitness_obj.subfunction(subfunction_index,npvars)
    return subf

cdef public double gomea_pyfitness_subfunction_discrete(obj, int subfunction_index, vector[char] &variables ):
    fitness_obj = <PythonFitnessFunction?>obj
    
    cdef void *vec_ptr = &variables[0]
    cdef np.npy_intp shape[1]
    shape[0] = <np.npy_intp> variables.size()
    npvars = np.PyArray_SimpleNewFromData(1, shape, np.NPY_BYTE, vec_ptr)
    
    cdef double subf = fitness_obj.subfunction(subfunction_index,npvars)
    return subf

cdef public vector[int] gomea_pyfitness_inputsToSubfunction(obj, int subfunction_index ):
    fitness_obj = <PythonFitnessFunction?>obj
    cdef np.ndarray indices = fitness_obj.inputs_to_subfunction(subfunction_index)
    cdef vector[int] vec
    for i in indices:
        vec.push_back(i)
    return vec

cdef public int gomea_pyfitness_numberOfSubfunctions(obj):
    fitness_obj = <PythonFitnessFunction?>obj
    cdef int n = fitness_obj.number_of_subfunctions()
    return n
