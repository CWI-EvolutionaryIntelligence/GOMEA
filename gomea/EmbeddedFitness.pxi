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

cdef public double gomea_pyfitness_objective_function(obj, int objective_index, vector[double] &fitness_buffers ):
    fitness_obj = <PythonFitnessFunction?>obj
    
    cdef void *vec_ptr = &fitness_buffers[0]
    cdef np.npy_intp shape[1]
    shape[0] = <np.npy_intp> fitness_buffers.size()
    npvars = np.PyArray_SimpleNewFromData(1, shape, np.NPY_FLOAT64, vec_ptr)
    
    cdef double result = fitness_obj.objective_function(objective_index,npvars)
    return result

cdef public double gomea_pyfitness_constraint_function(obj, vector[double] &fitness_buffers ):
    fitness_obj = <PythonFitnessFunction?>obj
    
    cdef void *vec_ptr = &fitness_buffers[0]
    cdef np.npy_intp shape[1]
    shape[0] = <np.npy_intp> fitness_buffers.size()
    npvars = np.PyArray_SimpleNewFromData(1, shape, np.NPY_FLOAT64, vec_ptr)
    
    cdef double result = fitness_obj.constraint_function(npvars)
    return result 

cdef public vector[int] gomea_pyfitness_inputsToSubfunction(obj, int subfunction_index ):
    fitness_obj = <PythonFitnessFunction?>obj
    cdef np.ndarray indices = np.array(fitness_obj.inputs_to_subfunction(subfunction_index), np.int32)
    cdef vector[int] vec
    for i in indices:
        vec.push_back(i)
    return vec

cdef public int gomea_pyfitness_index_of_fitness_buffer(obj, int subfunction_index):
    fitness_obj = <PythonFitnessFunction?>obj
    cdef int result = fitness_obj.fitness_buffer_index_for_subfunction(subfunction_index)
    return result

cdef public int gomea_pyfitness_number_of_fitness_buffers(obj):
    fitness_obj = <PythonFitnessFunction?>obj
    cdef int n = fitness_obj.number_of_fitness_buffers()
    return n

cdef public int gomea_pyfitness_numberOfSubfunctions(obj):
    fitness_obj = <PythonFitnessFunction?>obj
    cdef int n = fitness_obj.number_of_subfunctions()
    return n

cdef public double gomea_pyfitness_similarity_measure(obj, size_t var_a, size_t var_b):
    fitness_obj = <PythonFitnessFunction?>obj
    cdef double result = fitness_obj.similarity_measure(var_a,var_b)
    return result
