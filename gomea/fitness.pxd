from cpython cimport PyObject
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string
from libcpp.map cimport map
from libc.math cimport INFINITY
from libcpp.memory cimport shared_ptr
import numpy as np
cimport numpy as np

cdef extern from "gomea/src/fitness/fitness.hpp" namespace "gomea::fitness":
    cdef cppclass fitness_t[T]:
        fitness_t() except +
        int getNumberOfVariables() except +
        double getVTR() except +

cdef extern from "gomea/src/fitness/gbo_fitness.hpp" namespace "gomea::fitness":
    cdef cppclass GBOFitnessFunction_t[T](fitness_t[T]):
        GBOFitnessFunction_t() except +

cdef extern from "gomea/src/fitness/bbo_fitness.hpp" namespace "gomea::fitness":
    cdef cppclass BBOFitnessFunction_t[T](fitness_t[T]):
        BBOFitnessFunction_t() except +

cdef extern from "gomea/src/fitness/py_gbo_fitness.hpp" namespace "gomea::fitness":
    cdef cppclass pyGBOFitnessFunction_t[T](GBOFitnessFunction_t[T]):
        pyGBOFitnessFunction_t() except +
        pyGBOFitnessFunction_t(int,PyObject*) except +
        pyGBOFitnessFunction_t(int,double,PyObject*) except +

cdef extern from "gomea/src/fitness/py_bbo_fitness.hpp" namespace "gomea::fitness":
    cdef cppclass pyBBOFitnessFunction_t[T](fitness_t[T]):
        pyBBOFitnessFunction_t() except +
        pyBBOFitnessFunction_t(int,PyObject*) except +
        pyBBOFitnessFunction_t(int,double,PyObject*) except +

cdef extern from "gomea/src/fitness/your_fitness_discrete.hpp" namespace "gomea::fitness":
    cdef cppclass yourFitnessFunctionDiscrete(GBOFitnessFunction_t[char]):
        yourFitnessFunctionDiscrete() except +
        yourFitnessFunctionDiscrete(int,double) except +

cdef extern from "gomea/src/fitness/your_fitness_realvalued.hpp" namespace "gomea::fitness":
    cdef cppclass yourFitnessFunctionRealValued(GBOFitnessFunction_t[double]):
        yourFitnessFunctionRealValued() except +
        yourFitnessFunctionRealValued(int,double) except +

cdef extern from "gomea/src/fitness/benchmarks-rv.hpp" namespace "gomea::fitness":
    cdef cppclass sphereFunction_t(GBOFitnessFunction_t[double]):
        sphereFunction_t() except +
        sphereFunction_t(int,double) except +

    cdef cppclass sphereFunctionBBO_t(BBOFitnessFunction_t[double]):
        sphereFunctionBBO_t() except +
        sphereFunctionBBO_t(int,double) except +

    cdef cppclass rosenbrockFunction_t(GBOFitnessFunction_t[double]):
        rosenbrockFunction_t() except +
        rosenbrockFunction_t(int,double) except +

    cdef cppclass SOREBChainStrong_t(GBOFitnessFunction_t[double]):
        SOREBChainStrong_t() except +
        SOREBChainStrong_t(int,double) except +
    
    cdef cppclass rosenbrockFunctionBBO_t(BBOFitnessFunction_t[double]):
        rosenbrockFunctionBBO_t() except +
        rosenbrockFunctionBBO_t(int,double) except +

    cdef cppclass SOREBChainStrongBBO_t(BBOFitnessFunction_t[double]):
        SOREBChainStrongBBO_t() except +
        SOREBChainStrongBBO_t(int,double) except +

    cdef cppclass circlesInASquareBBO_t(BBOFitnessFunction_t[double]):
        circlesInASquareBBO_t() except +
        circlesInASquareBBO_t(int,double) except +

cdef extern from "gomea/src/fitness/benchmarks-discrete.hpp" namespace "gomea::fitness":
    cdef cppclass oneMax_t(GBOFitnessFunction_t[char]):
        oneMax_t() except +
        oneMax_t(int) except +

    cdef cppclass deceptiveTrap_t(GBOFitnessFunction_t[char]):
        deceptiveTrap_t() except +
        deceptiveTrap_t(int,int) except +

    cdef cppclass deceptiveTrapBBO_t(BBOFitnessFunction_t[char]):
        deceptiveTrapBBO_t() except +
        deceptiveTrapBBO_t(int,int) except +

    cdef cppclass maxCut_t(GBOFitnessFunction_t[char]):
        maxCut_t() except +
        maxCut_t(string,string) except +

    cdef cppclass maxCutBBO_t(BBOFitnessFunction_t[char]):
        maxCutBBO_t() except +
        maxCutBBO_t(string,string) except +

    cdef cppclass NKlandscapes_t(GBOFitnessFunction_t[char]):
        NKlandscapes_t() except +
        NKlandscapes_t(int,int,long long) except +

    cdef cppclass NKlandscapesBBO_t(BBOFitnessFunction_t[char]):
        NKlandscapesBBO_t() except +
        NKlandscapesBBO_t(int,int,long long) except +

cdef class FitnessFunction:
    cdef public int number_of_variables
    cdef public double value_to_reach
    cdef map[pair[int,double],PyObject*] rotation_matrices
    cdef np.ndarray rotation_matrix

    cdef shared_ptr[fitness_t[char]] c_inst_discrete
    cdef fitness_t[double] *c_inst_realvalued
    
    cpdef initialize_rotation_matrix(self, int rotation_block_size, float rotation_angle)
    cpdef rotate_variables(self, np.ndarray variables, float rotation_angle)

cdef class GBOFitnessFunction(FitnessFunction):
    cpdef double subfunction( self, int subfunction_index, np.ndarray variables ) except? INFINITY
    cpdef double objective_function( self, int objective_index, np.ndarray fitness_buffers ) except? INFINITY
    cpdef double constraint_function( self, np.ndarray fitness_buffers ) except? INFINITY
    cpdef vector[int] inputs_to_subfunction( self, int ) except *
    cpdef int number_of_subfunctions( self ) except -1

    cpdef int number_of_fitness_buffers( self ) except -1
    cpdef int fitness_buffer_index_for_subfunction( self, int ) except -1
    
    cpdef double similarity_measure( self, size_t, size_t ) except? INFINITY

cdef class GBOFitnessFunctionDiscrete(GBOFitnessFunction):
    pass

cdef class GBOFitnessFunctionRealValued(GBOFitnessFunction):
    pass

cdef class BBOFitnessFunction(FitnessFunction):
    cpdef double objective_function( self, int objective_index, np.ndarray variables ) except? INFINITY
    cpdef double constraint_function( self, np.ndarray variables ) except? INFINITY

cdef class BBOFitnessFunctionDiscrete(BBOFitnessFunction):
    pass

cdef class BBOFitnessFunctionRealValued(BBOFitnessFunction):
    pass

