cdef int fitness_embedded() except -1

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

cdef class FitnessFunction:
    cdef fitness_t *c_inst