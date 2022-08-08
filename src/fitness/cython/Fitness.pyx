# distutils: language = c++

from Fitness cimport fitness_t, sphereFunction_t, rosenbrockFunction_t, yourFitnessFunction_t

cdef public int fitness_embedded_pub() except -1:
    print("pyfitness_pub - EMBEDDED TEST")
    return 1

cdef int fitness_embedded() except -1:
    print("pyfitness - EMBEDDED TEST")
    return 1

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
