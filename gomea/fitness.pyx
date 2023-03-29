from cpython cimport PyObject

include "gomea/EmbeddedFitness.pxi"
    
cdef class FitnessFunction:
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : double = 1e308,
        *args,
        **kwargs
    ):
        self.number_of_variables = number_of_variables
        self.value_to_reach = value_to_reach

cdef class YourFitnessFunctionDiscrete(FitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : double = 1e308
    ):
        self.c_inst_discrete = new yourFitnessFunctionDiscrete(number_of_variables,value_to_reach)
    
    def __dealloc__(self):
        del self.c_inst_discrete

cdef class YourFitnessFunctionRealValued(FitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : double = 0.0
    ):
        self.c_inst_realvalued = new yourFitnessFunctionRealValued(number_of_variables,value_to_reach)
    
    def __dealloc__(self):
        del self.c_inst_realvalued

cdef class GBOFitnessFunction(FitnessFunction):
    cpdef int number_of_subfunctions( self ) except +:
        return -1
    
    cpdef vector[int] inputs_to_subfunction( self, int subfunction_index ) except +:
        return []
    
    cpdef double subfunction( self, int subfunction_index, np.ndarray variables ) except +:
        return 1e308
    
    cpdef double objective_function( self, int objective_index, np.ndarray fitness_buffers ) except +:
        return fitness_buffers[objective_index]
    
    cpdef double constraint_function( self, np.ndarray fitness_buffers ) except +:
        return 0
    
    cpdef int number_of_fitness_buffers( self ) except +:
        return 1

    cpdef int fitness_buffer_index_for_subfunction( self, int subfunction_index ) except +: 
        return 0

    cpdef double similarity_measure( self, size_t var_a, size_t var_b ) except +:
        return -1

cdef class GBOFitnessFunctionDiscrete(GBOFitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : double = 1e308
    ):
        self.c_inst_discrete = new pyGBOFitnessFunction_t[char](number_of_variables,value_to_reach,<PyObject*>self)
    
    def __dealloc__(self):
        del self.c_inst_discrete

cdef class GBOFitnessFunctionRealValued(GBOFitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : double = 0.0
    ):
        self.c_inst_realvalued = new pyGBOFitnessFunction_t[double](number_of_variables,value_to_reach,<PyObject*>self)
    
    def __dealloc__(self):
        del self.c_inst_realvalued

cdef class BBOFitnessFunction(FitnessFunction):
    cpdef double objective_function( self, int objective_index, np.ndarray variables ) except +:
        return 1e308
    
    cpdef double constraint_function( self, np.ndarray variables ) except +:
        return 0

cdef class BBOFitnessFunctionDiscrete(BBOFitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : double = 1e308
    ):
        self.c_inst_discrete = new pyBBOFitnessFunction_t[char](number_of_variables,value_to_reach,<PyObject*>self)
    
    def __dealloc__(self):
        del self.c_inst_discrete

cdef class BBOFitnessFunctionRealValued(BBOFitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : double = 0.0
    ):
        self.c_inst_realvalued = new pyBBOFitnessFunction_t[double](number_of_variables,value_to_reach,<PyObject*>self)
    
    def __dealloc__(self):
        del self.c_inst_realvalued


cdef class SphereFunction(GBOFitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : double = 0.0
    ):
        self.c_inst_realvalued = new sphereFunction_t(number_of_variables,value_to_reach)

cdef class RosenbrockFunction(GBOFitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : double = 0.0
    ):
        self.c_inst_realvalued = new rosenbrockFunction_t(number_of_variables,value_to_reach)

cdef class OneMaxFunction(GBOFitnessFunction):
    def __cinit__(self, 
        number_of_variables : int
    ):
        self.c_inst_discrete = new oneMax_t(number_of_variables)

cdef class DeceptiveTrapFunction(GBOFitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        trap_size : int
    ):
        self.c_inst_discrete = new deceptiveTrap_t(number_of_variables,trap_size)

cdef class DeceptiveTrapFunctionBBO(BBOFitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        trap_size : int
    ):
        self.c_inst_discrete = new deceptiveTrapBBO_t(number_of_variables,trap_size)
