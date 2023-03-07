from cpython cimport PyObject

include "gomea/EmbeddedFitness.pxi"

cdef class YourFitnessFunctionDiscrete(FitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : double = 0.0
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

cdef class PythonFitnessFunction(FitnessFunction):
    cpdef int number_of_subfunctions( self ) except +:
        return -1
    
    cpdef np.ndarray inputs_to_subfunction( self, int subfunction_index ) except +:
        return np.ndarray()
    
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

cdef class PythonFitnessFunctionDiscrete(PythonFitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        *args,
        **kwargs
    ):
        self.number_of_variables = number_of_variables
        self.c_inst_discrete = new pyFitnessFunction_t[char](number_of_variables,<PyObject*>self)

    def __dealloc__(self):
        del self.c_inst_discrete

cdef class PythonFitnessFunctionRealValued(PythonFitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : double = 0.0
    ):
        self.number_of_variables = number_of_variables
        self.c_inst_realvalued = new pyFitnessFunction_t[double](number_of_variables,value_to_reach,<PyObject*>self)
    
    def __dealloc__(self):
        del self.c_inst_realvalued

cdef class SphereFunction(FitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : double = 0.0
    ):
        self.c_inst_realvalued = new sphereFunction_t(number_of_variables,value_to_reach)

cdef class RosenbrockFunction(FitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : double = 0.0
    ):
        self.c_inst_realvalued = new rosenbrockFunction_t(number_of_variables,value_to_reach)

cdef class OneMaxFunction(FitnessFunction):
    def __cinit__(self, 
        number_of_variables : int
    ):
        self.c_inst_discrete = new oneMax_t(number_of_variables)
