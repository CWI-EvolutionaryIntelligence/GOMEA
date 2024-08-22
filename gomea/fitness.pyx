from cpython cimport PyObject
from libcpp.string cimport string
from libcpp.pair cimport pair

include "EmbeddedFitness.pxi"
    
cdef class FitnessFunction:
    def __cinit__(self, 
        *args,
        **kwargs
    ):
        c_inst_discrete = NULL
        c_inst_realvalued = NULL
    
    def __init__(self, 
        *args,
        **kwargs
    ):
        self.rotation_matrix = np.ndarray(0)

    cpdef initialize_rotation_matrix(self, int rotation_block_size, float rotation_angle):
        rotation_matrix : np.ndarray = np.identity(rotation_block_size)
        theta : float = np.radians(rotation_angle)
        cos_theta : float = np.cos(theta)
        sin_theta : float = np.sin(theta)
        for index0 in range(rotation_block_size-1):
            for index1 in range(index0+1,rotation_block_size):
                matrix = np.identity(rotation_block_size)
                matrix[index0][index0] = cos_theta;
                matrix[index0][index1] = -sin_theta;
                matrix[index1][index0] = sin_theta;
                matrix[index1][index1] = cos_theta;
                rotation_matrix = np.matmul(matrix, rotation_matrix)
        self.rotation_matrices[(rotation_block_size,rotation_angle)] = <PyObject*> rotation_matrix
        self.rotation_matrix = rotation_matrix

    cpdef rotate_variables(self, np.ndarray variables, float rotation_angle):
        key = (len(variables),rotation_angle)
        if(self.rotation_matrices.count(key) == 0):
            self.initialize_rotation_matrix(len(variables), rotation_angle)
        rotation_matrix = <np.ndarray?> self.rotation_matrices[(len(variables),rotation_angle)]
        rotated_variables = np.matmul(rotation_matrix, variables)
        return rotated_variables
        

cdef class YourFitnessFunctionDiscrete(FitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : float = 1e308
    ):
        self.number_of_variables = number_of_variables
        self.value_to_reach = value_to_reach
        self.c_inst_discrete = new yourFitnessFunctionDiscrete(number_of_variables,value_to_reach)
    
    def __dealloc__(self):
        del self.c_inst_discrete

cdef class YourFitnessFunctionRealValued(FitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : float = 0.0
    ):
        self.number_of_variables = number_of_variables
        self.value_to_reach = value_to_reach
        self.c_inst_realvalued = new yourFitnessFunctionRealValued(number_of_variables,value_to_reach)
    
    def __dealloc__(self):
        del self.c_inst_realvalued

cdef class GBOFitnessFunction(FitnessFunction):
    cpdef int number_of_subfunctions( self ) except -1:
        return -1
    
    cpdef vector[int] inputs_to_subfunction( self, int subfunction_index ) except *:
        return []
    
    cpdef double subfunction( self, int subfunction_index, np.ndarray variables ) except? INFINITY:
        return INFINITY
    
    cpdef double objective_function( self, int objective_index, np.ndarray fitness_buffers ) except? INFINITY:
        return fitness_buffers[objective_index]
    
    cpdef double constraint_function( self, np.ndarray fitness_buffers ) except? INFINITY:
        return 0
    
    cpdef int number_of_fitness_buffers( self ) except -1:
        return 1

    cpdef int fitness_buffer_index_for_subfunction( self, int subfunction_index ) except -1: 
        return 0

    cpdef double similarity_measure( self, size_t var_a, size_t var_b ) except? INFINITY:
        return -1

cdef class GBOFitnessFunctionDiscrete(GBOFitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : float = 1e308
    ):
        self.number_of_variables = number_of_variables
        self.value_to_reach = value_to_reach
        self.c_inst_discrete = new pyGBOFitnessFunction_t[char](number_of_variables,value_to_reach,<PyObject*>self)
    
    def __dealloc__(self):
        del self.c_inst_discrete

cdef class GBOFitnessFunctionRealValued(GBOFitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : float = 0.0
    ):
        self.number_of_variables = number_of_variables
        self.value_to_reach = value_to_reach
        self.c_inst_realvalued = new pyGBOFitnessFunction_t[double](number_of_variables,value_to_reach,<PyObject*>self)
    
    def __dealloc__(self):
        del self.c_inst_realvalued

cdef class BBOFitnessFunction(FitnessFunction):
    cpdef double objective_function( self, int objective_index, np.ndarray variables ) except? INFINITY:
        return INFINITY
    
    cpdef double constraint_function( self, np.ndarray variables ) except? INFINITY:
        return 0

cdef class BBOFitnessFunctionDiscrete(BBOFitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : float = 1e308
    ):
        self.number_of_variables = number_of_variables
        self.value_to_reach = value_to_reach
        self.c_inst_discrete = new pyBBOFitnessFunction_t[char](number_of_variables,value_to_reach,<PyObject*>self)
    
    def __dealloc__(self):
        del self.c_inst_discrete

cdef class BBOFitnessFunctionRealValued(BBOFitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : float = 0.0
    ):
        self.number_of_variables = number_of_variables
        self.value_to_reach = value_to_reach
        self.c_inst_realvalued = new pyBBOFitnessFunction_t[double](number_of_variables,value_to_reach,<PyObject*>self)
    
    def __dealloc__(self):
        del self.c_inst_realvalued


cdef class SphereFunction(GBOFitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : float = 0.0
    ):
        self.c_inst_realvalued = new sphereFunction_t(number_of_variables,value_to_reach)
        self.number_of_variables = self.c_inst_realvalued.getNumberOfVariables()
        self.value_to_reach = self.c_inst_realvalued.getVTR()

cdef class SphereFunctionBBO(BBOFitnessFunction):
    def __cinit__(self,
        number_of_variables : int,
        value_to_reach : float = 0.0
    ):
        self.c_inst_realvalued = new sphereFunctionBBO_t(number_of_variables,value_to_reach)
        self.number_of_variables = self.c_inst_realvalued.getNumberOfVariables()
        self.value_to_reach = self.c_inst_realvalued.getVTR()

cdef class RosenbrockFunction(GBOFitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : float = 0.0
    ):
        self.c_inst_realvalued = new rosenbrockFunction_t(number_of_variables,value_to_reach)
        self.number_of_variables = self.c_inst_realvalued.getNumberOfVariables()
        self.value_to_reach = self.c_inst_realvalued.getVTR()

cdef class RosenbrockFunctionBBO(BBOFitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : float = 0.0
    ):
        self.c_inst_realvalued = new rosenbrockFunctionBBO_t(number_of_variables,value_to_reach)
        self.number_of_variables = self.c_inst_realvalued.getNumberOfVariables()
        self.value_to_reach = self.c_inst_realvalued.getVTR()

cdef class SOREBChainStrong(GBOFitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : float = 0.0
    ):
        self.c_inst_realvalued = new SOREBChainStrong_t(number_of_variables,value_to_reach)
        self.number_of_variables = self.c_inst_realvalued.getNumberOfVariables()
        self.value_to_reach = self.c_inst_realvalued.getVTR()

cdef class SOREBChainStrongBBO(BBOFitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : float = 0.0
    ):
        self.c_inst_realvalued = new SOREBChainStrongBBO_t(number_of_variables,value_to_reach)
        self.number_of_variables = self.c_inst_realvalued.getNumberOfVariables()
        self.value_to_reach = self.c_inst_realvalued.getVTR()

cdef class CirclesInASquareBBO(BBOFitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        value_to_reach : float = 0.0
    ):
        self.c_inst_realvalued = new circlesInASquareBBO_t(number_of_variables,value_to_reach)
        self.number_of_variables = self.c_inst_realvalued.getNumberOfVariables()
        self.value_to_reach = self.c_inst_realvalued.getVTR()


cdef class OneMaxFunction(GBOFitnessFunction):
    def __cinit__(self, 
        number_of_variables : int
    ):
        self.c_inst_discrete = new oneMax_t(number_of_variables)
        self.number_of_variables = self.c_inst_discrete.getNumberOfVariables()
        self.value_to_reach = self.c_inst_discrete.getVTR()

cdef class DeceptiveTrapFunction(GBOFitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        trap_size : int
    ):
        self.c_inst_discrete = new deceptiveTrap_t(number_of_variables,trap_size)
        self.number_of_variables = self.c_inst_discrete.getNumberOfVariables()
        self.value_to_reach = self.c_inst_discrete.getVTR()

cdef class DeceptiveTrapFunctionBBO(BBOFitnessFunction):
    def __cinit__(self, 
        number_of_variables : int,
        trap_size : int
    ):
        self.c_inst_discrete = new deceptiveTrapBBO_t(number_of_variables,trap_size)
        self.number_of_variables = self.c_inst_discrete.getNumberOfVariables()
        self.value_to_reach = self.c_inst_discrete.getVTR()

cdef class MaxCut(GBOFitnessFunction):
    def __cinit__(self,
        input_file : str,
        vtr_file : str = ""
    ):
        self.c_inst_discrete = new maxCut_t(str.encode(input_file),str.encode(vtr_file))
        self.number_of_variables = self.c_inst_discrete.getNumberOfVariables()
        self.value_to_reach = self.c_inst_discrete.getVTR()

cdef class MaxCutBBO(BBOFitnessFunction):
    def __cinit__(self,
        input_file : str,
        vtr_file : str = ""
    ):
        self.c_inst_discrete = new maxCutBBO_t(str.encode(input_file),str.encode(vtr_file))
        self.number_of_variables = self.c_inst_discrete.getNumberOfVariables()
        self.value_to_reach = self.c_inst_discrete.getVTR()
