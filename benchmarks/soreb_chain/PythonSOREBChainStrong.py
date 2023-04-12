import gomea
import numpy as np
from functools import cache

# Custom fitness function resembling the concatenated deceptive trap function of size k
class SOREBChainStrongGBO(gomea.fitness.GBOFitnessFunctionRealValued):
    # Any members must be assigned in __new__ to make them accessible during instantiation of superclass
    def __new__(self, number_of_variables, value_to_reach ):
        self.rotation_angle = 45
        self.rotation_block_size = 2
        return super().__new__(self,number_of_variables,value_to_reach)

    @cache
    def number_of_subfunctions( self ):
        return self.number_of_variables - 1
    
    @cache
    def inputs_to_subfunction( self, subfunction_index ):
        return [subfunction_index, subfunction_index+1]

    def subfunction(self, subfunction_index, variables):
        rotated_variables = self.rotate_variables(variables[self.inputs_to_subfunction(subfunction_index)],self.rotation_angle)
        result = 0.0
        for i in range(self.rotation_block_size):
            result += 10.0**(6.0*(i/(self.rotation_block_size-1) ))*rotated_variables[i]*rotated_variables[i]
        return result

class SOREBChainStrongBBO(gomea.fitness.BBOFitnessFunctionRealValued):
    # Any members must be assigned in __new__ to make them accessible during instantiation of superclass
    def __new__(self, number_of_variables, value_to_reach ):
        self.rotation_angle = 45
        self.rotation_block_size = 2
        return super().__new__(self,number_of_variables,value_to_reach)

    def objective_function(self, objective_index, variables):
        result = 0.0
        for i in range(self.number_of_variables-1):
            rotated_variables = self.rotate_variables(variables[i:i+2],self.rotation_angle)
            for i in range(self.rotation_block_size):
                result += 10.0**(6.0*(i/(self.rotation_block_size-1) ))*rotated_variables[i]*rotated_variables[i]
        return result
    
    def constraint_function(self, variables):
        return 0.0