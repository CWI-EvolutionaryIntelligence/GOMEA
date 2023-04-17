import gomea
import numpy as np
from functools import cache

# Custom fitness function resembling the concatenated deceptive trap function of size k
class RosenbrockFunctionGBO(gomea.fitness.GBOFitnessFunctionRealValued):
    # Any members must be assigned in __new__ to make them accessible during instantiation of superclass
    @cache
    def number_of_subfunctions( self ):
        return self.number_of_variables - 1
    
    @cache
    def inputs_to_subfunction( self, subfunction_index ):
        return [subfunction_index, subfunction_index+1]

    def subfunction(self, subfunction_index, variables):
        x = variables[subfunction_index]
        y = variables[subfunction_index+1]
        return( 100*(y-x*x)*(y-x*x) + (1.0-x)*(1.0-x) )

class RosenbrockFunctionBBO(gomea.fitness.BBOFitnessFunctionRealValued):
    # Any members must be assigned in __new__ to make them accessible during instantiation of superclass
    def objective_function(self, objective_index, variables):
        result = 0.0
        for i in range(self.number_of_variables-1):
            x = variables[i]
            y = variables[i+1]
            result += 100*(y-x*x)*(y-x*x) + (1.0-x)*(1.0-x)
        return result
    
    def constraint_function(self, variables):
        return 0.0