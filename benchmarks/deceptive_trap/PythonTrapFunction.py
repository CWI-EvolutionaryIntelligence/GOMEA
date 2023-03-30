import gomea
import numpy as np
from functools import cache

# Custom fitness function resembling the concatenated deceptive trap function of size k
class DeceptiveTrapFunctionGBO(gomea.fitness.GBOFitnessFunctionDiscrete):
    # Any members must be assigned in __new__ to make them accessible during instantiation of superclass
    def __new__(self, number_of_variables, k ):
        assert( number_of_variables % k == 0 )
        self.k = k
        value_to_reach = number_of_variables
        return super().__new__(self,number_of_variables,value_to_reach)

    @cache
    def number_of_subfunctions( self ):
        return self.number_of_variables // self.k
    
    @cache
    def inputs_to_subfunction( self, subfunction_index ):
        return range(self.k*subfunction_index,self.k*subfunction_index+self.k)

    def subfunction(self, subfunction_index, variables):
        trap_variables = variables[self.inputs_to_subfunction(subfunction_index)]
        unitation = np.sum(trap_variables)
        if unitation == self.k:
            return unitation
        else:
            return self.k - unitation - 1

class DeceptiveTrapFunctionBBO(gomea.fitness.BBOFitnessFunctionDiscrete):
    # Any members must be assigned in __new__ to make them accessible during instantiation of superclass
    def __new__(self, number_of_variables, k ):
        assert( number_of_variables % k == 0 )
        self.k = k
        value_to_reach = number_of_variables
        return super().__new__(self,number_of_variables,value_to_reach)

    def objective_function(self, objective_index, variables):
        def trap_function(unitation):
            if unitation == self.k:
                return unitation
            else:
                return self.k - unitation - 1
        
        start_indices = np.arange(0,self.number_of_variables,self.k)
        unitations = [np.sum(variables[ind:ind+self.k]) for ind in start_indices]
        f = np.sum([trap_function(u) for u in unitations])
        return f
    
    def constraint_function(self, variables):
        return 0.0