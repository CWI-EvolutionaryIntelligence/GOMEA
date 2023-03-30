import gomea
import numpy as np
import cython
from functools import cache

# Custom fitness function resembling the concatenated deceptive trap function of size k
class DeceptiveTrapFunction(gomea.fitness.GBOFitnessFunctionDiscrete):
    k : cython.int
    
    # Any members must be assigned in __new__ to make them accessible during instantiation of superclass
    def __new__(self, number_of_variables, k ):
        assert( number_of_variables % k == 0 )
        self.k = k
        value_to_reach = number_of_variables
        return super().__new__(self,number_of_variables,value_to_reach)

    @cache
    def number_of_subfunctions( self ) -> cython.int:
        return self.number_of_variables // self.k
    
    @cache
    def inputs_to_subfunction( self, subfunction_index : cython.int ) -> np.ndarray:
        return range(self.k*subfunction_index,self.k*subfunction_index+self.k)

    def subfunction(self, subfunction_index : cython.int, variables : np.ndarray ):
        trap_variables : np.ndarray = variables[self.inputs_to_subfunction(subfunction_index)]
        unitation : cython.int = np.sum(trap_variables)
        if unitation == self.k:
            return unitation
        else:
            return self.k - unitation - 1