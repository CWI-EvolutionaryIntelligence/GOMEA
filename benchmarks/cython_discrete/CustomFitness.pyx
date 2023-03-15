import cython
import numpy as np
from gomea.fitness import GBOFitnessFunctionDiscrete

class CustomTrapFunction(GBOFitnessFunctionDiscrete):
    # Any members must be assigned in __new__ to make them accessible during instantiation of superclass
    def __new__(self, number_of_variables : cython.int, k : cython.int, value_to_reach : cython.double ):
        assert( number_of_variables % k == 0 )
        self.k = k
        return super().__new__(self,number_of_variables,value_to_reach)

    def number_of_subfunctions( self ) -> cython.int:
        return self.number_of_variables // self.k
    
    def inputs_to_subfunction( self, subfunction_index : cython.int ) -> np.ndarray:
        return np.arange(self.k*subfunction_index,self.k*subfunction_index+self.k)

    def subfunction(self, subfunction_index, variables) -> cython.double:
        trap_variables : np.ndarray
        unitation : cython.int
        trap_variables = variables[self.inputs_to_subfunction(subfunction_index)]
        unitation = np.sum(trap_variables)
        if unitation == self.k:
            return unitation
        else:
            return self.k - unitation - 1

