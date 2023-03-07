import cython
import numpy as np
from gomea.fitness import PythonFitnessFunctionRealValued

class YourCythonFitnessFunction(PythonFitnessFunctionRealValued):
    def number_of_subfunctions( self ) -> cython.int:
        return self.number_of_variables-1
    
    def inputs_to_subfunction( self, subfunction_index : cython.int ) -> np.ndarray:
        arr : np.ndarray
        arr = np.array([subfunction_index, subfunction_index+1], np.int32)
        return arr

    def subfunction( self, subfunction_index : cython.int, variables: np.ndarray ) -> cython.double:
        x : cython.double
        y : cython.double
        x = variables[subfunction_index]
        y = variables[subfunction_index+1]
        return 100*(y-x*x)*(y-x*x) + (1.0-x)*(1.0-x);
