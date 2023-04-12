import gomea
import numpy as np
import cython
from functools import cache

# Custom fitness function resembling the concatenated deceptive trap function of size k
class MaxCut(gomea.fitness.GBOFitnessFunctionDiscrete):
    edges : np.ndarray
    def __new__(self, input_file, vtr_file ):
        number_of_variables : cython.int = 0
        value_to_reach : cython.double = 1e308
        self.edges : np.ndarray = []
        with open(input_file,"r") as f:
            lines = [line.split() for line in f.readlines()]
            number_of_variables = int(lines[0][0])
            number_of_edges = int(lines[0][1])
            for i in range(number_of_edges):
                v1 = int(lines[i+1][0])-1
                v2 = int(lines[i+1][1])-1
                w = float(lines[i+1][2])
                self.edges.append((v1,v2,w))
        with open(vtr_file,"r") as f:
            lines = f.readlines()
            value_to_reach = float(lines[0].strip())
        return super().__new__(self,number_of_variables,value_to_reach)

    @cache
    def number_of_subfunctions( self ) -> cython.int:
        return len(self.edges)
    
    @cache
    def inputs_to_subfunction( self, subfunction_index ) -> np.ndarray:
        return list(self.edges[subfunction_index])

    def subfunction(self, subfunction_index, variables) -> cython.double:
        v1 : cython.int = self.edges[subfunction_index][0]
        v2 : cython.int = self.edges[subfunction_index][1]
        w : cython.double = self.edges[subfunction_index][2]
        if variables[v1] == variables[v2]:
            return 0.0
        else:
            return w