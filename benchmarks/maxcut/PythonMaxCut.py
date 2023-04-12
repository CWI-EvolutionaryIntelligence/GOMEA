import gomea
import numpy as np
from functools import cache

# Custom fitness function resembling the concatenated deceptive trap function of size k
class MaxCutGBO(gomea.fitness.GBOFitnessFunctionDiscrete):
    # Any members must be assigned in __new__ to make them accessible during instantiation of superclass
    def __new__(self, input_file, vtr_file ):
        number_of_variables = 0
        value_to_reach = 1e308
        self.edges = []
        with open(input_file,"r") as f:
            lines = [line.split() for line in f.readlines()]
            number_of_variables = int(lines[0][0])
            number_of_edges = int(lines[0][1])
            for i in range(number_of_edges):
                v1 = int(lines[i+1][0])-1
                v2 = int(lines[i+1][1])-1
                w = float(lines[i+1][2])
                self.edges.append((v1,v2,w))
            assert(number_of_edges == len(self.edges))
        with open(vtr_file,"r") as f:
            lines = f.readlines()
            value_to_reach = float(lines[0].strip())
        return super().__new__(self,number_of_variables,value_to_reach)

    @cache
    def number_of_subfunctions( self ):
        return len(self.edges)
    
    @cache
    def inputs_to_subfunction( self, subfunction_index ):
        return list(self.edges[subfunction_index][0:2])

    def subfunction(self, subfunction_index, variables):
        v1,v2,w = self.edges[subfunction_index]
        if variables[v1] == variables[v2]:
            return 0.0
        else:
            return w

class MaxCutBBO(gomea.fitness.BBOFitnessFunctionDiscrete):
    def __new__(self, input_file, vtr_file ):
        number_of_variables = 0
        value_to_reach = 1e308
        self.edges = []
        with open(input_file,"r") as f:
            lines = [line.split() for line in f.readlines()]
            number_of_variables = int(lines[0][0])
            number_of_edges = int(lines[0][1])
            for i in range(number_of_edges):
                v1 = int(lines[i+1][0])-1
                v2 = int(lines[i+1][1])-1
                w = float(lines[i+1][2])
                self.edges.append((v1,v2,w))
            assert(number_of_edges == len(self.edges))
        with open(vtr_file,"r") as f:
            lines = f.readlines()
            value_to_reach = float(lines[0].strip())
        return super().__new__(self,number_of_variables,value_to_reach)

    def objective_function(self, objective_index, variables):
        cut = 0.0
        for v1,v2,w in self.edges:
            if variables[v1] != variables[v2]:
                cut += w
        return cut