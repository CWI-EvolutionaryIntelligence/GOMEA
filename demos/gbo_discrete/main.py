import gomea
import numpy as np

# Custom fitness function resembling the concatenated deceptive trap function of size k
class CustomTrapFunction(gomea.fitness.GBOFitnessFunctionDiscrete):
    # Any members must be assigned in __new__ to make them accessible during instantiation of superclass
    def __new__(self, number_of_variables, k, **kwargs):
        assert( number_of_variables % k == 0 )
        self.k = k
        return super().__new__(self,number_of_variables,**kwargs)

    def number_of_subfunctions( self ):
        return self.number_of_variables // self.k
    
    def inputs_to_subfunction( self, subfunction_index ):
        return range(self.k*subfunction_index,self.k*subfunction_index+self.k)

    def subfunction(self, subfunction_index, variables):
        trap_variables = variables[self.inputs_to_subfunction(subfunction_index)]
        unitation = np.sum(trap_variables)
        if unitation == self.k:
            return unitation
        else:
            return self.k - unitation - 1

dim = 20
lm = gomea.linkage.StaticLinkageTree(maximum_set_size=5)
fd = CustomTrapFunction(dim,k=5,value_to_reach=dim)
dgom = gomea.DiscreteGOMEA(fitness=fd,linkage_model=lm,max_number_of_evaluations=100000)
result = dgom.run()
result.printAllStatistics()
