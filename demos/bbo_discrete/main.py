import gomea
import numpy as np

# Custom fitness function resembling the concatenated deceptive trap function of size k
class CustomTrapFunction(gomea.fitness.BBOFitnessFunctionDiscrete):
    # Any members must be assigned in __new__ to make them accessible during instantiation of superclass
    def __new__(self, number_of_variables, k, value_to_reach):
        assert( number_of_variables % k == 0 )
        self.k = k
        return super().__new__(self,number_of_variables,value_to_reach)

    def objective_function(self, objective_index, variables):
        f = 0
        for i in range(0,self.number_of_variables,self.k):
            trap_variables = variables[i:i+self.k]
            unitation = np.sum(trap_variables)
            if unitation == self.k:
                f += unitation
            else:
                f += self.k - unitation - 1
        return f

dim = 20
lm = gomea.linkage.BlockMarginalProduct(block_size=5)
fd = CustomTrapFunction(dim,k=5,value_to_reach=dim)
dgom = gomea.DiscreteGOMEA(fitness=fd,linkage_model=lm,max_number_of_evaluations=100000)
result = dgom.run()
result.printAllStatistics()
