import gomea
import numpy as np

# Custom fitness function resembling the Rosenbrock function
class CustomRosenbrockFunction(gomea.fitness.PythonFitnessFunctionRealValued):
    def number_of_subfunctions( self ):
        return self.number_of_variables-1
    
    def inputs_to_subfunction( self, subfunction_index ):
        return [subfunction_index, subfunction_index+1]

    def subfunction(self, subfunction_index, variables):
        x = variables[subfunction_index]
        y = variables[subfunction_index+1]
        return 100*(y-x*x)*(y-x*x) + (1.0-x)*(1.0-x)

# Custom fitness function resembling the concatenated deceptive trap function of size k
class CustomTrapFunction(gomea.fitness.PythonFitnessFunctionDiscrete):
    # Any members must be assigned in __new__ to make them accessible during instantiation of superclass
    def __new__(self, number_of_variables, k):
        assert( number_of_variables % k == 0 )
        self.k = k
        return super().__new__(self,number_of_variables)

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

#frv = gomea.fitness.SphereFunction(10000,value_to_reach=1e-10)
#frv = gomea.fitness.RosenbrockFunction(20,value_to_reach=1e-10)
frv = CustomRosenbrockFunction(20,value_to_reach=1e-10)

#lm = gomea.linkage.MarginalProductModel(1)
lm = gomea.linkage.Univariate()
#lm = gomea.linkage.LinkageTree(maximum_set_size=10)
#lm = gomea.linkage.LinkageTree()
#lm = gomea.linkage.StaticLinkageTree(maximum_set_size=10)
#lm = gomea.linkage.UCondHG()
#lm = gomea.linkage.FromFile("FOS.in")
#lm = gomea.linkage.Full()
#lm = gomea.linkage.StaticLinkageTree()
rvgom = gomea.RealValuedGOMEA(fitness=frv,linkage_model=lm,lower_init_range=-115,upper_init_range=-100,random_seed=100, maximum_number_of_populations=1, base_population_size=20,max_evals=10000)
result = rvgom.run()
result.printFinalStatistics()

#fd = gomea.fitness.OneMaxFunction(1000)
lm = gomea.linkage.StaticLinkageTree(maximum_set_size=5)
fd = CustomTrapFunction(20,k=5)
dgom = gomea.DiscreteGOMEA(fitness=fd,linkage_model=lm,max_evals=1000)
result = dgom.run()
result.printFinalStatistics()
result.printAllStatistics()
