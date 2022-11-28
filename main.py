import gomea
import numpy as np

class YourPythonFitnessFunction(gomea.fitness.PythonFitnessFunctionRealValued):
    def number_of_subfunctions( self ):
        return self.number_of_variables-1
    
    def inputs_to_subfunction( self, subfunction_index ):
        arr = np.array([subfunction_index, subfunction_index+1], np.int32)
        return arr

    def subfunction(self, subfunction_index, variables):
        x = variables[subfunction_index]
        y = variables[subfunction_index+1]
        return 100*(y-x*x)*(y-x*x) + (1.0-x)*(1.0-x)
    

#frv = gomea.fitness.SphereFunction(10000,value_to_reach=1e-10)
#frv = gomea.fitness.RosenbrockFunction(20,value_to_reach=1e-10)
frv = YourPythonFitnessFunction(20,value_to_reach=1e-10)
lm = gomea.linkage.Univariate()
#lm = gomea.linkage.MarginalProductModel(1)
rvgom = gomea.RealValuedGOMEA(fitness=frv,linkage_model=lm,lower_init_range=-115,upper_init_range=-100,random_seed=100) #, maximum_number_of_populations=1, base_population_size=20)
result = rvgom.run()

fd = gomea.fitness.OneMaxFunction(1000)
dgom = gomea.DiscreteGOMEA(fitness=fd,linkage_model=lm) 
result = dgom.run()
