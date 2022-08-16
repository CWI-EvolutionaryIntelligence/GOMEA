import gomea
import numpy as np

class YourPythonFitnessFunction(gomea.PythonFitnessFunction):
    def number_of_subfunctions( self ):
        return self.number_of_variables-1
    
    def inputs_to_subfunction( self, subfunction_index ):
        arr = np.array([subfunction_index, subfunction_index+1], np.int32)
        return arr

    def subfunction(self, subfunction_index, variables):
        x = variables[subfunction_index]
        y = variables[subfunction_index+1]
        return 100*(y-x*x)*(y-x*x) + (1.0-x)*(1.0-x);
    

#f = gomea.SphereFunction(10000,value_to_reach=1e-10)
#f = gomea.RosenbrockFunction(40,value_to_reach=1e-10)
#f = gomea.YourFitnessFunction(10000,value_to_reach=1e-10)
f = YourPythonFitnessFunction(40,value_to_reach=1e-10)
rvgom = gomea.RealValuedGOMEA(fitness=f,lower_init_range=-115,upper_init_range=-100) #, maximum_number_of_populations=1, base_population_size=20)
#rvgom = gomea.RealValuedGOMEA(fitness=f, maximum_number_of_populations=1, base_population_size=20)
result = rvgom.run()

