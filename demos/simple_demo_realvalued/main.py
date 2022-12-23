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

frv = CustomRosenbrockFunction(20,value_to_reach=1e-10)
lm = gomea.linkage.Univariate()
rvgom = gomea.RealValuedGOMEA(fitness=frv, linkage_model=lm, lower_init_range=-115, upper_init_range=-100, maximum_number_of_populations=1, base_population_size=20, max_evals=100000)
result = rvgom.run()
result.printFinalStatistics()
