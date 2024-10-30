import gomea
import numpy as np

# Custom BBO fitness function resembling the Rosenbrock function
class CustomRosenbrockFunction(gomea.fitness.BBOFitnessFunctionRealValued):
    def objective_function( self, objective_index, variables ):
        f = 0
        for i in range(len(variables)-1):
            x = variables[i]
            y = variables[i+1]
            f += 100*(y-x*x)*(y-x*x) + (1.0-x)*(1.0-x)
        return f

dim = 10
frv = CustomRosenbrockFunction(dim,value_to_reach=1e-6)
lm = gomea.linkage.Univariate()
rvgom = gomea.RealValuedGOMEA(fitness=frv, linkage_model=lm, lower_init_range=-115, upper_init_range=-100, max_number_of_populations=1, base_population_size=100, max_number_of_evaluations=1000000)
result = rvgom.run()
result.printAllStatistics()
