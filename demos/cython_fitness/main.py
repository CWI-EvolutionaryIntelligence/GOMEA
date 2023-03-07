import cython
import numpy as np
import gomea
from CustomFitness import YourCythonFitnessFunction

frv = YourCythonFitnessFunction(20,value_to_reach=1e-10)
lm = gomea.linkage.Univariate()
rvgom = gomea.RealValuedGOMEA(fitness=frv, linkage_model=lm, lower_init_range=-115, upper_init_range=-100, maximum_number_of_populations=1, base_population_size=40, max_evals=100000)
result = rvgom.run()
result.printFinalStatistics()
