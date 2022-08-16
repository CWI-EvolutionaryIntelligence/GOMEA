import cython
import numpy as np
import gomea
from CustomFitness import YourCythonFitnessFunction

f = YourCythonFitnessFunction(40,value_to_reach=1e-10)
rvgom = gomea.RealValuedGOMEA(fitness=f,lower_init_range=-115,upper_init_range=-100) #, maximum_number_of_populations=1, base_population_size=20)
result = rvgom.run()

