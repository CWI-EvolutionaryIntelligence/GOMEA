import gomea
import Fitness 

f = Fitness.SphereFunction(50,value_to_reach=1e-10)
rvgom = gomea.RealValuedGOMEA(fitness=f, problem_index=0)
result = rvgom.run()

