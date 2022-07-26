#import DiscreteGOMEA
#import RealValuedGOMEA
import GOMEA

#rvgom = RealValuedGOMEA.pyRealValuedGOMEA(problem_index=0, number_of_variables=100)
rvgom = GOMEA.RealValuedGOMEA(problem_index=0, number_of_variables=100)
rvgom.run()

#gom = DiscreteGOMEA.pyDiscreteGOMEA(problem_index=0, number_of_variables=100, max_time=60)
#gom.run()


