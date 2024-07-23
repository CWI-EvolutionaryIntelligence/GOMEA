import numpy as np
import gc
import gomea

# Custom fitness function resembling the Rosenbrock function
class CustomRastriginFunction(gomea.fitness.GBOFitnessFunctionRealValued):
    def number_of_subfunctions( self ):
        return self.number_of_variables
    
    def inputs_to_subfunction( self, subfunction_index ):
        return [subfunction_index]

    def objective_function( self, objective_index, fitness_buffers ):
        return 10.0 * self.number_of_variables + fitness_buffers[0]

    def subfunction(self, subfunction_index, variables):
        x = variables[subfunction_index]
        return x*x - 10.0*np.cos(2.0*np.pi*x)
    
class CustomSphere(gomea.fitness.BBOFitnessFunctionRealValued):
    def objective_function( self, objective_index, variables ):
        return np.sum(variables**2)
    
    def constraint_function( self, variables ):
        return 0


# Custom fitness function resembling the concatenated deceptive trap function of size k
class CustomTrapFunction(gomea.fitness.GBOFitnessFunctionDiscrete):
    # Any members must be assigned in __new__ to make them accessible during instantiation of superclass
    def __new__(self, number_of_variables, k, value_to_reach):
        assert( number_of_variables % k == 0 )
        self.k = k
        return super().__new__(self,number_of_variables,value_to_reach)

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

from benchmarks.soreb_chain.PythonSOREBChainStrong import SOREBChainStrongGBO as PythonSOREBChainStrongGBO
from benchmarks.rosenbrock.PythonRosenbrock import RosenbrockFunctionGBO as PythonRosenbrockFunctionGBO
dim = 20
vtr = 1e-10
#frv = gomea.fitness.SphereFunction(10000,value_to_reach=vtr)
frv = gomea.fitness.RosenbrockFunction(dim,value_to_reach=vtr)
#frv = CustomSphere(dim,value_to_reach=vtr)
#frv = gomea.fitness.SOREBChainStrong(dim,value_to_reach=vtr)
#frv = PythonSOREBChainStrongGBO(dim,value_to_reach=vtr)

#lm = gomea.linkage.Univariate()
#lm = gomea.linkage.BlockMarginalProduct(2)
#lm = gomea.linkage.Full()
#lm = gomea.linkage.LinkageTree(maximum_set_size=4)
#lm = gomea.linkage.LinkageTree()
#lm = gomea.linkage.StaticLinkageTree(maximum_set_size=10)
lm = gomea.linkage.UCondFG()
#lm = gomea.linkage.FromFile("FOS.in")
#lm = gomea.linkage.Full()
#lm = gomea.linkage.StaticLinkageTree()
nsucc = 0
nruns = 0
#rvgom = gomea.RealValuedGOMEA(fitness=frv,linkage_model=lm,max_number_of_evaluations=10000000)
rvgom = gomea.RealValuedGOMEA(fitness=frv,linkage_model=lm, max_number_of_evaluations=5000000,base_population_size=40,max_number_of_populations=1,
                              lower_init_range=-115,upper_init_range=-100,fitness_variance_tolerance=0.0)
#rvgom = gomea.RealValuedGOMEA(fitness=frv,linkage_model=lm,lower_init_range=0,upper_init_range=1, max_number_of_populations=1, base_population_size=100, max_number_of_evaluations=10000000)
#rvgom = gomea.RealValuedGOMEA(fitness=frv,linkage_model=lm,lower_init_range=0,upper_init_range=1, max_number_of_populations=1, base_population_size=100, max_number_of_evaluations=10000000)
#rvgom = gomea.RealValuedGOMEA(fitness=frv,linkage_model=lm,lower_init_range=-115,upper_init_range=-100, max_number_of_evaluations=10000000)
print("ObjVal\tNumEvaluations\tGenerations\tTime(s)")
evals = []
for i in range(nruns):
    result = rvgom.run()
    print(result['best_obj_val'][-1],result['evaluations'][-1],result['generation'][-1],result['time'][-1])
    if result['best_obj_val'][-1] < frv.value_to_reach:
        nsucc += 1
    evals.append(result['evaluations'][-1])
print(nsucc,"/",nruns," successes")
#print(np.median(evals))

#result.printFinalStatistics()
#result.printAllStatistics()

#from benchmarks.maxcut.PythonMaxCut import MaxCutGBO as PythonMaxCutGBO
#fd = PythonMaxCutGBO("benchmarks/problem_data/maxcut/set0c/n0000784i00.txt","benchmarks/problem_data/maxcut/set0c/n0000784i00.bkv")

#fd = gomea.fitness.OneMaxFunction(1000)
#lm = gomea.linkage.Custom(file="FOS.in")
#lm = gomea.linkage.Custom(fos=[range(0,5),range(5,10),range(10,15),range(15,20)])
#lm = gomea.linkage.StaticLinkageTree()
lm = gomea.linkage.UCondHG()
#lm = gomea.linkage.LinkageTree()
#lm = gomea.linkage.BlockMarginalProduct(block_size=5)
#fd = CustomTrapFunction(20,k=5,value_to_reach=20)
#dim = 16
#fd = gomea.fitness.DeceptiveTrapFunction(dim,trap_size=5)
#fd = gomea.fitness.MaxCut("benchmarks/problem_data/maxcut/set0c/n0000049i00.txt","benchmarks/problem_data/maxcut/set0c/n0000049i00.bkv")
fd = gomea.fitness.MaxCut("benchmarks/problem_data/maxcut/set0c/n0000016i00.txt","benchmarks/problem_data/maxcut/set0c/n0000016i00.bkv")
#fd = CustomTrapFunction(dim,k=5,value_to_reach=dim)
nsucc = 0
nruns = 1
print("ObjVal\tFEV(G)\tFEV(B)\tTime(s)")
for i in range(nruns):
    v = True 
    #dgom = gomea.DiscreteGOMEA(fitness=fd,linkage_model=lm,max_number_of_evaluations=1000000,random_seed=i,base_population_size=20,max_number_of_populations=1,verbose=v)
    dgom = gomea.DiscreteGOMEA(fitness=fd,linkage_model=lm,max_number_of_evaluations=1000000,base_population_size=20,max_number_of_populations=1,verbose=v)
    result = dgom.run()
    #result.printFinalStatistics()
    print(result['best_obj_val'][-1],result['evaluations'][-1],result['evaluations_black_box'][-1],result['time'][-1])
#result.printFinalStatistics()
#result.printAllStatistics()
