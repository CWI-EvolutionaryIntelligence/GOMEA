import gomea
import numpy as np
import gc

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

dim = 20
vtr = 1e-10
#frv = gomea.fitness.SphereFunction(10000,value_to_reach=vtr)
#frv = gomea.fitness.RosenbrockFunction(20,value_to_reach=vtr)
frv = CustomRastriginFunction(dim,value_to_reach=vtr)

lm = gomea.linkage.Univariate()
nsucc = 0
nruns = 5
rvgom = gomea.RealValuedGOMEA(fitness=frv,linkage_model=lm,lower_init_range=5,upper_init_range=10, max_number_of_populations=1, base_population_size=200, max_number_of_evaluations=100000)
print("ObjVal\tNumEvaluations\tTime(s)")
for i in range(nruns):
    result = rvgom.run()
    print(result['best_obj_val'][-1],result['evaluations'][-1],result['time'][-1])
    if result['best_obj_val'][-1] < vtr:
        nsucc += 1
print(nsucc,"/",nruns," successes")
#result.printFinalStatistics()
result.printAllStatistics()

#fd = gomea.fitness.OneMaxFunction(1000)
#lm = gomea.linkage.StaticLinkageTree(maximum_set_size=5)
#lm = gomea.linkage.Custom(file="FOS.in")
#lm = gomea.linkage.Custom(fos=[range(0,5),range(5,10),range(10,15),range(15,20)])
#lm = gomea.linkage.BlockMarginalProduct(block_size=5)
lm = gomea.linkage.LinkageTree()
#fd = CustomTrapFunction(20,k=5,value_to_reach=20)
dim = 320
fd = gomea.fitness.DeceptiveTrapFunction(dim,trap_size=5)
#fd = CustomTrapFunction(dim,k=5,value_to_reach=dim)
nsucc = 0
print("ObjVal\tNumEvaluations\tTime(s)")
for i in range(nruns):
    dgom = gomea.DiscreteGOMEA(fitness=fd,linkage_model=lm,max_number_of_evaluations=1000000)
    result = dgom.run()
    #result.printFinalStatistics()
    print(result['best_obj_val'][-1],result['evaluations'][-1],result['time'][-1])
    if result['best_obj_val'][-1] == dim:
        nsucc += 1
print(nsucc,"/",nruns," successes")
result.printAllStatistics()
