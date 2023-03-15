import cython
import numpy as np
import gomea
from CustomFitness import CustomTrapFunction

nruns = 30

dim = 10
lm = gomea.linkage.BlockMarginalProduct(block_size=5)
fev_results = {}
time_results = {}
succ_results = {}
while dim < 400:
    fd = CustomTrapFunction(dim,k=5,value_to_reach=dim)
    nsucc = 0
    fev_results[dim] = []
    time_results[dim] = []
    succ_results[dim] = 0
    for i in range(nruns):
        dgom = gomea.DiscreteGOMEA(fitness=fd,linkage_model=lm,max_evals=1000000)
        result = dgom.run()
        #print(result['best_obj_val'],result['evaluations'],result['time'])
        if result['best_obj_val'] == dim:
            succ_results[dim] = succ_results[dim] + 1
        fev_results[dim].append(result['evaluations'])
        time_results[dim].append(result['time'])
    print("D",dim," ",succ_results[dim]," suc ",np.median(time_results[dim]),"sec",np.median(fev_results[dim]),"fev")
    dim = 2*dim