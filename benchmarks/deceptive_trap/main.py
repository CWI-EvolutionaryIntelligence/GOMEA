import numpy as np
import gomea
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import gc 
import json
import pathos.multiprocessing as mp 

def writeResultsToJson(input,tag):
    with open("results/"+tag+".json",'w') as fp:
        json.dump(input,fp)

def readResultsFromJson(tag):
    filename = "results/"+tag+".json"
    try:
        with open(filename,'r') as fp:
            res = json.load(fp)
            for k,sub_dict in res.items():
                res[k] = {int(dim):[float(r) for r in res] for dim,res in sub_dict.items()}
            return res
    except Exception as e:
        print(e)
        return {}

def resultsExist(tag):
    filename = "results/"+tag+".json"
    try:
        with open(filename,'r') as fp:
            return True
    except FileNotFoundError:
        return False

from PythonTrapFunction import DeceptiveTrapFunctionGBO as PythonTrapFunctionGBO
from PythonTrapFunction import DeceptiveTrapFunctionBBO as PythonTrapFunctionBBO
from CythonTrapFunction import DeceptiveTrapFunction as CythonTrapFunction
from CythonTrapFunctionCDEF import DeceptiveTrapFunction as CythonTrapFunctionCDEF
from gomea.fitness import DeceptiveTrapFunction as CppTrapFunctionGBO
from gomea.fitness import DeceptiveTrapFunctionBBO as CppTrapFunctionBBO

def getFitnessFunctionTrap5(prob_tag,lang_tag,dim):
    if prob_tag == 'trap5':
        if lang_tag == 'py':
            return PythonTrapFunctionGBO(dim,5)
        elif lang_tag == 'py_bbo':
            return PythonTrapFunctionBBO(dim,5)
        elif lang_tag == 'cy':
            return CythonTrapFunction(dim,5)
        elif lang_tag == 'cyc':
            return CythonTrapFunctionCDEF(dim,5)
        elif lang_tag == 'cpp':
            return CppTrapFunctionGBO(dim,5)
        elif lang_tag == 'cpp_bbo':
            return CppTrapFunctionBBO(dim,5)
    return None

def getLinkageModel(lm_tag):
    if lm_tag == 'bm5':
        return gomea.linkage.BlockMarginalProduct(block_size=5)
    elif lm_tag == 'lt':
        return gomea.linkage.LinkageTree(maximum_set_size=100)
    elif lm_tag == 'slt':
        return gomea.linkage.StaticLinkageTree()

def runDiscreteGOMEA(prob_tag,lang_tag,lm_tag,dim,max_evals):
    ff = getFitnessFunctionTrap5(prob_tag,lang_tag,dim)
    lm = getLinkageModel(lm_tag)
    dgom = gomea.DiscreteGOMEA(fitness=ff,linkage_model=lm,max_evals=max_evals)
    result = dgom.run()
    del ff
    del lm
    del dgom
    gc.collect()
    return result

def runExperiment(prob_tag,lang_tag,lm_tag,max_dim=1000,rerun=False,nthreads=1,nruns=10):
    full_tag = prob_tag+"_"+lang_tag+"_"+lm_tag
    dim = 10
    cur_results = readResultsFromJson(full_tag)
    fev_budget = 100000000
    fev_results = {}
    time_results = {}
    bestf_results = {}
    if not rerun and 'fev' in cur_results.keys():
        fev_results = cur_results['fev']
    if not rerun and 't' in cur_results.keys():
        time_results = cur_results['t']
    if not rerun and 'bestf' in cur_results.keys():
        bestf_results = cur_results['bestf']
    while dim < max_dim:
        if rerun or (dim not in fev_results.keys()) or (dim not in time_results.keys()) or (dim not in bestf_results.keys()):
            nsucc = 0
            fev_results[dim] = []
            time_results[dim] = []
            bestf_results[dim] = []
            
            thr_pool = mp.Pool(nthreads)
            tup = (prob_tag,lang_tag,lm_tag,dim,fev_budget)
            args = ((tup,) * nruns )
            print(dim,full_tag)
            results = thr_pool.starmap(runDiscreteGOMEA,args)
            for res in results:
                fev_results[dim].append(res['evaluations'])
                time_results[dim].append(max(0.001,res['time']))
                bestf_results[dim].append(res['best_obj_val'])
                if res['best_obj_val'] == dim:
                    nsucc += 1
            results = {"t" : time_results, "fev" : fev_results, "bestf" : bestf_results}
            writeResultsToJson(results,full_tag)
            print("Dim",dim," ",nsucc/nruns," suc ",np.median(time_results[dim]),"sec",np.median(fev_results[dim]),"fev")
        dim = 2*dim

def runExperimentsTrap5BlockMarginal(max_dim=1000,rerun=False,nthreads=1,nruns=10):
    runExperiment('trap5','cy','bm5',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns)
    #runExperiment('trap5','cyc','bm5',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns)
    runExperiment('trap5','cpp','bm5',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns)
    runExperiment('trap5','cpp_bbo','bm5',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns)
    runExperiment('trap5','py','bm5',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns)
    runExperiment('trap5','py_bbo','bm5',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns)

def runExperimentsTrap5LT(max_dim=1000,rerun=False,nthreads=1,nruns=10):
    runExperiment('trap5','cy','lt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns)
    #runExperiment('trap5','cyc','lt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns)
    runExperiment('trap5','cpp','lt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns)
    runExperiment('trap5','cpp_bbo','lt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns)
    runExperiment('trap5','py','lt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns)
    runExperiment('trap5','py_bbo','lt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns)

def runExperimentsTrap5StaticLT(max_dim=1000,rerun=False,nthreads=1,nruns=10):
    runExperiment('trap5','py','slt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns)
    runExperiment('trap5','cy','slt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns)
    #runExperiment('trap5','cyc','slt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns)
    runExperiment('trap5','cpp','slt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns)

def addResultToPlot(results,label):
    dimensions = np.sort(np.int_(list(results.keys())))
    medians = [np.median(results[dim]) for dim in dimensions]
    err_high = [np.percentile(results[dim],90)-np.median(results[dim]) for dim in dimensions] 
    err_low = [np.median(results[dim])-np.percentile(results[dim],10) for dim in dimensions] 
    plt.errorbar(dimensions,medians,yerr=[err_low,err_high],capsize=3,label=label)

def addTimePerEvalResultToPlot(time_results,fev_results,label):
    dimensions = np.sort(list(time_results.keys()))
    tpf_results = {dim:np.divide(time_results[dim],fev_results[dim]) for dim in dimensions}
    medians = [np.median(tpf_results[dim]) for dim in dimensions]
    err_high = [np.percentile(tpf_results[dim],90)-np.median(tpf_results[dim]) for dim in dimensions] 
    err_low = [np.median(tpf_results[dim])-np.percentile(tpf_results[dim],10) for dim in dimensions] 
    plt.errorbar(dimensions,medians,yerr=[err_low,err_high],capsize=3,label=label)

def plotFEVScalability(results,tag='',max_dim=1000):
    fig,ax = plt.subplots(dpi=150)
    for label,res in results:
        addResultToPlot(res['fev'],label=label)
    fig.patch.set_facecolor('white')
    plt.yscale('log')
    plt.xscale('log',base=10)
    plt.xlabel("Dimensionality")
    plt.ylabel("FEV")
    plt.xticks([10*2**i for i in range(20)])
    #plt.ylim(100,50000)
    plt.xlim(8,max_dim)
    plt.grid(linewidth=0.5)
    formatter = plt.ScalarFormatter()
    formatter.set_scientific(False)
    ax.xaxis.set_major_formatter(formatter)
    ax.legend(loc='upper left')
    plt.savefig('plots/fev_'+tag+'.pdf')

def plotTimeScalability(results,tag='',max_dim=1000):
    fig,ax = plt.subplots(dpi=150)
    for label,res in results:
        addResultToPlot(res['t'],label=label)
    fig.patch.set_facecolor('white')
    plt.yscale('log')
    plt.xscale('log',base=10)
    plt.xlabel("Dimensionality")
    plt.ylabel("Time (s)")
    plt.xticks([10*2**i for i in range(20)])
    #plt.ylim(0.0005,30)
    plt.xlim(8,max_dim)
    plt.grid(linewidth=0.5)
    formatter = plt.ScalarFormatter()
    formatter.set_scientific(False)
    ax.xaxis.set_major_formatter(formatter)
    ax.legend(loc='upper left')
    plt.savefig('plots/t_'+tag+'.pdf')

def plotTimePerFEV(results,tag='',max_dim=1000):
    fig,ax = plt.subplots(dpi=150)
    for label,res in results:
        #print(label,res)
        addTimePerEvalResultToPlot(res['t'],res['fev'],label=label)
    fig.patch.set_facecolor('white')
    plt.yscale('log')
    plt.xscale('log',base=10)
    plt.xlabel("Dimensionality")
    plt.ylabel("Time (s) / FEV")
    plt.xticks([10*2**i for i in range(20)])
    #plt.ylim(0.0005,30)
    plt.xlim(8,max_dim)
    plt.grid(linewidth=0.5)
    formatter = plt.ScalarFormatter()
    formatter.set_scientific(False)
    ax.xaxis.set_major_formatter(formatter)
    ax.legend(loc='upper left')
    plt.savefig('plots/tperfev_'+tag+'.pdf')

def runExperiments():
    rerun_experiments = False 
    nthreads = 15
    nruns = 30
    runExperimentsTrap5LT(max_dim=10000,rerun=rerun_experiments,nthreads=nthreads,nruns=nruns)
    runExperimentsTrap5BlockMarginal(max_dim=100000,rerun=rerun_experiments,nthreads=nthreads,nruns=nruns)
    runExperimentsTrap5StaticLT(max_dim=10000,rerun=rerun_experiments,nthreads=nthreads,nruns=nruns)
    #runExperiment('trap5','cpp','slt',max_dim=5000,rerun=True,nthreads=nthreads)
    #runExperiment('trap5','cpp','lt',max_dim=1000,rerun=True,nthreads=nthreads)
    #runExperiment('trap5','py','lt',max_dim=5000,rerun=True,nthreads=nthreads)
    #runExperiment('trap5','cpp_bbo','lt',max_dim=5000,rerun=True,nthreads=nthreads)

def makePlots():
    res_py = readResultsFromJson('trap5_py_bm5')
    res_cy = readResultsFromJson('trap5_cy_bm5')
    #res_cyc = readResultsFromJson('trap5_cyc_bm5')
    res_cpp = readResultsFromJson('trap5_cpp_bm5')
    res_pybbo = readResultsFromJson('trap5_py_bbo_bm5')
    res_cppbbo = readResultsFromJson('trap5_cpp_bbo_bm5')
    results = [("Python",res_py),("Cython",res_cy),("C++",res_cpp),("Python(BBO)",res_pybbo),("C++(BBO)",res_cppbbo)]
    plotTimePerFEV(results,tag='trap5_bm5',max_dim=10000)
    plotTimeScalability(results,tag='trap5_bm5',max_dim=10000)
    plotFEVScalability(results,tag='trap5_bm5',max_dim=10000)

    res_py = readResultsFromJson('trap5_py_lt')
    res_cy = readResultsFromJson('trap5_cy_lt')
    #res_cyc = readResultsFromJson('trap5_cyc_lt')
    res_cpp = readResultsFromJson('trap5_cpp_lt')
    res_pybbo = readResultsFromJson('trap5_py_bbo_lt')
    res_cppbbo = readResultsFromJson('trap5_cpp_bbo_lt')
    results = [("Python",res_py),("Cython",res_cy),("C++",res_cpp),("Python(BBO)",res_pybbo),("C++(BBO)",res_cppbbo)]
    plotTimePerFEV(results,tag='trap5_lt',max_dim=10000)
    plotTimeScalability(results,tag='trap5_lt',max_dim=10000)
    plotFEVScalability(results,tag='trap5_lt',max_dim=10000)
    
    res_py = readResultsFromJson('trap5_py_slt')
    res_cy = readResultsFromJson('trap5_cy_slt')
    res_cpp = readResultsFromJson('trap5_cpp_slt')
    res_pybbo = readResultsFromJson('trap5_py_bbo_lt')
    res_cppbbo = readResultsFromJson('trap5_cpp_bbo_lt')
    results = [("Python",res_py),("Cython",res_cy),("C++",res_cpp),("Python(BBO)",res_pybbo),("C++(BBO)",res_cppbbo)]
    plotTimePerFEV(results,tag='trap5_slt',max_dim=10000)
    plotTimeScalability(results,tag='trap5_slt',max_dim=10000)
    plotFEVScalability(results,tag='trap5_slt',max_dim=10000)

if __name__ == '__main__':
    runExperiments()
    makePlots()