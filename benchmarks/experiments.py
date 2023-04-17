import gomea
import gc
import numpy as np
from output import readResultsFromJson, writeResultsToJson
import pathos.multiprocessing as mp 

from deceptive_trap.PythonTrapFunction import DeceptiveTrapFunctionGBO as PythonTrapFunctionGBO
from deceptive_trap.PythonTrapFunction import DeceptiveTrapFunctionBBO as PythonTrapFunctionBBO
from deceptive_trap.CythonTrapFunction import DeceptiveTrapFunction as CythonTrapFunction
from deceptive_trap.CythonTrapFunctionCDEF import DeceptiveTrapFunction as CythonTrapFunctionCDEF
from gomea.fitness import DeceptiveTrapFunction as CppTrapFunctionGBO
from gomea.fitness import DeceptiveTrapFunctionBBO as CppTrapFunctionBBO

from maxcut.PythonMaxCut import MaxCutGBO as PythonMaxCutGBO
from maxcut.PythonMaxCut import MaxCutBBO as PythonMaxCutBBO
from maxcut.CythonMaxCut import MaxCut as CythonMaxCut
from gomea.fitness import MaxCut as CppMaxCutGBO
from gomea.fitness import MaxCutBBO as CppMaxCutBBO

from soreb_chain.PythonSOREBChainStrong import SOREBChainStrongGBO as PythonSOREBChainStrongGBO
from soreb_chain.PythonSOREBChainStrong import SOREBChainStrongBBO as PythonSOREBChainStrongBBO
from gomea.fitness import SOREBChainStrong as CppSOREBChainStrongGBO
from gomea.fitness import SOREBChainStrongBBO as CppSOREBChainStrongBBO

from rosenbrock.PythonRosenbrock import RosenbrockFunctionGBO as PythonRosenbrockFunctionGBO
from rosenbrock.PythonRosenbrock import RosenbrockFunctionBBO as PythonRosenbrockFunctionBBO
from gomea.fitness import RosenbrockFunction as CppRosenbrockFunctionGBO
from gomea.fitness import RosenbrockFunctionBBO as CppRosenbrockFunctionBBO

def getFitnessFunction(prob_tag,lang_tag,dim=-1,input_file="",vtr_file=""):
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
    elif prob_tag == 'maxcut':
        if lang_tag == 'py':
            return PythonMaxCutGBO(input_file,vtr_file)
        elif lang_tag == 'py_bbo':
            return PythonMaxCutBBO(input_file,vtr_file)
        elif lang_tag == 'cy':
            return CythonMaxCut(input_file,vtr_file)
        elif lang_tag == 'cpp':
            return CppMaxCutGBO(input_file,vtr_file)
        elif lang_tag == 'cpp_bbo':
            return CppMaxCutBBO(input_file,vtr_file)
    elif prob_tag == 'rosenbrock':
        if lang_tag == 'py':
            return PythonRosenbrockFunctionGBO(dim,1e-10)
        elif lang_tag == 'py_bbo':
            return PythonRosenbrockFunctionBBO(dim,1e-10)
        elif lang_tag == 'cpp':
            return CppRosenbrockFunctionGBO(dim,1e-10)
        elif lang_tag == 'cpp_bbo':
            return CppRosenbrockFunctionBBO(dim,1e-10)
    elif prob_tag == 'sorebchainstrong':
        if lang_tag == 'py':
            return PythonSOREBChainStrongGBO(dim,1e-10)
        elif lang_tag == 'py_bbo':
            return PythonSOREBChainStrongBBO(dim,1e-10)
        elif lang_tag == 'cpp':
            return CppSOREBChainStrongGBO(dim,1e-10)
        elif lang_tag == 'cpp_bbo':
            return CppSOREBChainStrongBBO(dim,1e-10)
    return None

def getLinkageModel(lm_tag):
    if lm_tag == 'bm5':
        return gomea.linkage.BlockMarginalProduct(block_size=5)
    elif lm_tag == 'uni':
        return gomea.linkage.Univariate()
    elif lm_tag == 'full':
        return gomea.linkage.Full()
    elif lm_tag == 'ucondhg':
        return gomea.linkage.UCondHG()
    elif lm_tag == 'lt':
        return gomea.linkage.LinkageTree()
    elif lm_tag == 'lt10':
        return gomea.linkage.LinkageTree(maximum_set_size=10)
    elif lm_tag == 'lt50':
        return gomea.linkage.LinkageTree(maximum_set_size=50)
    elif lm_tag == 'lt100':
        return gomea.linkage.LinkageTree(maximum_set_size=100)
    elif lm_tag == 'slt':
        return gomea.linkage.StaticLinkageTree()
    elif lm_tag == 'slt10':
        return gomea.linkage.StaticLinkageTree(maximum_set_size=10)
    elif lm_tag == 'slt50':
        return gomea.linkage.StaticLinkageTree(maximum_set_size=50)
    elif lm_tag == 'slt100':
        return gomea.linkage.StaticLinkageTree(maximum_set_size=100)

def getDomain(prob_tag):
    if prob_tag == 'trap5':
        return 'd'
    elif prob_tag == 'maxcut':
        return 'd'
    elif prob_tag == 'rosenbrock':
        return 'rv'
    elif prob_tag == 'sorebchainstrong':
        return 'rv'
    return None

def runGOMEA(prob_tag,lang_tag,lm_tag,dim,max_evals,max_sec,input_file,vtr_file):
    ff = getFitnessFunction(prob_tag,lang_tag,dim,input_file,vtr_file)
    lm = getLinkageModel(lm_tag)
    result = None
    if getDomain(prob_tag) == 'd':
        dgom = gomea.DiscreteGOMEA(fitness=ff,linkage_model=lm,max_number_of_evaluations=max_evals,max_number_of_seconds=max_sec)
        result = dgom.run()
        del dgom
    elif getDomain(prob_tag) == 'rv':
        rvgom = gomea.RealValuedGOMEA(fitness=ff,linkage_model=lm,max_number_of_evaluations=max_evals,max_number_of_seconds=max_sec)
        result = rvgom.run()
        del rvgom
    del ff
    del lm
    gc.collect()
    return result

def runExperiment(prob_tag,lang_tag,lm_tag,max_dim=1000,max_sec=-1,rerun=False,nthreads=1,nruns=10,input_file="",vtr_file=""):
    dim = 10
    cur_results = readResultsFromJson(prob_tag,lang_tag,lm_tag)
    fev_budget = 10000000
    fev_results = {}
    time_results = {}
    eval_time_results = {}
    bestf_results = {}
    succ_rate_results = {}
    if not rerun and 'fev' in cur_results.keys():
        fev_results = cur_results['fev']
    if not rerun and 't' in cur_results.keys():
        time_results = cur_results['t']
    if not rerun and 'bestf' in cur_results.keys():
        bestf_results = cur_results['bestf']
    if not rerun and 'eval_t' in cur_results.keys():
        eval_time_results = cur_results['eval_t']
    if not rerun and 'succ_rate' in cur_results.keys():
        succ_rate_results = cur_results['succ_rate']
    if len(input_file) > 0:
        dim = int(input_file[-14:-7])
    while dim <= max_dim or max_dim == -1:
        if rerun or (dim not in fev_results.keys()) or (dim not in time_results.keys()) or (dim not in eval_time_results.keys()) or (dim not in bestf_results.keys()) or (dim not in succ_rate_results.keys()):
            nsucc = 0
            fev_results[dim] = []
            time_results[dim] = []
            eval_time_results[dim] = []
            bestf_results[dim] = []
            succ_rate_results[dim] = 0

            thr_pool = mp.Pool(nthreads)
            tup = (prob_tag,lang_tag,lm_tag,dim,fev_budget,max_sec,input_file,vtr_file)
            args = ((tup,) * nruns )
            print(dim,prob_tag,lang_tag,lm_tag)
            results = thr_pool.starmap(runGOMEA,args)
            nsucc = 0
            vtr = getFitnessFunction(prob_tag,lang_tag,dim,input_file,vtr_file).value_to_reach
            domain = getDomain(prob_tag)
            for res in results:
                fev_results[dim].append(res['evaluations'][-1])
                time_results[dim].append(res['time'][-1])
                eval_time_results[dim].append(res['eval_time'][-1])
                bestf_results[dim].append(res['best_obj_val'][-1])
                if (domain == 'd' and res['best_obj_val'][-1] >= vtr) or (domain == 'rv' and res['best_obj_val'][-1] <= vtr):
                    nsucc += 1
            succ_rate_results[dim] = nsucc / nruns
            results = {"t" : time_results, "fev" : fev_results, "bestf" : bestf_results, "eval_t" : eval_time_results, "succ_rate" : succ_rate_results }
            writeResultsToJson(results,prob_tag,lang_tag,lm_tag)
            print("Dim",dim," ",nsucc/nruns," suc ",np.median(time_results[dim]),"sec",np.median(fev_results[dim]),"fev")
        dim = 2*dim
        if len(input_file) > 0:
            break
