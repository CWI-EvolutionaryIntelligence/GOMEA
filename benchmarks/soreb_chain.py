from experiments import runExperiment
from output import readResultsFromJson
from plotting import plotEvalTimeScalability, plotTimePerFEV, plotTimeScalability, plotFEVScalability

def runExperimentsSOREBChainFull(max_dim=1000,rerun=False,nthreads=1,nruns=10,max_sec=-1):
    runExperiment('sorebchainstrong','cpp','full',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('sorebchainstrong','cpp_bbo','full',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('sorebchainstrong','py','full',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('sorebchainstrong','py_bbo','full',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)

def runExperimentsSOREBChainLT(max_dim=1000,rerun=False,nthreads=1,nruns=10,max_sec=-1):
    runExperiment('sorebchainstrong','cpp','lt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('sorebchainstrong','cpp_bbo','lt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('sorebchainstrong','py','lt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('sorebchainstrong','py_bbo','lt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('sorebchainstrong','py_bbo','lt50',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('sorebchainstrong','cpp_bbo','lt50',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)

def runExperimentsSOREBChainSLT(max_dim=1000,rerun=False,nthreads=1,nruns=10,max_sec=-1):
    runExperiment('sorebchainstrong','py','slt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('sorebchainstrong','cpp','slt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)

def runExperimentsSOREBChainUCond(max_dim=1000,rerun=False,nthreads=1,nruns=10,max_sec=-1):
    runExperiment('sorebchainstrong','cpp','ucondhg',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('sorebchainstrong','py','ucondhg',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)

def runExperiments():
    rerun_experiments = False 
    nthreads = 15
    nruns = 30
    max_sec = 3600
    runExperimentsSOREBChainFull(max_dim=500,rerun=rerun_experiments,nthreads=nthreads,nruns=nruns)
    runExperimentsSOREBChainLT(max_dim=500,rerun=rerun_experiments,nthreads=nthreads,nruns=nruns)
    runExperimentsSOREBChainSLT(max_dim=500,rerun=rerun_experiments,nthreads=nthreads,nruns=nruns)
    runExperimentsSOREBChainUCond(max_dim=500,rerun=rerun_experiments,nthreads=nthreads,nruns=nruns,max_sec=max_sec)

def makePlots(out_dir="plots/"):
    res_py_lt = readResultsFromJson('sorebchainstrong','py','lt')
    res_py_lt100 = readResultsFromJson('sorebchainstrong','py','lt100')
    res_cpp_lt = readResultsFromJson('sorebchainstrong','cpp','lt')
    res_cpp_lt50 = readResultsFromJson('sorebchainstrong','cpp','lt50')
    res_cpp_lt100 = readResultsFromJson('sorebchainstrong','cpp','lt100')
    res_pybbo_lt100 = readResultsFromJson('sorebchainstrong','py_bbo','lt100')
    res_cppbbo_lt = readResultsFromJson('sorebchainstrong','cpp_bbo','lt')
    res_cppbbo_lt50 = readResultsFromJson('sorebchainstrong','cpp_bbo','lt50')
    res_cppbbo_lt100 = readResultsFromJson('sorebchainstrong','cpp_bbo','lt100')
    results = [("Python",res_py_lt100),("Cython",res_cy_lt100),("C++",res_cpp_lt100),("Python(BBO)",res_pybbo_lt100),("C++(BBO)",res_cppbbo_lt100)]
    plotTimePerFEV(results,out_dir=out_dir,tag='sorebchainstrong_lt100',max_dim=10000,xticks=xticks)
    plotTimeScalability(results,out_dir=out_dir,tag='sorebchainstrong_lt100',max_dim=10000,xticks=xticks)
    plotEvalTimeScalability(results,out_dir=out_dir,tag='sorebchainstrong_lt100',max_dim=10000,xticks=xticks)
    plotFEVScalability(results,out_dir=out_dir,tag='sorebchainstrong_lt100',max_dim=10000,xticks=xticks)
    
    res_py_slt = readResultsFromJson('sorebchainstrong','py','slt')
    #res_cy_slt = readResultsFromJson('sorebchainstrong','cy','slt')
    res_cpp_slt = readResultsFromJson('sorebchainstrong','cpp','slt')
    res_pybbo_lt = readResultsFromJson('sorebchainstrong','py_bbo','lt')
    res_cppbbo_lt = readResultsFromJson('sorebchainstrong','cpp_bbo','lt')
    #results = [("Python",res_py_slt),("Cython",res_cy_slt),("C++",res_cpp_slt),("Python(BBO)",res_pybbo_lt),("C++(BBO)",res_cppbbo_lt)]
    #plotTimePerFEV(results,out_dir=out_dir,tag='sorebchainstrong_slt',max_dim=10000,xticks=xticks)
    #plotEvalTimeScalability(results,out_dir=out_dir,tag='sorebchainstrong_slt',max_dim=10000,xticks=xticks)
    #plotTimeScalability(results,out_dir=out_dir,tag='sorebchainstrong_slt',max_dim=10000,xticks=xticks)
    #plotFEVScalability(results,out_dir=out_dir,tag='sorebchainstrong_slt',max_dim=10000,xticks=xticks)
    
    results = [("BBO-LT",res_cppbbo_lt),("GBO-SLT",res_cpp_slt),("GBO-LT",res_cpp_lt)]
    #results = [("BBO-LT",res_cppbbo_lt),("BBO-LT(50)",res_cppbbo_lt50),("BBO-LT(100)",res_cppbbo_lt100),("GBO-SLT",res_cpp_slt),("GBO-LT",res_cpp_lt),("GBO-LT(50)",res_cpp_lt50),("GBO-LT(100)",res_cpp_lt100)]
    plotEvalTimeScalability(results,out_dir=out_dir,tag='sorebchainstrong_cpp_lms',max_dim=10000,xticks=xticks)
    plotTimeScalability(results,out_dir=out_dir,tag='sorebchainstrong_cpp_lms',max_dim=10000,xticks=xticks)
    plotFEVScalability(results,out_dir=out_dir,tag='sorebchainstrong_cpp_lms',max_dim=10000,xticks=xticks)
    
    results = [("BBO-LT",res_pybbo_lt),("GBO-SLT",res_py_slt),("GBO-LT",res_py_lt)]
    #results = [("BBO-LT",res_cppbbo_lt),("BBO-LT(50)",res_cppbbo_lt50),("BBO-LT(100)",res_cppbbo_lt100),("GBO-SLT",res_cpp_slt),("GBO-LT",res_cpp_lt),("GBO-LT(50)",res_cpp_lt50),("GBO-LT(100)",res_cpp_lt100)]
    plotEvalTimeScalability(results,out_dir=out_dir,tag='sorebchainstrong_py_lms',max_dim=10000,xticks=xticks)
    plotTimeScalability(results,out_dir=out_dir,tag='sorebchainstrong_py_lms',max_dim=10000,xticks=xticks)
    plotFEVScalability(results,out_dir=out_dir,tag='sorebchainstrong_py_lms',max_dim=10000,xticks=xticks)

if __name__ == '__main__':
    runExperiments()
    makePlots("plots/sorebchainstrong/")