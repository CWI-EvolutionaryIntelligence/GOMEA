from experiments import runExperiment
from output import readResultsFromJson
from plotting import plotEvalTimeScalability, plotTimePerFEV, plotTimeScalability, plotFEVScalability

def runExperimentsTrap5BlockMarginal(max_dim=1000,rerun=False,nthreads=1,nruns=10,max_sec=-1):
    runExperiment('trap5','cy','bm5',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    #runExperiment('trap5','cyc','bm5',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('trap5','cpp','bm5',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('trap5','cpp_bbo','bm5',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('trap5','py','bm5',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('trap5','py_bbo','bm5',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)

def runExperimentsTrap5LT(max_dim=1000,rerun=False,nthreads=1,nruns=10,max_sec=-1):
    runExperiment('trap5','cy','lt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    #runExperiment('trap5','cyc','lt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('trap5','cpp','lt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('trap5','cpp_bbo','lt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('trap5','py','lt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('trap5','py_bbo','lt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)

def runExperimentsTrap5StaticLT(max_dim=1000,rerun=False,nthreads=1,nruns=10,max_sec=-1):
    runExperiment('trap5','py','slt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('trap5','cy','slt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    #runExperiment('trap5','cyc','slt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('trap5','cpp','slt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)

def runExperiments():
    rerun_experiments = False 
    nthreads = 15
    nruns = 30
    max_sec = 3600
    #runExperimentsTrap5LT(max_dim=5000,rerun=rerun_experiments,nthreads=nthreads,nruns=nruns)
    #runExperimentsTrap5BlockMarginal(max_dim=10000,rerun=rerun_experiments,nthreads=nthreads,nruns=nruns)
    #runExperimentsTrap5StaticLT(max_dim=5000,rerun=rerun_experiments,nthreads=nthreads,nruns=nruns)
    #runExperimentsTrap5LT(max_dim=10000,rerun=rerun_experiments,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    #runExperimentsTrap5BlockMarginal(max_dim=10000,rerun=rerun_experiments,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    #runExperimentsTrap5StaticLT(max_dim=10000,rerun=rerun_experiments,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('trap5','py','slt',max_dim=10000,rerun=False,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('trap5','py','lt',max_dim=10000,rerun=False,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('trap5','py','lt50',max_dim=10000,rerun=False,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('trap5','py_bbo','lt',max_dim=2000,rerun=False,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('trap5','cpp','slt',max_dim=10000,rerun=False,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('trap5','cpp','lt50',max_dim=5000,rerun=False,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('trap5','cpp','lt',max_dim=10000,rerun=False,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('trap5','cpp_bbo','lt',max_dim=5000,rerun=False,nthreads=nthreads,nruns=nruns,max_sec=max_sec)

def makePlots(out_dir="plots/"):
    res_py_bm5 = readResultsFromJson('trap5','py','bm5')
    res_cy_bm5 = readResultsFromJson('trap5','cy','bm5')
    #res_cyc = readResultsFromJson('trap5','cyc','bm5')
    res_cpp_bm5 = readResultsFromJson('trap5','cpp','bm5')
    res_pybbo_bm5 = readResultsFromJson('trap5','py_bbo','bm5')
    res_cppbbo_bm5 = readResultsFromJson('trap5','cpp_bbo','bm5')
    #results = [("GBO-Python",res_py_bm5),("GBO-C++",res_cpp_bm5),("BBO-Python",res_pybbo_bm5),("BBO-C++",res_cppbbo_bm5)]
    #xticks = [10*2**i for i in range(20)]
    xticks = []
    #plotTimePerFEV(results,out_dir=out_dir,tag='trap5_bm5',max_dim=10000,xticks=xticks)
    #plotEvalTimeScalability(results,out_dir=out_dir,tag='trap5_bm5',max_dim=10000,xticks=xticks)
    #plotTimeScalability(results,out_dir=out_dir,tag='trap5_bm5',max_dim=10000,xticks=xticks)
    #plotFEVScalability(results,out_dir=out_dir,tag='trap5_bm5',max_dim=10000,xticks=xticks)

    res_py_lt = readResultsFromJson('trap5','py','lt')
    res_py_lt100 = readResultsFromJson('trap5','py','lt100')
    res_cy_lt100 = readResultsFromJson('trap5','cy','lt100')
    #res_cyc = readResultsFromJson('trap5','cyc','lt')
    res_cpp_lt = readResultsFromJson('trap5','cpp','lt')
    res_cpp_lt50 = readResultsFromJson('trap5','cpp','lt50')
    res_cpp_lt100 = readResultsFromJson('trap5','cpp','lt100')
    res_pybbo_lt100 = readResultsFromJson('trap5','py_bbo','lt100')
    res_cppbbo_lt = readResultsFromJson('trap5','cpp_bbo','lt')
    res_cppbbo_lt50 = readResultsFromJson('trap5','cpp_bbo','lt50')
    res_cppbbo_lt100 = readResultsFromJson('trap5','cpp_bbo','lt100')
    #results = [("Python",res_py_lt100),("C++",res_cpp_lt100),("Python(BBO)",res_pybbo_lt100),("C++(BBO)",res_cppbbo_lt100)]
    #plotTimePerFEV(results,out_dir=out_dir,tag='trap5_lt100',max_dim=10000,xticks=xticks)
    #plotTimeScalability(results,out_dir=out_dir,tag='trap5_lt100',max_dim=10000,xticks=xticks)
    #plotEvalTimeScalability(results,out_dir=out_dir,tag='trap5_lt100',max_dim=10000,xticks=xticks)
    #plotFEVScalability(results,out_dir=out_dir,tag='trap5_lt100',max_dim=10000,xticks=xticks)
    
    res_py_slt = readResultsFromJson('trap5','py','slt')
    #res_cy_slt = readResultsFromJson('trap5','cy','slt')
    res_cpp_slt = readResultsFromJson('trap5','cpp','slt')
    res_pybbo_lt = readResultsFromJson('trap5','py_bbo','lt')
    res_cppbbo_lt = readResultsFromJson('trap5','cpp_bbo','lt')
    results = [("GBO-Python",res_py_lt),("GBO-C++",res_cpp_lt),("BBO-Python",res_pybbo_lt),("BBO-C++",res_cppbbo_lt)]
    #plotTimePerFEV(results,out_dir=out_dir,tag='trap5_slt',max_dim=10000,xticks=xticks)
    plotEvalTimeScalability(results,out_dir=out_dir,tag='trap5_lt',max_dim=10000,xticks=xticks)
    #plotTimeScalability(results,out_dir=out_dir,tag='trap5_slt',max_dim=10000,xticks=xticks)
    #plotFEVScalability(results,out_dir=out_dir,tag='trap5_slt',max_dim=10000,xticks=xticks)
    
    results = [("GBO-LT",res_cpp_lt),("GBO-SLT",res_cpp_slt),("BBO-LT",res_cppbbo_lt)]
    #results = [("BBO-LT",res_cppbbo_lt),("BBO-LT(50)",res_cppbbo_lt50),("BBO-LT(100)",res_cppbbo_lt100),("GBO-SLT",res_cpp_slt),("GBO-LT",res_cpp_lt),("GBO-LT(50)",res_cpp_lt50),("GBO-LT(100)",res_cpp_lt100)]
    #plotEvalTimeScalability(results,out_dir=out_dir,tag='trap5_cpp_lms',max_dim=10000,xticks=xticks)
    #plotTimeScalability(results,out_dir=out_dir,tag='trap5_cpp_lms',max_dim=10000,xticks=xticks)
    plotFEVScalability(results,out_dir=out_dir,tag='trap5_cpp_lms',max_dim=10000,xticks=xticks)
    
    results = [("GBO-LT",res_py_lt),("GBO-SLT",res_py_slt),("BBO-LT",res_pybbo_lt)]
    #results = [("BBO-LT",res_cppbbo_lt),("BBO-LT(50)",res_cppbbo_lt50),("BBO-LT(100)",res_cppbbo_lt100),("GBO-SLT",res_cpp_slt),("GBO-LT",res_cpp_lt),("GBO-LT(50)",res_cpp_lt50),("GBO-LT(100)",res_cpp_lt100)]
    #plotEvalTimeScalability(results,out_dir=out_dir,tag='trap5_py_lms',max_dim=10000,xticks=xticks)
    plotTimeScalability(results,out_dir=out_dir,tag='trap5_py_lms',max_dim=10000,xticks=xticks)
    #plotFEVScalability(results,out_dir=out_dir,tag='trap5_py_lms',max_dim=10000,xticks=xticks)

if __name__ == '__main__':
    #runExperiments()
    makePlots("plots/trap5/")