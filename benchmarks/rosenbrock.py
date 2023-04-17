from experiments import runExperiment
from output import readResultsFromJson
from plotting import plotEvalTimeScalability, plotTimePerFEV, plotTimeScalability, plotFEVScalability

def runExperimentsRosenbrockUni(max_dim=1000,rerun=False,nthreads=1,nruns=10,max_sec=-1):
    runExperiment('rosenbrock','cpp','uni',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    #runExperiment('rosenbrock','cpp_bbo','uni',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    #runExperiment('rosenbrock','py','uni',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    #runExperiment('rosenbrock','py_bbo','uni',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)

def runExperimentsRosenbrockLT(max_dim=1000,rerun=False,nthreads=1,nruns=10,max_sec=-1):
    runExperiment('rosenbrock','cpp','lt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('rosenbrock','cpp_bbo','lt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('rosenbrock','py','lt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('rosenbrock','py_bbo','lt',max_dim=160,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('rosenbrock','py_bbo','lt50',max_dim=160,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('rosenbrock','cpp_bbo','lt50',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)

def runExperimentsRosenbrockSLT(max_dim=1000,rerun=False,nthreads=1,nruns=10,max_sec=-1):
    #runExperiment('rosenbrock','py','slt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    #runExperiment('rosenbrock','cpp','slt',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    #runExperiment('rosenbrock','cpp','slt50',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperiment('rosenbrock','py','slt50',max_dim=max_dim,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec)

def runExperiments():
    rerun_experiments = False 
    nthreads = 15
    nruns = 30
    max_sec = 3600
    runExperimentsRosenbrockUni(max_dim=1000,rerun=rerun_experiments,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperimentsRosenbrockLT(max_dim=400,rerun=rerun_experiments,nthreads=nthreads,nruns=nruns,max_sec=max_sec)
    runExperimentsRosenbrockSLT(max_dim=1000,rerun=rerun_experiments,nthreads=nthreads,nruns=nruns,max_sec=max_sec)

def makePlots(out_dir="plots/"):
    xticks = []
    yrange = [1e3,1e8]
    #res_py_lt = readResultsFromJson('rosenbrock','py','lt')
    #res_cpp_lt = readResultsFromJson('rosenbrock','cpp','lt')
    #res_cppbbo_lt = readResultsFromJson('rosenbrock','cpp_bbo','lt')
    
    #res_py_slt = readResultsFromJson('rosenbrock','py','slt')
    res_py_slt50 = readResultsFromJson('rosenbrock','py','slt50')
    res_cpp_slt50 = readResultsFromJson('rosenbrock','cpp','slt50')
    #res_cy_slt = readResultsFromJson('rosenbrock','cy','slt')
    #res_cpp_slt = readResultsFromJson('rosenbrock','cpp','slt')
    #res_pybbo_lt = readResultsFromJson('rosenbrock','py_bbo','lt')
    #res_cppbbo_lt = readResultsFromJson('rosenbrock','cpp_bbo','lt')
    res_cppbbo_lt50 = readResultsFromJson('rosenbrock','cpp_bbo','lt50')
    res_pybbo_lt50 = readResultsFromJson('rosenbrock','py_bbo','lt50')
    res_cppbbo_uni = readResultsFromJson('rosenbrock','cpp_bbo','uni')
    res_cpp_uni = readResultsFromJson('rosenbrock','cpp','uni')
    res_py_uni = readResultsFromJson('rosenbrock','py','uni')
    res_pybbo_uni = readResultsFromJson('rosenbrock','py_bbo','uni')
    results = [("GBO-Python",res_py_uni),("GBO-C++",res_cpp_uni),("BBO-Python",res_pybbo_uni),("BBO-C++",res_cppbbo_uni)]
    #plotTimePerFEV(results,out_dir=out_dir,tag='rosenbrock_uni',max_dim=1000,xticks=xticks)
    plotEvalTimeScalability(results,out_dir=out_dir,tag='rosenbrock_uni',max_dim=1000,xticks=xticks)
    #plotTimeScalability(results,out_dir=out_dir,tag='rosenbrock_uni',max_dim=1000,xticks=xticks)
    #plotFEVScalability(results,out_dir=out_dir,tag='rosenbrock_uni',max_dim=1000,xticks=xticks,yrange=yrange)

    results = [("GBO-Uni",res_cpp_uni),("GBO-SLT(50)",res_cpp_slt50),("BBO-Uni",res_cppbbo_uni),("BBO-LT(50)",res_cppbbo_lt50)]
    #results = [("BBO-LT",res_cppbbo_lt),("BBO-LT(50)",res_cppbbo_lt50),("BBO-LT(100)",res_cppbbo_lt100),("GBO-SLT",res_cpp_slt),("GBO-LT",res_cpp_lt),("GBO-LT(50)",res_cpp_lt50),("GBO-LT(100)",res_cpp_lt100)]
    #plotEvalTimeScalability(results,out_dir=out_dir,tag='rosenbrock_cpp_lms',max_dim=1000,xticks=xticks)
    #plotTimeScalability(results,out_dir=out_dir,tag='rosenbrock_cpp_lms',max_dim=1000,xticks=xticks)
    plotFEVScalability(results,out_dir=out_dir,tag='rosenbrock_cpp_lms',max_dim=1000,xticks=xticks,yrange=yrange)
    
    #results = [("BBO-LT",res_pybbo_lt),("GBO-SLT",res_py_slt),("GBO-LT",res_py_lt)]
    results = [("GBO-Uni",res_py_uni),("GBO-SLT(50)",res_py_slt50),("BBO-Uni",res_pybbo_uni),("BBO-LT(50)",res_pybbo_lt50)]
    #results = [("BBO-LT",res_pybbo_lt),("GBO-LT",res_py_lt)]
    #results = [("BBO-LT",res_cppbbo_lt),("BBO-LT(50)",res_cppbbo_lt50),("BBO-LT(100)",res_cppbbo_lt100),("GBO-SLT",res_cpp_slt),("GBO-LT",res_cpp_lt),("GBO-LT(50)",res_cpp_lt50),("GBO-LT(100)",res_cpp_lt100)]
    #plotEvalTimeScalability(results,out_dir=out_dir,tag='rosenbrock_py_lms',max_dim=1000,xticks=xticks)
    plotTimeScalability(results,out_dir=out_dir,tag='rosenbrock_py_lms',max_dim=1000,xticks=xticks)
    #plotFEVScalability(results,out_dir=out_dir,tag='rosenbrock_py_lms',max_dim=1000,xticks=xticks,yrange=yrange)

if __name__ == '__main__':
    #runExperiments()
    makePlots("plots/rosenbrock/")