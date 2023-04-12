from experiments import runExperiment
from output import readResultsFromJson
from plotting import plotEvalTimeScalability, plotTimePerFEV, plotTimeScalability, plotFEVScalability

def runExperimentsMaxCutBFLT(bound,max_dim=-1,rerun=False,nthreads=1,nruns=10,max_sec=-1,input_file="",vtr_file=""):
    runExperiment('maxcut','cy','lt'+str(bound),max_dim=-1,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec,input_file=input_file,vtr_file=vtr_file)
    #runExperiment('maxcut','cyc','lt'+str(bound),max_dim=-1,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec,input_file=input_file,vtr_file=vtr_file)
    runExperiment('maxcut','cpp','lt'+str(bound),max_dim=-1,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec,input_file=input_file,vtr_file=vtr_file)
    runExperiment('maxcut','cpp_bbo','lt'+str(bound),max_dim=-1,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec,input_file=input_file,vtr_file=vtr_file)
    runExperiment('maxcut','py','lt'+str(bound),max_dim=-1,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec,input_file=input_file,vtr_file=vtr_file)
    runExperiment('maxcut','py_bbo','lt'+str(bound),max_dim=-1,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec,input_file=input_file,vtr_file=vtr_file)

def runExperimentsMaxCutLT(max_dim=-1,rerun=False,nthreads=1,nruns=10,max_sec=-1,input_file="",vtr_file=""):
    runExperiment('maxcut','cy','lt',max_dim=-1,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec,input_file=input_file,vtr_file=vtr_file)
    #runExperiment('maxcut','cyc','lt',max_dim=-1,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec,input_file=input_file,vtr_file=vtr_file)
    runExperiment('maxcut','cpp','lt',max_dim=-1,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec,input_file=input_file,vtr_file=vtr_file)
    runExperiment('maxcut','cpp_bbo','lt',max_dim=-1,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec,input_file=input_file,vtr_file=vtr_file)
    runExperiment('maxcut','py','lt',max_dim=-1,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec,input_file=input_file,vtr_file=vtr_file)
    runExperiment('maxcut','py_bbo','lt',max_dim=-1,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec,input_file=input_file,vtr_file=vtr_file)

def runExperimentsMaxCutStaticLT(max_dim=-1,rerun=False,nthreads=1,nruns=10,max_sec=-1,input_file="",vtr_file=""):
    runExperiment('maxcut','py','slt',max_dim=-1,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec,input_file=input_file,vtr_file=vtr_file)
    runExperiment('maxcut','cy','slt',max_dim=-1,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec,input_file=input_file,vtr_file=vtr_file)
    #runExperiment('maxcut','cyc','slt',max_dim=-1,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec,input_file=input_file,vtr_file=vtr_file)
    runExperiment('maxcut','cpp','slt',max_dim=-1,rerun=rerun,nthreads=nthreads,nruns=nruns,max_sec=max_sec,input_file=input_file,vtr_file=vtr_file)

def runExperiments(input_dir):
    rerun_experiments = False 
    nthreads = 15
    nruns = 30
    max_sec = 600
    dims = [9,16,25,49,100,196,400,784,1600]
    for d in dims:
        input_file = f'{input_dir}n{d:07d}i00.txt'
        vtr_file = f'{input_dir}n{d:07d}i00.bkv'
        runExperimentsMaxCutLT(max_dim=-1,rerun=rerun_experiments,nthreads=nthreads,nruns=nruns,input_file=input_file,vtr_file=vtr_file)
        runExperimentsMaxCutStaticLT(max_dim=-1,rerun=rerun_experiments,nthreads=nthreads,nruns=nruns,input_file=input_file,vtr_file=vtr_file)
        #runExperimentsMaxCutBFLT(10,max_dim=-1,rerun=rerun_experiments,nthreads=nthreads,nruns=nruns,input_file=input_file,vtr_file=vtr_file)
        runExperimentsMaxCutBFLT(50,max_dim=-1,rerun=rerun_experiments,nthreads=nthreads,nruns=nruns,input_file=input_file,vtr_file=vtr_file)
        runExperimentsMaxCutBFLT(100,max_dim=-1,rerun=rerun_experiments,nthreads=nthreads,nruns=nruns,input_file=input_file,vtr_file=vtr_file)

def makePlots(out_dir="plots/"):
    max_dim = 2500
    res_py_lt = readResultsFromJson('maxcut','py','lt')
    res_cy_lt = readResultsFromJson('maxcut','cy','lt')
    #res_cyc = readResultsFromJson('maxcut_cyc_lt')
    res_cpp_lt = readResultsFromJson('maxcut','cpp','lt')
    res_pybbo_lt = readResultsFromJson('maxcut','py_bbo','lt')
    res_cppbbo_lt = readResultsFromJson('maxcut','cpp_bbo','lt')
    results = [("Python",res_py_lt),("Cython",res_cy_lt),("C++",res_cpp_lt),("Python(BBO)",res_pybbo_lt),("C++(BBO)",res_cppbbo_lt)]
    plotTimePerFEV(results,out_dir=out_dir,tag='maxcut_lt',max_dim=max_dim)
    plotTimeScalability(results,out_dir=out_dir,tag='maxcut_lt',max_dim=max_dim)
    plotEvalTimeScalability(results,out_dir=out_dir,tag='maxcut_lt',max_dim=max_dim)
    plotFEVScalability(results,out_dir=out_dir,tag='maxcut_lt',max_dim=max_dim)
    
    res_py_slt = readResultsFromJson('maxcut','py','slt')
    res_cy_slt = readResultsFromJson('maxcut','cy','slt')
    res_cpp_slt = readResultsFromJson('maxcut','cpp','slt')
    results = [("Python",res_py_slt),("Cython",res_cy_slt),("C++",res_cpp_slt),("Python(BBO)",res_pybbo_lt),("C++(BBO)",res_cppbbo_lt)]
    plotTimePerFEV(results,out_dir=out_dir,tag='maxcut_slt',max_dim=max_dim)
    plotEvalTimeScalability(results,out_dir=out_dir,tag='maxcut_slt',max_dim=max_dim)
    plotTimeScalability(results,out_dir=out_dir,tag='maxcut_slt',max_dim=max_dim)
    plotFEVScalability(results,out_dir=out_dir,tag='maxcut_slt',max_dim=max_dim)
    
    res_py_lt50 = readResultsFromJson('maxcut','py','lt50')
    res_py_lt100 = readResultsFromJson('maxcut','py','lt100')
    res_pybbo_lt = readResultsFromJson('maxcut','py_bbo','lt')
    res_pybbo_lt50 = readResultsFromJson('maxcut','py_bbo','lt50')
    res_pybbo_lt100 = readResultsFromJson('maxcut','py_bbo','lt100')
    results = [("BBO-LT",res_pybbo_lt),("BBO-LT(50)",res_pybbo_lt50),("BBO-LT(100)",res_pybbo_lt100),("GBO-SLT",res_py_slt),("GBO-LT",res_py_lt),("GBO-LT(50)",res_py_lt50),("GBO-LT(100)",res_py_lt100)]
    plotEvalTimeScalability(results,out_dir=out_dir,tag='maxcut_py_lms',max_dim=max_dim)
    plotTimeScalability(results,out_dir=out_dir,tag='maxcut_py_lms',max_dim=max_dim)
    plotFEVScalability(results,out_dir=out_dir,tag='maxcut_py_lms',max_dim=max_dim)
    
    res_cpp_lt50 = readResultsFromJson('maxcut','cpp','lt50')
    res_cpp_lt100 = readResultsFromJson('maxcut','cpp','lt100')
    res_cppbbo_lt = readResultsFromJson('maxcut','cpp_bbo','lt')
    res_cppbbo_lt50 = readResultsFromJson('maxcut','cpp_bbo','lt50')
    res_cppbbo_lt100 = readResultsFromJson('maxcut','cpp_bbo','lt100')
    results = [("BBO-LT",res_cppbbo_lt),("BBO-LT(50)",res_cppbbo_lt50),("BBO-LT(100)",res_cppbbo_lt100),("GBO-SLT",res_cpp_slt),("GBO-LT",res_cpp_lt),("GBO-LT(50)",res_cpp_lt50),("GBO-LT(100)",res_cpp_lt100)]
    plotEvalTimeScalability(results,out_dir=out_dir,tag='maxcut_cpp_lms',max_dim=max_dim)
    plotTimeScalability(results,out_dir=out_dir,tag='maxcut_cpp_lms',max_dim=max_dim)
    plotFEVScalability(results,out_dir=out_dir,tag='maxcut_cpp_lms',max_dim=max_dim)

if __name__ == '__main__':
    runExperiments("problem_data/maxcut/set0c/")
    makePlots("plots/set0c/")