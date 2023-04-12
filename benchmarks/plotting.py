import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

def addResultToPlot(results,succ_rate_results,label):
    dimensions = np.sort(np.int_(list(results.keys())))
    dimensions = [d for d in dimensions if succ_rate_results[d] > 0.99]
    medians = [np.median(results[dim]) for dim in dimensions]
    err_high = [np.percentile(results[dim],90)-np.median(results[dim]) for dim in dimensions] 
    err_low = [np.median(results[dim])-np.percentile(results[dim],10) for dim in dimensions] 
    plt.errorbar(dimensions,medians,yerr=[err_low,err_high],capsize=3,label=label)

def addTimePerEvalResultToPlot(time_results,fev_results,succ_rate_results,label):
    dimensions = np.sort(list(time_results.keys()))
    dimensions = [d for d in dimensions if succ_rate_results[d] > 0.99]
    tpf_results = {dim:np.divide(time_results[dim],fev_results[dim]) for dim in dimensions}
    medians = [np.median(tpf_results[dim]) for dim in dimensions]
    err_high = [np.percentile(tpf_results[dim],90)-np.median(tpf_results[dim]) for dim in dimensions] 
    err_low = [np.median(tpf_results[dim])-np.percentile(tpf_results[dim],10) for dim in dimensions] 
    plt.errorbar(dimensions,medians,yerr=[err_low,err_high],capsize=3,label=label)

def plotFEVScalability(results,out_dir="plots/",tag='',min_dim=8,max_dim=1000,xticks=[]):
    fig,ax = plt.subplots(dpi=150)
    for label,res in results:
        addResultToPlot(res['fev'],res['succ_rate'],label=label)
    fig.patch.set_facecolor('white')
    plt.yscale('log')
    plt.xscale('log',base=10)
    plt.xlabel("Dimensionality")
    plt.ylabel("Number of evaluations")
    if len(xticks)>0:
        plt.xticks(xticks)
    #plt.ylim(100,50000)
    plt.xlim(min_dim,max_dim)
    plt.grid(linewidth=0.5)
    formatter = plt.ScalarFormatter()
    formatter.set_scientific(False)
    ax.xaxis.set_major_formatter(formatter)
    ax.legend(loc='upper left')
    plt.savefig(out_dir+'/fev_'+tag+'.pdf')

def plotEvalTimeScalability(results,out_dir="plots/",tag='',min_dim=8,max_dim=1000,xticks=[]):
    fig,ax = plt.subplots(dpi=150)
    for label,res in results:
        addResultToPlot(res['eval_t'],res['succ_rate'],label=label)
    fig.patch.set_facecolor('white')
    plt.yscale('log')
    plt.xscale('log',base=10)
    plt.xlabel("Dimensionality")
    plt.ylabel("Evaluation time (s)")
    if len(xticks)>0:
        plt.xticks(xticks)
    #plt.ylim(0.0005,30)
    plt.xlim(min_dim,max_dim)
    plt.grid(linewidth=0.5)
    formatter = plt.ScalarFormatter()
    formatter.set_scientific(False)
    ax.xaxis.set_major_formatter(formatter)
    ax.legend(loc='upper left')
    plt.savefig(out_dir+'/evalt_'+tag+'.pdf')

def plotTimeScalability(results,out_dir="plots/",tag='',min_dim=8,max_dim=1000,xticks=[]):
    fig,ax = plt.subplots(dpi=150)
    for label,res in results:
        addResultToPlot(res['t'],res['succ_rate'],label=label)
    fig.patch.set_facecolor('white')
    plt.yscale('log')
    plt.xscale('log',base=10)
    plt.xlabel("Dimensionality")
    plt.ylabel("Time (s)")
    if len(xticks)>0:
        plt.xticks(xticks)
    #plt.ylim(0.0005,30)
    plt.xlim(min_dim,max_dim)
    plt.grid(linewidth=0.5)
    formatter = plt.ScalarFormatter()
    formatter.set_scientific(False)
    ax.xaxis.set_major_formatter(formatter)
    ax.legend(loc='upper left')
    plt.savefig(out_dir+'/t_'+tag+'.pdf')

def plotTimePerFEV(results,out_dir="plots/",tag='',min_dim=8,max_dim=1000,xticks=[]):
    fig,ax = plt.subplots(dpi=150)
    for label,res in results:
        #print(label,res)
        addTimePerEvalResultToPlot(res['t'],res['fev'],res['succ_rate'],label=label)
    fig.patch.set_facecolor('white')
    plt.yscale('log')
    plt.xscale('log',base=10)
    plt.xlabel("Dimensionality")
    plt.ylabel("Time (s) / Evaluations ")
    if len(xticks)>0:
        plt.xticks(xticks)
    #plt.ylim(0.0005,30)
    plt.xlim(min_dim,max_dim)
    plt.grid(linewidth=0.5)
    formatter = plt.ScalarFormatter()
    formatter.set_scientific(False)
    ax.xaxis.set_major_formatter(formatter)
    ax.legend(loc='upper left')
    plt.savefig(out_dir+'/tperfev_'+tag+'.pdf')
