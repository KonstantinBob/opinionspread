#!/usr/bin/env python

# imports
import numpy as np
import matplotlib.pyplot as plt
import sys

# Setup Latex for Plots
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'font.family': 'lmodern'})
plt.rcParams.update({'font.size': 23})
plt.rcParams.update({'backend': 'PDF'})

# ========= constants
errormessage = "usage: "+sys.argv[0]+" inputfile  plotname doFit"

# global variables 
argc = len(sys.argv)

# ========= functions

# ========= main =======================
if argc == 4:
    
    # get the arguments
    filename = sys.argv[1]
    plotname = sys.argv[2]
    fitFlag = sys.argv[3]
    
    # check if fitting is to be done
    doFit = False  
    if fitFlag == "y":
        doFit = True     
    
    # set the legend    
    title = ""
    xLabel = "Nodes"
    yLabel = "Edges"    
    
    
         
    #perform the plot
    ax = plt.subplot(111)
    plt.title(title)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)     
    
    plt.xlim([1e2,1e5])
    plt.ylim([1e2,1e6])
    plt.ticklabel_format(style='sci',axis='x',scilimits=(0,0))
    plt.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
        
    if doFit:
        # get the data
        nodes, edges, error, fit = np.loadtxt(filename,unpack=True,comments="#")
        plt.loglog(nodes,fit,'g',alpha=0.9,label="Fit",zorder = 2)
    else:
        ax.set_xscale("log")
        ax.set_yscale("log")
        # get the data
        nodes, edges, error = np.loadtxt(filename,unpack=True,comments="#")
    plt.errorbar(nodes,edges,yerr=error ,fmt='o',alpha=0.9,label="Data", zorder = 1)
    
    
    ax.spines["top"].set_visible(False) 
    ax.spines["right"].set_visible(False)  
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    legend=plt.legend(loc='best',numpoints=1,frameon=False)    
    plt.setp(plt.gca().get_legend().get_texts(), fontsize='18') #legend 'list' fontsize
    plt.tight_layout()
    plt.savefig(plotname+".pdf")
    plt.show()
    
else:
    print errormessage