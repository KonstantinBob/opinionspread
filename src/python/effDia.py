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
errormessage = "usage: "+sys.argv[0]+" inputfile  plotname"

# global variables 
argc = len(sys.argv)

# ========= functions

# ========= main =======================
if argc == 3:
    
    # get the names
    filename = sys.argv[1]
    plotname = sys.argv[2]
    
    # set the legend    
    title = ""
    xLabel = "Nodes"
    yLabel = "Effective diameter"    
    
    # get the data
    size,dia,deltaDia,var,DeltaVar = np.loadtxt(filename,unpack=True,comments="#")
    
    # calculate the standard deviation and its error
    std = np.sqrt(var)
    errorStd = DeltaVar/(2.0*std)
    
    #perform the plot
    plt.title(title)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    
    
    plt.xlim([0,3.05e4])
    plt.ylim([0,10])
    plt.ticklabel_format(style='sci',axis='x',scilimits=(0,0))
    
    plt.errorbar(size,dia,yerr=deltaDia,fmt='o',alpha=1.0)
    #plt.errorbar(size,dia,yerr=errorStd,fmt='o',alpha=0.9,label="std")
    
    ax = plt.subplot(111)
    ax.spines["top"].set_visible(False) 
    ax.spines["right"].set_visible(False)  
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')    
    plt.tight_layout()
    
    plt.savefig(plotname+".pdf")
    plt.show()
    
else:
    print errormessage