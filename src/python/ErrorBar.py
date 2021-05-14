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
if argc == 2:
    filename = sys.argv[1]
    
    # get the legend
    file = open(filename)
    title = file.readline().replace('\n', '')
    xLabel = file.readline().replace('\n', '')
    yLabel = file.readline().replace('\n', '')
    plotname = file.readline().replace('\n', '')
    file.close()
    
    # get the data
    x,y,deltaY = np.loadtxt(filename,unpack=True,skiprows=4,comments="#")
    
    #perform the plot
    plt.title(title)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    
    plt.xlim([0,3.5e3])
    plt.ylim([0,7])
    #plt.ticklabel_format(style='sci',axis='x',scilimits=(0,0))
    
    plt.errorbar(x,y,yerr=deltaY, fmt='o',alpha=0.9,markersize = 0.3)
   
    plt.tight_layout()
    plt.savefig(plotname)
    plt.show()
    
else:
    print errormessage