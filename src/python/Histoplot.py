#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

# ========= constants
errormessage = "usage: "+sys.argv[0]+" inputfile"

# global variables 
argc = len(sys.argv)

# ========= functions
def pot(x,a,b,c,d,e):    
    return a*np.power(x,4.)+b*np.power(x,3.)+c*np.power(x,2.)+d*x+e
    
def theory(x):
    return np.exp(-x)
# ========= main =======================
if argc == 2:
    # function
    
    
    # Setup Latex for Plots
    plt.rcParams.update({'text.usetex': True})
    plt.rcParams.update({'font.family': 'lmodern'})
    plt.rcParams.update({'font.size': 23})
    plt.rcParams.update({'backend': 'PDF'})
    
    filename = sys.argv[1]
    
    # get the legend
    file = open(filename)
    title = file.readline().replace('\n', '')
    xLabel = file.readline().replace('\n', '')
    yLabel = file.readline().replace('\n', '')
    plotname = file.readline().replace('\n', '')
    file.close()
    
    # get the data
    x,y = np.loadtxt(filename,unpack=True,skiprows=4)
    
    # calculate the bin with
    binWidth = x[1]-x[0]
    
    #normalization
    normalize = False
    if normalize:
        s = np.sum(y) * binWidth
        y = y / s
        yLabel = "normalized "+ yLabel
            
    #====== perform the plot =======
        
    # for two axes
    fig, ax1 = plt.subplots()
    
    # labels        
    plt.title(title+" ($"+filename+"$)")
    ax1.set_xlabel(xLabel)
        
    # histo plot
    ax1.set_ylabel(yLabel, color='r')
    ax1.tick_params('y', colors='r')
    ax1.bar(x,y,width=binWidth,align='center',alpha=0.6,color='r')
    
    # free energy plot
    ax2 = ax1.twinx()
    z = np.linspace(0.0,1.0,100)
    ax2.plot(z,pot(z,160./3.,-104.,188./3.,-12.,1.0),)
    ax2.set_ylabel('free energy', color='b')
    ax2.tick_params('y', colors='b')
    
    # save and show
    plt.tight_layout()
    plt.savefig(plotname+".jpg")
    plt.show()

else:
    print errormessage
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    