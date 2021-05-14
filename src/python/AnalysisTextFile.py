#!/usr/bin/env python


# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 12:05:13 2016

@author: konbob

"""
import numpy as np
import matplotlib.pyplot as plt
import sys

# ========= constants
errormessage = "usage: "+sys.argv[0]+" inputfile"

# global variables 
argc = len(sys.argv)

# ========= functions

# ========= main =======================
if argc == 2:

    # read the data
    filename = sys.argv[1]
    
    passes, clusters,fitness,op = np.loadtxt(filename,unpack=True)
    
    # for two axes
    fig, ax1 = plt.subplots()         
       
    # free energy plot
    ax2 = ax1.twinx()
    ax2.set_ylim([0,1])
    ax2.plot( passes,op,'o-')
    ax2.set_ylabel('op', color='b')
    ax2.tick_params('y', colors='b')
    
    # cluster plot    
    ax1.set_ylabel("clusters", color='r')
    ax1.tick_params('y', colors='r')
    ax1.plot(passes, clusters,alpha=0.6,color='r')    
            
    plt.title("Plot of "+filename)   
    plt.savefig(filename+".png")
    plt.show()
    
else:
    print errormessage