#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 14:47:30 2017

@author: konbob
"""

import numpy as np
import sys
import matplotlib.pyplot as plt

# Setup Latex for Plots
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'font.family': 'lmodern'})
plt.rcParams.update({'font.size': 26})
plt.rcParams.update({'backend': 'PDF'})

#plt.axes().set_aspect('equal', 'datalim')   

# ========= constants =
errormessage = "usage: python "+sys.argv[0] + " plotname datafilePrefix"

# global variables 
argc = len(sys.argv)


# ========= functions

            
# ========= main =======================
# parse the arguments
if argc < 2:    
    print errormessage
    exit
    
plotname =sys.argv[1]
        
timeSeries = {}    
    
# get data
for i in np.arange(2,argc):
    
    # get the values of q and e 
    values = sys.argv[i].split("_")
    net = values[0]
    timeStep = int(values[1])
         
    
    # store Data
    timeSeries[timeStep] = np.loadtxt(sys.argv[i]) 

ax = plt.subplot(111)
ax.spines["top"].set_visible(False) 
ax.spines["right"].set_visible(False)  
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

#plot results
for key in timeSeries:
    #labelString = "p="+str(key[0])+"; q="+str(key[1])
    labelString = str(key[1])
    plt.plot(timeSeries[key],label=labelString,alpha=0.5,lw=1)
    

plt.savefig(str(plotname)+'.pdf', bbox_inches='tight')
plt.show()