#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 14:12:32 2017

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
    net = values[1]
    
         
    
    # store Data
    timeSeries[net] = np.loadtxt(sys.argv[i]) 

ax = plt.subplot(111)
ax.spines["top"].set_visible(False) 
ax.spines["right"].set_visible(False)  
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

#plot results in detail 
for key in timeSeries:
    #labelString = "p="+str(key[0])+"; q="+str(key[1])
    labelString = str(key)
    plt.plot(timeSeries[key][:,0],timeSeries[key][:,1],label=labelString,alpha=0.5,lw=1)
    
    print key, np.mean(timeSeries[key][:,1]), np.std(timeSeries[key][:,1])
    
legend=plt.legend(numpoints=1,frameon=False,ncol=3,loc='upper center', bbox_to_anchor=(0.5, -0.25))
legend.get_title().set_fontsize('21') 
plt.setp(plt.gca().get_legend().get_texts(), fontsize='26')

plt.xlabel("Node number")
plt.ylabel("OutDegree-InDegree")

plt.savefig(str(plotname)+'_Detail.pdf', bbox_inches='tight')
plt.show()

#plot results in histogram
for key in timeSeries:
    
    labelString = str(key)
    plt.hist(timeSeries[key][:,1],bins=50,normed=1,label=labelString,alpha=0.5,lw=1)  
    plt.yscale('log', nonposy='clip')
    
legend=plt.legend(numpoints=1,frameon=False,ncol=3,loc='upper center', bbox_to_anchor=(0.5, -0.25))
legend.get_title().set_fontsize('21') 
plt.setp(plt.gca().get_legend().get_texts(), fontsize='26')

plt.xlabel("OutDegree-InDegree")
plt.ylabel("PMF")

plt.savefig(str(plotname)+'_Histo.pdf', bbox_inches='tight')
plt.show()