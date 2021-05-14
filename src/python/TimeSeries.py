#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 13:52:13 2017

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
errormessage = "usage: python "+sys.argv[0] + "plotname datafile1  datafile2 datafile3"

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
pSet = set()
qSet = set()    
    
# get data
for i in np.arange(2,argc):
    
    # get the values of q and e 
    values = sys.argv[i].split("_")
    net = values[0]
    e = float(values[1])
    q = float(values[2])   

    # store values for ordering
    pSet.add(e)
    qSet.add(q)      
    
    # store Data
    timeSeries[(e,q,net)] = np.loadtxt(sys.argv[i]) 

ax = plt.subplot(111)
ax.spines["top"].set_visible(False) 
ax.spines["right"].set_visible(False)  
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

# create dictionary for ordering
qList = list(qSet)
pList = list(pSet)
qDict = {}
pDict = {}
for i,item in enumerate(qList):
    qDict[item] = i+1
for i,item in enumerate(pList):
    pDict[item] = i+1
print "p ",pDict,'\nq ',qDict

#plot results
for key in sorted(timeSeries.keys(), key =  lambda x : float(x[1])):
    p,q,net = key
    labelString = "p="+str(key[0])+"; q="+str(key[1])
    #labelString = "q="+str(key[1])    
    orderNumber = qDict[q]
    print labelString, orderNumber
    plt.plot(timeSeries[key],label=labelString,alpha=0.5,lw=1,zorder=orderNumber)
    
plt.ylim([0.0,1.0])


#handles, labels = ax.get_legend_handles_labels()
# sort both labels and handles by labels
#labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
#ax.legend(handles, labels)


legend=plt.legend(numpoints=1,frameon=False,ncol=1,loc='upper center', bbox_to_anchor=(0.5, -0.25))
legend.get_title().set_fontsize('21') 
plt.setp(plt.gca().get_legend().get_texts(), fontsize='26')

plt.xlabel('Time (Cycles)')
plt.ylabel('$\\frac{N_{\\uparrow}}{N}$')
#text = plt.figtext(0.8,0.15,'$p=0.95$',fontsize=12)
#plt.savefig(str(plotname)+'.pdf',bbox_extra_artists=(legend,), bbox_inches='tight')
plt.savefig(str(plotname)+'.pdf', bbox_inches='tight')
plt.show()
