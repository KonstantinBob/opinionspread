#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 10:35:37 2017

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
errormessage = "usage: python "+sys.argv[0] + " plotname orderParameter datafiles"

# global variables 
argc = len(sys.argv)


# ========= functions

            
# ========= main =======================
# parse the arguments
if argc < 3:    
    print errormessage
    sys.exit()
    
plotname =sys.argv[1]
orderByP = True 

if str(sys.argv[2]) == "q":
    orderByP = False
        
timeSeries = {}    
pSet = set()
qSet = set()    

# get data
for i in np.arange(3,argc):
    
    # get the values of q and e 
    values = sys.argv[i].split("_")
    net = values[0]
    p = float(values[1])
    q = float(values[2])  

    # store values for ordering
    pSet.add(p)
    qSet.add(q)      
    
    # store Data
    timeSeries[(p,q,net)] = np.loadtxt(sys.argv[i]) 
    
    
# ============= plotting ===============

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
for key in sorted(timeSeries.keys()):
    p,q,net = key
    labelString = ""
    orderNumber = 0
    if orderByP:        
        labelString = str(p)
        orderNumber = pDict[p]
    else:
        labelString = str(q)
        orderNumber = qDict[q]
    
    print 'p ',p,' q ',q,' label ',labelString, ' orderNumber ', orderNumber
    
    plt.plot(timeSeries[key],label=labelString,alpha=0.5,lw=1,zorder=orderNumber)
    
plt.ylim([0.0,1.0])


legend=plt.legend(numpoints=1,frameon=False,ncol=3,loc='upper center', bbox_to_anchor=(0.5, -0.25))
legend.get_title().set_fontsize('21') 
plt.setp(plt.gca().get_legend().get_texts(), fontsize='26')

plt.xlabel('Time (Cycles)')
plt.ylabel('$\\frac{N_{\\uparrow}}{N}$')
#text = plt.figtext(0.8,0.15,'$p=0.95$',fontsize=12)
#plt.savefig(str(plotname)+'.pdf',bbox_extra_artists=(legend,), bbox_inches='tight')
plt.savefig(str(plotname)+'.pdf', bbox_inches='tight')
plt.show()
