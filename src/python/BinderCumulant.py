#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 15:40:55 2017

@author: konbob
"""
import numpy as np
import sys
import matplotlib.pyplot as plt
import cPickle as pickle
import pymc3 as pm
import pandas as pd
from collections import OrderedDict
import os
from scipy.stats import norm
import warnings
warnings.filterwarnings("ignore")

from multiprocessing import Pool

# Setup Latex for Plots
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'font.family': 'lmodern'})
plt.rcParams.update({'font.size': 26})
plt.rcParams.update({'backend': 'PDF'})

# ========= constants =
errormessage = "usage: python "+sys.argv[0] + " peakfile  debug flag (-d for true, -n for false) plotname"

# global variables 
argc = len(sys.argv)


# ========= functions

#  ========= parse the arguments ========= 
if argc < 4:    
    print errormessage
    sys.exit()
    
# store filename for storage   
filenameStore = str(sys.argv[1])+"Cumulants"

#debugging
debug = False
if str(sys.argv[2]) == '-d':
    debug = True
elif str(sys.argv[2]) != '-n':
    print errormessage
    sys.exit()
    
plotname = str(sys.argv[3])
    
# ============= MAIN =================
    
# load file
cumulants = {}
if os.path.isfile(filenameStore):   
    fileCumu = open(filenameStore,"rb")
    cumulants = pickle.load(fileCumu)    
    fileCumu.close()

ax = plt.subplot(111)
ax.spines["top"].set_visible(False) 
ax.spines["right"].set_visible(False)  
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

# for colors
cmap = plt.cm.get_cmap('plasma')

# find max size in keys
sizeMax = float(max(cumulants,key= lambda x: x[2])[2])

# find qmin and qmax
qMin = float(min(cumulants,key= lambda x: x[1])[1])
qMax = float(max(cumulants,key= lambda x: x[1])[1])

# plot them
for key in cumulants:
    p,q,size = key    
    cumu = cumulants[key]
    
    color = cmap(float(size)/sizeMax)                
    plt.plot(q,cumu,'x',color=color,label=str(size))


plt.xlabel('q')
plt.ylabel('$\\frac{M_{4}}{M_{2}^{2}}$')   

# remove duplicates in legend
handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
legend=plt.legend(by_label.values(), by_label.keys(),numpoints=1,frameon=False,ncol=3,loc='upper center', bbox_to_anchor=(0.5, -0.25))
legend.get_title().set_fontsize('21') 
plt.setp(plt.gca().get_legend().get_texts(), fontsize='26')

plt.savefig(str(plotname)+'.pdf', bbox_inches='tight')   
plt.show()