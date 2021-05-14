#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 11:26:39 2017

@author: konbob
"""

import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy import stats
from scipy.signal import find_peaks_cwt
import pandas as pd
import cPickle as pickle
from collections import OrderedDict

# Setup Latex for Plots
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'font.family': 'lmodern'})
plt.rcParams.update({'font.size': 26})
plt.rcParams.update({'backend': 'PDF'})
params = {'text.latex.preamble' : [r'\usepackage{amsmath}',r'\usepackage{amssymb}']}
plt.rcParams.update(params)
  
# ========= constants =
errormessage = "usage: python "+sys.argv[0] + " plotname  debug flag (-d for true, -n for false) simulationFiles"

# global variables 
argc = len(sys.argv)

# ========= functions


# ========= main =======================
# *                                    *
# *                                    *
# ======================================

#  ========= parse the arguments ========= 
if argc < 3:    
    print errormessage
    sys.exit()
    
# store plotname
plotname = str(sys.argv[1])

#debugging
debug = False
if str(sys.argv[2]) == '-d':
    debug = True
elif str(sys.argv[2]) != '-n':
    print 'Could not parse ',str(sys.argv[2])
    print errormessage
    sys.exit()

# dict for results
stdValues = dict()

#  =========  file loop ========= 
for filename in sys.argv[3:]:
           
    # get the values of q and e 
    values = filename.split("_")
    net = str(values[0])
    size = int(values[1])
    p = float(values[2])
    q = float(values[3])  
    
     # print information
    if debug: 
        print "processing: p ",p,"q ",q,"size ",size       
        
    # store Data    
    df = pd.read_csv(filename,dtype='float',comment='#',header=None)    
    timeSeries = df.as_matrix()[:,0]
    
    # plot it
    time = np.arange(500)
    plt.plot(time,timeSeries,label=str(size),alpha = 0.5 + (float(size))/1.6e4)
    
    # store variance
    stdValues[size] = np.std(timeSeries)
#  =========  end file loop =========     
    
#  ========= plotting ===============

#sizeArray = np.fromiter(iter(stdValues.keys()), dtype=float)
#stdArray = np.fromiter(iter(stdValues.values()), dtype=float)

#plt.plot(sizeArray,stdArray,'x')

# Setup Latex for Plots
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'font.family': 'lmodern'})
plt.rcParams.update({'font.size': 26})
plt.rcParams.update({'backend': 'PDF'})
params = {'text.latex.preamble' : [r'\usepackage{amsmath}',r'\usepackage{amssymb}']}
plt.rcParams.update(params)    

# setup axes
ax = plt.subplot(111)
ax.spines["top"].set_visible(False) 
ax.spines["right"].set_visible(False)  
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
    
# information
plt.xlabel("Time (Cycles)")
plt.ylabel('$\\frac{N_{\\uparrow}}{N}$ ')

# legend
handles, labels = plt.gca().get_legend_handles_labels()    
by_label = OrderedDict(sorted(zip(labels, handles), key = lambda x : float(x[0])  )) 
legend=plt.legend(by_label.values(), by_label.keys(),numpoints=1,frameon=False,ncol=4,loc='upper center', bbox_to_anchor=(0.5, -0.25))
legend.get_title().set_fontsize('21') 
plt.setp(plt.gca().get_legend().get_texts(), fontsize='26')

# show it
if debug:                   
    plt.show()
else:
    plt.savefig(plotname+'.pdf',bbox_inches='tight')
    plt.title(net+' with p = '+str(p)+' q = '+str(q))
    plt.savefig(plotname+'.png',bbox_inches='tight')
    plt.close()    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    