#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 15:28:20 2017

@author: konbob
"""

import numpy as np
import sys
import matplotlib.pyplot as plt

# ========= constants =
errormessage = "usage: python "+sys.argv[0] + "datafile plotname"

# global variables 
argc = len(sys.argv)

# ========= main =======================
# parse the arguments
if argc != 3:    
    print errormessage
    exit
    
# get plotname
plotname =  sys.argv[2]
    
# store Data
data = np.loadtxt(sys.argv[1]) 

# plot the data
size = np.size(data,axis=1)

#for index in np.arange(size):
#    plt.plot(data[:,index])
 
for index in [0,13]:
    plt.plot(data[:,index],label=str(index))

plt.legend(loc=0)
plt.show()