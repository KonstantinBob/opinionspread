#!/usr/bin/env python

# imports
import numpy as np
import sys

# ========= constants
errormessage = "usage: "+sys.argv[0]+" inputfile"

# global variables 
argc = len(sys.argv)

# ========= functions

# ========= main =======================
if argc == 2:
    
    # get the names
    filename = sys.argv[1]
       
    # get the data
    left,right = np.loadtxt(filename,unpack=True,comments="#")
    
    data = np.append(left,right)
    
    # get the node ids
    ids = np.unique(data)
    size = len(ids)
    
    # make a map
    mapIds = dict({})
    
    for index in np.arange(size):
        mapIds[ids[index]] = index
          
    # reprint the data
    for index in np.arange(len(left)):
        print mapIds[left[index]], mapIds[right[index]]
    
else:
    print errormessage