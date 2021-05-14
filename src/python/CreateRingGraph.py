#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 12:14:56 2017

@author: konbob
"""

"""
Script for creating a ring graph in SNAP file format.
"""

import numpy as np
import sys


# ========= constants
errormessage = "usage: python "+sys.argv[0] + "number of nodes"

# global variables 
argc = len(sys.argv)


# ========= functions
            
# ========= main =======================
# parse the arguments
if argc != 2:    
    print errormessage
    
# parameters
numberNodes=int(sys.argv[1])
        
#open the file
outFile=file("network.txt",'w')   
    
# write header
outFile.write("# Directed graph \n# Ring graph\n")

# write nodes and edges
outFile.write("# Nodes: " +  str(numberNodes)+" Edges: "+ str(numberNodes)+"\n# FromNodeId	ToNodeId\n")

# write edges to file
for i in np.arange(numberNodes):
    outFile.write(str(i)+" "+str((i+1)%numberNodes)+"\n")    

#close file
outFile.close()        