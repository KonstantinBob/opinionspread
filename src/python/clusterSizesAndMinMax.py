#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 13:59:19 2017

@author: konbob
"""
import numpy as np
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



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
    
# get data
for i in np.arange(2,argc):
    
    # get the values of q and e 
    values = sys.argv[i].split("_")
    net = values[0]
    e = float(values[1])
    q = float(values[2])
        
    
    # store Data
    timeSeries[(e,q,net)] = np.loadtxt(sys.argv[i]) 

#plot results
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

X = np.array([]) 
Y = np.array([]) 
Z = np.array([]) 
Colors =  np.array([]) 

for key in timeSeries:
    X = np.append(X,key[0])
    Y = np.append(Y,key[1])
    Z = np.append(Z,np.max(timeSeries[key])-np.min(timeSeries[key]))
    Colors = np.append(Colors,key[0])
    
    
   
ax.scatter(X,Y,Z,c=Colors,cmap=plt.cm.viridis,s=10)


ax.set_xlabel('\n p')
ax.set_ylabel('\n q')
ax.set_zlabel('\n Max Dev')


plt.title(str(net))
plt.show()




























