#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 14:30:36 2017

@author: konbob
"""

"""
Script for exploring the parameter space.
Requires naming scheme Prefix_EXP_Q_.log
"""

import numpy as np
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



# Setup Latex for Plots
#plt.rcParams.update({'text.usetex': True})
#plt.rcParams.update({'font.family': 'lmodern'})
#plt.rcParams.update({'font.size': 26})
#plt.rcParams.update({'backend': 'PDF'})

#plt.axes().set_aspect('equal', 'datalim')   

# ========= constants =
errormessage = "usage: python "+sys.argv[0] + " datafile1 datafile2 ...."

# global variables 
argc = len(sys.argv)


# ========= functions


            
# ========= main =======================
# parse the arguments
if argc < 2:    
    print errormessage
    exit
       
    
# build up storage for data
timeSeries ={} 

# get data
for i in np.arange(1,argc):
    
    # get the values of q and e 
    values = sys.argv[i].split("_")
    e = float(values[1])
    q = float(values[2])
        
    
    # store Data
    timeSeries[(e,q)] =np.loadtxt(sys.argv[i]) 



#print key[0],key[1], np.mean(timeSeries[key]), np.std(timeSeries[key])


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
    Z = np.append(Z,np.std(timeSeries[key]))
    Colors = np.append(Colors,key[0])
    
    
   
ax.scatter(X,Y,Z,c=Colors,cmap=plt.cm.viridis,s=10)


ax.set_xlabel("\n MutStrength")
ax.set_ylabel("\n q")
ax.set_zlabel("\n Fluctuations")

plt.show()























    #data = np.loadtxt(sys.argv[i])
    #plt.plot(data,'-o',label=sys.argv[i])

#plt.legend(loc='best')
#plt.show()