#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 11:28:29 2016

@author: konbob
"""

# imports
import numpy as np
import sys
import matplotlib.pyplot as plt

# ========= constants
errormessage = "usage: "+sys.argv[0]+" inputfile"

# global variables 
argc = len(sys.argv)

# ========= functions
def fd(x,a,b):
   
   return (1.0)/(1.0 + np.exp(b*(x+a)))
  
        
def freeEnergy(x):
    xi = 2.*x+0.5;
    return 1. + (np.sin(10*np.pi*xi))/(2.*xi) + np.power(xi-1.0,4.0);
    
def freeEnergyLight(x):
    return np.power(x-0.5,2.)
    
# ========= main =======================
if argc == 1:
    
    # size
    n = 10;    
    # stepsize 
    s = 1./n
    
    # get the raw matrix
    K = np.ones([n,n]) 
    
    # set it up
    for i in np.arange(n):
        sumRow = 0
        for j in np.arange(n):
            # fill values
            fi = freeEnergy(i*s)
            fj = freeEnergy(j*s)
            K[i,j] = fd(fj-fi,0.,1.)            
            sumRow = sumRow + K[i,j]
        #normalization
        if sumRow > 0:
            for j in np.arange(n):
                K[i,j] = K[i,j] / sumRow
      
    # solve eigenproblem
    w,v = np.linalg.eig(K)
    
    print w[0]
    
    distr =  np.real(v[0,:])
    
    print distr
    
    distr = distr / np.sum(distr)
    
    # plotting
    x = s*np.arange(n)
    
    plt.plot(x,distr,label="prob")
    plt.plot(x,freeEnergyLight(x),label="energy")
    plt.legend(loc="best")
    plt.show()
           
else:
    print errormessage