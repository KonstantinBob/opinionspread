#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 11:05:05 2017

@author: konbob
"""

"""
Script for simulating quick and dirty metropolis.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.integrate as integrate

from matplotlib.ticker import NullFormatter

# ========= constants
errormessage = "usage: python "+sys.argv[0] + " filename"

# global variables 
argc = len(sys.argv)


# ========= functions    
def pot(x,a,b,c,d,e):    
    return a*pow(x,4.)+b*pow(x,3.)+c*pow(x,2.)+d*x+e
    
def u(x):
    return pot(x,160./3.,-104.,188./3.,-12.,1.0) 
    
def phi(beta):
    return integrate.quad(lambda x:np.exp(-beta*u(x)) ,0.,1.)[0]
    
def theory(betaEst,x):
    return np.exp(-betaEst*u(x))/phi(betaEst)

# ========= main =======================
# parse the arguments
if argc == 2:    

    # read the filename 
    filename = sys.argv[1]
    
    # set up coefficients aray
    X = np.random.uniform(low=0.0,high=1.0,size=400)
    
    # constants
    strength = 1e-1
    beta = 5.
    
    # time array
    time = np.arange(1e3)
    
    #record array
    record = np.ones_like(time)     
    
    for i in time:
        for j in np.arange(len(X)):
                        
            # store it 
            oldX = X[j]
            
            # store old f
            oldF = u(oldX)
            
            # move it
            increment = np.random.normal(0, strength, 1)
            	
            if increment < -1.0 * oldX:
                X[j]= increment + 1.0 + oldX
            elif increment > 1.0 - oldX:
                X[j]= increment - 1.0 + oldX
            else:
                X[j]= increment + oldX
            	
            # new energy
            newF = u(X[j])
             
            if newF > oldF:
                # apply metropolis
                if np.random.uniform(low=0.0,high=1.0,size=1) > np.exp(-beta*(newF-oldF)):
                    #discard it
                    X[j]=oldX
            
            # record 
            if j == 0:
                record[i] = X[j]   
                           
    
    # plot the distribution
    n, bins, patches = plt.hist(X, 30, normed=1, facecolor='blue', alpha=0.5)    
    plt.plot(bins,u(bins),'g',label="potential")
    plt.legend(loc="best")
    plt.show()
    
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02
    rect_histy = [left_h, bottom, 0.2, height]
    rect_scatter = [left, bottom, width, height]
    axHisty = plt.axes(rect_histy)
    axScatter = plt.axes(rect_scatter)
    
    # get the array for the potential
    xMin = np.min(record)
    xMax = np.max(record)
    plotX = np.linspace(0.0,1.0,num=100)
    
    # plot the potential
    nullfmt = NullFormatter()         # no labels
    axHisty.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    axHisty.plot(u(plotX),plotX,'g')
    
    axScatter.set_ylim([0.,1.])
    axScatter.plot(time,record)
    plt.show()
    
    # save it to file
    np.savetxt(filename,X)    
else:
    print errormessage