#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 11:04:01 2017

@author: konbob
"""

"""
Script for testing the 2 level 2 person game
"""

import numpy as np
import sys


# ========= constants
errormessage = "usage: python "+sys.argv[0]

# global variables 
argc = len(sys.argv)


# ========= functions

def process(y,q,xOther,fmp,fpm,wmp,wpm):
    # flip coin for mutation or communication
    if  np.random.uniform(0.,1.,1) < q:
        # do communication        
        if xOther == 1. and y == 0. and np.random.uniform(0.,1.,1) < fmp:
            return 1.                
        if xOther == 0. and y == 1. and np.random.uniform(0.,1.,1) < fpm:
            return 0.                
    else:
        # do mutation        
        if y == 0 and np.random.uniform(0.,1.,1) < wmp:
            return 1.              
        if y == 1 and np.random.uniform(0.,1.,1) < wpm:
            return 0.
    
    # if nothing happens
    return y
            
# ========= main =======================
# parse the arguments
if argc != 1:    
    print errormessage
        
# fix parmeters
q = 0.5
fmp = 0.487179
fpm = 0.539811
wmp = 1
wpm = 0.95
T = int(1e5)

# set up record
record = np.zeros([2,T])

# record loop
for t in np.arange(T):

    # start configuration
    x = np.array([0,1]) # 0 is minus/down and 1 is plus/up
    
    # time loop
    for time in np.arange(50):
        # process first one
        x[0] = process(x[0],q,x[1],fmp,fpm,wmp,wpm)
        
        # process the second one
        x[1] = process(x[1],q,x[0],fmp,fpm,wmp,wpm)
    
    # save for record
    record[:,t] = x
    
    
# matrix for counts
countMatrix = np.zeros([2,2])
    
# analysis of record
for t in np.arange(T): 
     if record[0,t] == 1 and record[1,t]==1:
         countMatrix[0,0]+=1
     if record[0,t] == 0 and record[1,t]==0:
         countMatrix[1,1]+=1
     if record[0,t] == 0 and record[1,t]==1:
         countMatrix[0,1]+=1
     if record[0,t] == 1 and record[1,t]==0:
         countMatrix[1,0]+=1
    
print countMatrix/T 
print "norm = " + str(np.sum(countMatrix/T))   
    
# error matrix
print  np.sqrt(countMatrix/T*(np.ones_like(countMatrix)-countMatrix/T))/np.sqrt(T)   
    
        