#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 14:46:52 2017

@author: konbob
"""
import numpy as np

# ========= functions
def getRawMoment(dataset,order):
    # sort data set
    np.ndarray.sort(dataset)
    
    # prepare the sum'
    res = 0.
    
    # compute powers
    for element in dataset:       
        res+= np.power(element,order)
        
    # normalize and return
    return res/(len(dataset))

# prepare test data
N = 1000
sigma = 1.0
data = np.random.normal(0.0,scale=sigma,size=N)

print getRawMoment(data,1.)
print getRawMoment(data,2.)
print getRawMoment(data,3.)
print getRawMoment(data,4.)