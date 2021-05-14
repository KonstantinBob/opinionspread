#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 14:20:15 2017

@author: konbob
"""

import numpy as np
import matplotlib.pyplot as plt
import sys


# Setup Latex for Plots
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'font.family': 'lmodern'})
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'backend': 'PDF'})

# ========= constants
errormessage = "usage: python "+sys.argv[0]+" filenames "

# global variables 
argc = len(sys.argv)

# ========= functions

def upup(m,q,wUp,wDwn,fUp,fDwn):
    return (wUp*(m*wUp+q*fUp))/(wUp*(m*wUp+q*fUp)+2.*m*wUp*wDwn+wDwn*(m*wDwn+q*fDwn))

def mix(m,q,wUp,wDwn,fUp,fDwn):
    return (2.*m*wUp*wDwn)/(wUp*(m*wUp+q*fUp)+2.*m*wUp*wDwn+wDwn*(m*wDwn+q*fDwn))

def downdown(m,q,wUp,wDwn,fUp,fDwn):
    return (wDwn*(m*wDwn+q*fDwn))/(wUp*(m*wUp+q*fUp)+2.*m*wUp*wDwn+wDwn*(m*wDwn+q*fDwn))

def f(p):
    return np.minimum(1./(p*(1.+p)),1.)
      
def g(p):
    return p/(1.+p)

def frac(p,q):
    return upup(1-q,q,1.,p,f(p),1-f(p))+mix(1-q,q,1.,p,f(p),1-f(p))

# ========= main =======================
# parse the arguments
if argc < 2:
    print errormessage
    exit
     
# read the file
for arg in sys.argv[1:]:
    filename = str(arg)
    print filename, np.mean(np.loadtxt(filename))

