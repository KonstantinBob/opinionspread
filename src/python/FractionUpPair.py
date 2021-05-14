#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 10:37:14 2017

@author: konbob
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
from collections import OrderedDict

# Setup Latex for Plots
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'font.family': 'lmodern'})
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'backend': 'PDF'})
params = {'text.latex.preamble' : [r'\usepackage{amsmath}',r'\usepackage{amssymb}']}
plt.rcParams.update(params)

# ========= constants
errormessage = "usage: python "+sys.argv[0]+" plotname "

# global variables 
argc = len(sys.argv)

# ========= functions

def upup(m,q,wUp,wDwn,fUp,fDwn):
    return (wUp*(m*wUp+q*fUp))/(wUp*(m*wUp+q*fUp)+2.*m*wUp*wDwn+wDwn*(m*wDwn+q*fDwn))

def mix(m,q,wUp,wDwn,fUp,fDwn):
    return (1.*m*wUp*wDwn)/(wUp*(m*wUp+q*fUp)+2.*m*wUp*wDwn+wDwn*(m*wDwn+q*fDwn))

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
if argc != 2:
    print errormessage
else:
    plotname =  sys.argv[1]   

# theory
q = np.linspace(0.,1.)

p1= 0.95
p2 = 0.5
p3 = 0.7 
 
# plot for fixed q
ax = plt.subplot(111)
ax.spines["top"].set_visible(False) 
ax.spines["right"].set_visible(False)  
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')


pl1 = plt.plot(q,frac(p1,q),ls='--',label="p=0.95",alpha=0.8)
color1 = pl1[0].get_color()

pl2 = plt.plot(q,frac(p2,q),ls='--',label="p=0.5",alpha=0.8)
color2 = pl2[0].get_color()

pl3 = plt.plot(q,frac(p3,q),ls='--',label="p=0.7",alpha=0.8)
color3 = pl3[0].get_color()

#plot data

### ==== begin data ER ===================
plt.plot(0.5,0.525831663327,'^',color=color1) # p=0.95
plt.plot(0.7,0.54248496994,'^',color=color1)   # p=0.95
plt.plot(0.95,0.688577154309,'^',color=color1)  # p=0.95

plt.plot(0.95,0.933082164329,'^',color=color3) # p = 0.7

plt.plot(0.95,0.965390781563,'^',color=color2) # p = 0.5

### ==== end data er ===================


plt.xlim(0.,1.)

handles, labels = plt.gca().get_legend_handles_labels()    
by_label = OrderedDict(sorted(zip(labels, handles), key = lambda x : str(x[0]) ))

legend=plt.legend(by_label.values(), by_label.keys(),numpoints=1,frameon=False,ncol=3,loc='upper center', bbox_to_anchor=(0.5, -0.25))
legend.get_title().set_fontsize('16') 
plt.setp(plt.gca().get_legend().get_texts(), fontsize='16')

plt.xlabel("$q$")
plt.ylabel("$\\mathbb{E}\\left [\\frac{N_{\\uparrow}}{N} \\right ]$")
plt.savefig(str(plotname)+'.pdf', bbox_inches='tight')
plt.show()