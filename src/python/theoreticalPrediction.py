#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 15:28:48 2017

@author: konbob

For plotting the theory
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
errormessage = "usage: python "+sys.argv[0]+" plotname "

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
    return 1./(p*(1.+p))

def g(p):
    return p/(1.+p)

# ========= main =======================
# parse the arguments
if argc != 2:
    print errormessage
else:
    plotname =  sys.argv[1]   
#prepare data
q = np.linspace(0.,1.)

p=0.95

theoryUp = upup(1-q,q,1.,p,f(p),1-f(p))
theoryMix = mix(1-q,q,1.,p,f(p),1-f(p))
theoryDwn = downdown(1-q,q,1.,p,f(p),1-f(p))

# plot for fixed p
ax = plt.subplot(111)
ax.spines["top"].set_visible(False) 
ax.spines["right"].set_visible(False)  
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')


plt.plot(q,theoryUp,label="$\\uparrow \\uparrow$")
plt.plot(q,theoryMix,label="$\\uparrow \\downarrow$ or $ \\downarrow\\uparrow$")
plt.plot(q,theoryDwn,label="$\\downarrow \\downarrow$")


plt.xlim(0.,1.)
#plt.ylim(0.,1.)

legend=plt.legend(numpoints=1,frameon=False,ncol=3,loc='upper center', bbox_to_anchor=(0.5, -0.25))
legend.get_title().set_fontsize('16') 
plt.setp(plt.gca().get_legend().get_texts(), fontsize='16')

plt.xlabel("$q$")
plt.ylabel('$P_{xy}$')
plt.savefig(str(plotname)+'.pdf', bbox_inches='tight')
plt.show()