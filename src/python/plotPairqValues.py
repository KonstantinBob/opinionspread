#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 14:32:01 2017

@author: konbob
"""

"""
Script for testing the 2 level 2 person game
"""

import numpy as np
import sys
import matplotlib.pyplot as plt

# Setup Latex for Plots
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'font.family': 'lmodern'})
plt.rcParams.update({'font.size': 26})
plt.rcParams.update({'backend': 'PDF'})

plt.axes().set_aspect('equal', 'datalim')   

# ========= constants =
errormessage = "usage: python "+sys.argv[0] + "pValue datafile1  datafile2 datafile3"

# global variables 
argc = len(sys.argv)


# ========= functions

def PlusPlus(q,wUp,wDwn,fUp,fDwn):
    m=1-q
    return (wUp*(m*wUp+q*fUp))/(wUp*(m*wUp+q*fUp)+2.*m*wUp*wDwn+wDwn*(m*wDwn+q*fDwn))
    
def PlusMinus(q,wUp,wDwn,fUp,fDwn):
    m=1-q
    return (2.*m*wUp*wDwn)/(wUp*(m*wUp+q*fUp)+2.*m*wUp*wDwn+wDwn*(m*wDwn+q*fDwn))
    
def MinusMinus(q,wUp,wDwn,fUp,fDwn):
    m=1-q
    return (wDwn*(m*wDwn+q*fDwn))/(wUp*(m*wUp+q*fUp)+2.*m*wUp*wDwn+wDwn*(m*wDwn+q*fDwn))
            
# ========= main =======================
# parse the arguments
if argc != 5:    
    print errormessage
    exit
        
p = float(sys.argv[1])
    
filename1=sys.argv[2]
filename2=sys.argv[3]
filename3=sys.argv[4]

# read the files
data1 = np.loadtxt(filename1,usecols=0)
data2 = np.loadtxt(filename2,usecols=0)
data3 = np.loadtxt(filename3,usecols=0)


# parameters
wUp = 1
wDwn = p
fUp = p/(1+p)
fDwn = np.max([1.,1./(p*(1+p))])

# simulation q
qSim = 0.1*np.arange(0,11)


# theory
qValues = np.linspace(0.0,1.0,50)
theoPP = PlusPlus(qValues,wUp,wDwn,fUp,fDwn)
theoMix = PlusMinus(qValues,wUp,wDwn,fUp,fDwn)
theoMM = MinusMinus(qValues,wUp,wDwn,fUp,fDwn)

# ++
ax = plt.subplot(111)
ax.spines["top"].set_visible(False) 
ax.spines["right"].set_visible(False)  
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.plot(qValues,theoPP,label="Theory",alpha=0.5,lw=2)
plt.plot(qSim,data1[0::4],'^',label='ER',ms=12,fillstyle='none')
plt.plot(qSim,data2[0::4],'v',label='FF',ms=12,fillstyle='none')
plt.plot(qSim,data3[0::4],'o',label='SK',ms=12,fillstyle='none')
legend=plt.legend(loc='best',numpoints=1,frameon=False)
legend.get_title().set_fontsize('21') #legend 'Title' fontsize
plt.setp(plt.gca().get_legend().get_texts(), fontsize='12') #legend 'list' fontsize
plt.ylim(0.15,1.1)
plt.xlabel('$q$')
plt.ylabel('$P_{\uparrow\uparrow}$')
plt.tight_layout()
plt.savefig('pp.pdf')
plt.show()


# --
ax = plt.subplot(111)
ax.spines["top"].set_visible(False) 
ax.spines["right"].set_visible(False)  
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.plot(qValues,theoMM,label="Theory",alpha=0.5,lw=2)
plt.plot(qSim,data1[1::4],'^',label='ER',ms=12,fillstyle='none')
plt.plot(qSim,data2[1::4],'v',label='FF',ms=12,fillstyle='none')
plt.plot(qSim,data3[1::4],'o',label='SK',ms=12,fillstyle='none')
legend=plt.legend(loc='best',numpoints=1,frameon=False)
legend.get_title().set_fontsize('21') #legend 'Title' fontsize
plt.setp(plt.gca().get_legend().get_texts(), fontsize='12') #legend 'list' fontsize
plt.xlim(-.05,1.05)
#plt.ylim(-.025,0.3)
plt.xlabel('$q$')
plt.ylabel('$P_{\downarrow\downarrow}$')
plt.tight_layout()
plt.savefig('mm.pdf')
plt.show()


# +-
ax = plt.subplot(111)
ax.spines["top"].set_visible(False) 
ax.spines["right"].set_visible(False)  
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.plot(qValues,theoMix,label="Theory",alpha=0.5,lw=2)
plt.plot(qSim,data1[2::4],'^',label='ER',ms=12,fillstyle='none')
plt.plot(qSim,data2[2::4],'v',label='FF',ms=12,fillstyle='none')
plt.plot(qSim,data3[2::4],'o',label='SK',ms=12,fillstyle='none')
legend=plt.legend(loc='best',numpoints=1,frameon=False)
legend.get_title().set_fontsize('21') #legend 'Title' fontsize
plt.setp(plt.gca().get_legend().get_texts(), fontsize='12') #legend 'list' fontsize
plt.xlim(-.05,1.05)
#3plt.ylim(-.025,0.3)
plt.xlabel('$q$')
plt.ylabel('$P_{\uparrow\downarrow}$')
plt.tight_layout()
plt.savefig('pm.pdf')
plt.show()

# -+
ax = plt.subplot(111)
ax.spines["top"].set_visible(False) 
ax.spines["right"].set_visible(False)  
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.plot(qValues,theoMix,label="Theory",alpha=0.5,lw=2)
plt.plot(qSim,data1[3::4],'^',label='ER',ms=12,fillstyle='none')
plt.plot(qSim,data2[3::4],'v',label='FF',ms=12,fillstyle='none')
plt.plot(qSim,data3[3::4],'o',label='SK',ms=12,fillstyle='none')
legend=plt.legend(loc='best',numpoints=1,frameon=False)
legend.get_title().set_fontsize('21') #legend 'Title' fontsize
plt.setp(plt.gca().get_legend().get_texts(), fontsize='12') #legend 'list' fontsize
plt.xlim(-.05,1.05)
#plt.ylim(-.025,0.3)
plt.xlabel('$q$')
plt.ylabel('$P_{\downarrow\uparrow}$')
plt.tight_layout()
plt.savefig('mp.pdf')
plt.show()







