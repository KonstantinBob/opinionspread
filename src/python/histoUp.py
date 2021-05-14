#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May  3 14:27:48 2017

@author: konbob

Plotting a histogram
"""


import numpy as np
import sys
import matplotlib.pyplot as plt

# Setup Latex for Plots
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'font.family': 'lmodern'})
plt.rcParams.update({'font.size': 26})
plt.rcParams.update({'backend': 'PDF'})

# ========= constants =
errormessage = "usage: python "+sys.argv[0] + " datafile1  datafile2 datafile3"

# global variables 
argc = len(sys.argv)


# ========= functions


            
# ========= main =======================
# parse the arguments
if argc != 2:    
    print errormessage
    exit
        
filename1=sys.argv[1]


# read the files
data1 = np.loadtxt(filename1)

# plot
ax = plt.subplot(111)
ax.spines["top"].set_visible(False) 
ax.spines["right"].set_visible(False)  
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.hist(data1,bins=np.arange(250),label="Theory",alpha=0.5,normed=True)
legend=plt.legend(loc='best',numpoints=1,frameon=False)
legend.get_title().set_fontsize('21') #legend 'Title' fontsize
plt.setp(plt.gca().get_legend().get_texts(), fontsize='12') #legend 'list' fontsize
#plt.ylim(0.15,1.1)
plt.xlabel('number Up')
plt.ylabel('PMF')
plt.tight_layout()
plt.savefig('histUp.pdf')
plt.show()


