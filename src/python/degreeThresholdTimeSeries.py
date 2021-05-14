#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 15:49:16 2017

@author: konbob
"""

import numpy as np
import sys
import matplotlib.pyplot as plt

# Setup Latex for Plots
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'font.family': 'lmodern'})
plt.rcParams.update({'font.size': 26})
plt.rcParams.update({'backend': 'PDF'})

#plt.axes().set_aspect('equal', 'datalim')   

# ========= constants =
errormessage = "usage: python "+sys.argv[0] + " plotname datafile1"

# global variables 
argc = len(sys.argv)


# ========= functions

            
# ========= main =======================
# parse the arguments
if argc != 3:    
    print errormessage
    exit
    
plotname =sys.argv[1]
        
filename = sys.argv[2]

above,below = np.loadtxt(filename,unpack=True)


#plot results

ax = plt.subplot(111)
ax.spines["top"].set_visible(False) 
ax.spines["right"].set_visible(False)  
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

plt.plot(above,label="above")
plt.plot(below,label="below")
    
plt.ylim([0.0,1.0])
legend=plt.legend(numpoints=1,frameon=False,ncol=3,loc='upper center', bbox_to_anchor=(0.5, -0.25))
legend.get_title().set_fontsize('21') 
plt.setp(plt.gca().get_legend().get_texts(), fontsize='26')

plt.xlabel('Time (Cycles)')
plt.ylabel('$\\frac{N_{\\uparrow}}{N}$')
#text = plt.figtext(0.8,0.15,'$p=0.95$',fontsize=12)
#plt.savefig(str(plotname)+'.pdf',bbox_extra_artists=(legend,), bbox_inches='tight')
plt.savefig(str(plotname)+'.pdf', bbox_inches='tight')
plt.show()