#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 10:04:20 2017

@author: konbob

Plots the phase diagram and computes the critical point.
Requires runs of finiteSizeScaling and findPeaksPhaseDiagram before!

"""
import numpy as np
import sys
import matplotlib.pyplot as plt
from collections import OrderedDict
from scipy.optimize import brentq 
import cPickle as pickle
import pymc3 as pm
import pandas as pd
import os.path

# Setup Latex for Plots
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'font.family': 'lmodern'})
plt.rcParams.update({'font.size': 26})
plt.rcParams.update({'backend': 'PDF'})
params = {'text.latex.preamble' : [r'\usepackage{amsmath}',r'\usepackage{amssymb}']}
plt.rcParams.update(params)
  
# ========= constants =
errormessage = "usage: python "+sys.argv[0] + " datafiles debug flag (-d for true, -n for false) plotname doFit -d for true, -n for false) isIsing (y/n)"

# global variables 
argc = len(sys.argv)

# ========= functions
def powerLaw(x,a,b,c):
    return c+a*(x**b)
    
   
def fitBoundary(dataSet,debug):
           
    # prepare numpy arra
    data = np.array(dataSet,dtype=float)
    
    # sort it by size
    data = data[data[:,0].argsort()]
    
    # extract x and y points
    x = data[:,0]
    y = data[:,1]        
    error =  data[:,2]    
       
    # reduce y to m
    m = y - 0.5*np.ones_like(y)
    
    if debug:
        df = pd.DataFrame(data=np.column_stack((x,m,error)),columns=['x','m','error'])
        print df
        
    # ==== prepare initial guesses  
         
    aGuess = 0.02   
    bGuess = 10.0
    cGuess = np.min(m)
    
    if m[-1] < m[0]:      
       aGuess = -1.0 * aGuess
      
        
    # call fit function
    basic_model = pm.Model()
    with basic_model:
    
        # Priors for unknown model parameters
        a = pm.Normal('a', mu=aGuess, sd=0.1)
        b = pm.Normal('b', mu=bGuess, sd=3.0)         
        c = pm.Normal('c', mu=cGuess, sd=0.1)   
        
        # Expected value of outcome
        mu = powerLaw(x,a,b,c)
    
        # Likelihood (sampling distribution) of observations
        Like = pm.Normal('like', mu=mu, sd=error, observed=m)
       
   
    # === sample and evaluate the model    
    with basic_model:             
                  
        trace = pm.sample(1000,njobs=2,progressbar=False,init='ADVI',step = pm.NUTS())    
                                          
        # give summary
        summary = pm.df_summary(trace)      
        
        if debug:
            dfGuess = pd.DataFrame(data=[aGuess,bGuess,cGuess],index=['a','b','c'],columns=['guess'])            
            print pd.concat([summary, dfGuess], axis=1)
            
            
                 
    # === calculate the infinite value
    aBest = summary.loc['a']['mean']
    bBest = summary.loc['b']['mean']
    cBest = summary.loc['c']['mean']
    
    aError = summary.loc['a']['sd']
    bError = summary.loc['b']['sd']
    cError = summary.loc['c']['sd']
    
    if debug:
        plt.errorbar(x,m,yerr=error,fmt='o')
        xPlot = np.linspace(np.min(x),np.max(x))
        plt.plot(xPlot,powerLaw(xPlot,aBest,bBest,cBest))
        plt.show()
        
    return [aBest,bBest,cBest], [aError,bError,cError]


# ========= main =========================================================================================================================
# *
# *
# ========================================================================================================================================

#  ========= parse the arguments ========= 
if argc != 6:    
    print errormessage
    sys.exit()

 
#get file with results
resultFile =  str(sys.argv[1])
    
#debugging
debug = False
if str(sys.argv[2]) == '-d':
    debug = True
elif str(sys.argv[2]) != '-n':
    print errormessage
    sys.exit()
    
# fit flag
doFit = False
if str(sys.argv[4]) == '-d':
    doFit = True
elif str(sys.argv[4]) != '-n':
    print errormessage
    sys.exit()  
    
# fit flag
isIsing = False
if str(sys.argv[5]) == 'y':
    isIsing = True
elif str(sys.argv[5]) != 'n':
    print errormessage
    sys.exit()  
    
# store plotname    
plotname = str(sys.argv[3])

#  =========  load results from file ========= 
infinitePeaksHigh = []
nameHigh = resultFile+"InfHigh.pickle"
if os.path.isfile(nameHigh):    
    fileInfPeaksHigh = open(nameHigh,"rb")
    infinitePeaksHigh = pickle.load(fileInfPeaksHigh)
    fileInfPeaksHigh.close()
    
infinitePeaksLow = []
nameLow = resultFile+"InfLow.pickle"
if os.path.isfile(nameLow):      
    fileInfPeaksLow = open(nameLow,"rb")
    infinitePeaksLow = pickle.load(fileInfPeaksLow)
    fileInfPeaksLow.close()
    
infinitePeaksMiddle = []
nameMiddle = resultFile+"InfMiddle.pickle"
if os.path.isfile(nameMiddle):
    fileInfPeaksMiddle = open(nameMiddle,"rb")
    infinitePeaksMiddle = pickle.load(fileInfPeaksMiddle)
    fileInfPeaksMiddle.close()
    

# ========= fit continous boundary with infinite size
fitBounds = False 

if len(infinitePeaksHigh) > 1 and len(infinitePeaksLow) > 1 and doFit == True:
    fitBounds = True

if fitBounds:
    data = np.concatenate((infinitePeaksMiddle,infinitePeaksHigh),axis = 0)
    
    paramHigh,errorHigh = fitBoundary(data,debug)
    def boundHigh(x):
        return powerLaw(x,paramHigh[0],paramHigh[1],paramHigh[2])

if fitBounds:
    data = np.concatenate((infinitePeaksMiddle,infinitePeaksLow),axis = 0)
    
    paramLow, errorLow = fitBoundary(data,debug)
    def boundLow(x):
        return powerLaw(x,paramLow[0],paramLow[1],paramLow[2])
        
# ========= final plot phase diagramm =====================
ax = plt.subplot(111)
ax.spines["top"].set_visible(False) 
ax.spines["right"].set_visible(False)  
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

# for colors
cmap = plt.cm.get_cmap('plasma_r')

# load the peaks
fileAllPeaks = open(resultFile+'.pickle','rb')
allpeaks = pickle.load(fileAllPeaks)
fileAllPeaks.close()

# find max size in keys
sizeMax = float(max(allpeaks,key= lambda x: x[2])[2])

# find qmin and qmax
qMin = float(min(allpeaks,key= lambda x: x[1])[1])
qMax = float(max(allpeaks,key= lambda x: x[1])[1])

# plot them
for key in sorted(allpeaks,key= lambda x: x[2]):
    p,q,size = key
    peaks,errors = allpeaks[key]
    
    color = cmap(float(size)/sizeMax)
    
    for index in np.arange(len(peaks)):            
            plt.errorbar(q,peaks[index]-0.5,yerr=errors[index],fmt='o',color=color,label=str(int(size)))

infinitePeaks = infinitePeaksMiddle
if len(infinitePeaksLow) > 0:
    infinitePeaks = np.concatenate((infinitePeaks,infinitePeaksLow),axis=0)
    
if len(infinitePeaksHigh) > 0:
    infinitePeaks = np.concatenate((infinitePeaks,infinitePeaksHigh),axis=0)
    
for infPeak in infinitePeaks:
    q,m,err = infPeak
    plt.errorbar(q,m-0.5,yerr=err,fmt='o',color='k',label="$\\infty$")
    
if fitBounds:    
    qLinspace = np.linspace(qMin,qMax)                      
    lowPlot = plt.plot(qLinspace,boundLow(qLinspace),label="$\\mathrm{Fit}$")
    plt.plot(qLinspace,boundHigh(qLinspace),'--',label="$\\mathrm{Fit}$",color =  lowPlot[0].get_color())      

plt.xlabel('$q$')
if isIsing:
    plt.xlabel('$\\beta$')
plt.ylabel('$m$')
if isIsing:
    plt.ylabel('$m_{\\text{Ising}}$')

# === calculate critical point by intersection
if fitBounds:
    try:
        qCrit = brentq(lambda x : boundHigh(x)-boundLow(x),0.,1.0)
        
        qCritY = boundHigh(qCrit)
        
        plt.plot(qCrit,qCritY,'x',markersize=15,zorder=5)
        
        fileCrit = open(plotname+'crit.dat','w')
        
        fileCrit.write('qCrit = '+str(qCrit)+'\n')
        fileCrit.write('betaUp = '+ str(paramHigh[1])+' +- '+str(errorHigh[1])+'\n')
        fileCrit.write('betaDwn = ' + str(paramLow[1])+' +- '+str(errorLow[1])+'\n')
        fileCrit.close()
        
    except ValueError as e:
        print e

# remove duplicates and sort labels 
handles, labels = plt.gca().get_legend_handles_labels()    
by_label = OrderedDict(sorted(zip(labels, handles), key = lambda x : float(x[0]) if x[0].isdigit() else 1e5 + len(x[0]) ))

columns = 3
if fitBounds:
    columns =4
    
legend=plt.legend(by_label.values(), by_label.keys(),numpoints=1,frameon=False,ncol=columns,loc='upper center', bbox_to_anchor=(0.5, -0.25))
legend.get_title().set_fontsize('21') 
plt.setp(plt.gca().get_legend().get_texts(), fontsize='26')

plt.savefig(str(plotname)+'.pdf', bbox_inches='tight')   
#plt.show()
plt.close()