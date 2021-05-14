#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 15:32:50 2017

@author: konbob
"""

import numpy as np
import sys
import matplotlib.pyplot as plt
import cPickle as pickle
import pymc3 as pm
import pandas as pd
from collections import OrderedDict
from scipy.stats import norm
import warnings
warnings.filterwarnings("ignore")

from multiprocessing import Pool


# Setup Latex for Plots
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'font.family': 'lmodern'})
plt.rcParams.update({'font.size': 26})
plt.rcParams.update({'backend': 'PDF'})
  
# ========= constants =
errormessage = "usage: python "+sys.argv[0] + "peaksfile debug flag (-d for true, -n for false) "

# global variables 
argc = len(sys.argv)

# ========= functions
def Model(x,a,b,c):
    
    return a*(x**(-b))+c

def fit(x,y,meanVec,stdVec,errors):
    
    aMu,bMu,cMu = meanVec
    aStd,bStd,cStd = stdVec
    
    model = pm.Model()
    
    if False:    
        df = pd.DataFrame(np.transpose([x,y,errors]),columns=['x','y','error'])
        print df
    
    with model:    
        # Priors for unknown model parameters
        a = pm.Normal('a', mu=aMu, sd=aStd)
        b = pm.Normal('b', mu=bMu, sd=bStd)
        c = pm.Normal('c', mu=cMu, sd=cStd)    
            
        # Expected value of outcome
        mu = Model(x,a,b,c)
    
        # Likelihood (sampling distribution) of observations
        Like = pm.Normal('Like', mu=mu, sd=errors, observed=y)
        
        # do sampling        
        trace = pm.sample(1000,progressbar=False,init='ADVI',step = pm.NUTS(),njobs=1)
        
        # give summary
        summary = pm.df_summary(trace)
        
        return summary
   
    
def findInfinitePeak(dataSet,debug):
       
    # prepare numpy arra
    data = np.array(dataSet,dtype=float)
       
    # sort it by size
    data = data[data[:,0].argsort()]
    
    # extract x and y points
    x = data[:,0]
    y = data[:,1]       
    err = data[:,2]
    
    # rescale for numerical stability
    xScale = x / np.max(x)
    
    # prepare initial guesses    
    aGuess = 0.02   
    bGuess = 0.5 
    cGuess = np.min(y)

    vecMean = [aGuess,bGuess,cGuess]
    vecStd = [.04,.25,0.2]            
    
    summary = fit(xScale,y,vecMean,vecStd,err)      
                 
    # === calculate the infinite value
    aBest = summary.loc['a']['mean']
    aError = summary.loc['a']['sd']
    bBest = summary.loc['b']['mean']
    bError = summary.loc['b']['sd']
    cBest = summary.loc['c']['mean']
    cError = summary.loc['c']['sd']
        
    if debug:
        print pd.DataFrame(data=np.column_stack((xScale,y,err)),columns=['xScale','y','yerr'])                
        guess =  pd.DataFrame(data=[aGuess,bGuess,cGuess],columns=['Guess'],index=['a','b','c'],)        
        print pd.concat([summary, guess], axis=1)
    
    # plot of results        
    plt.errorbar(xScale,y,yerr=err,fmt='o',label="Data")
    xPlot = np.linspace(0.,1.)
    plt.plot(xPlot,Model(xPlot,aBest,bBest,cBest),label="Fit")
    ax = plt.axhline(y = cBest,linestyle = '--',label="$\\phi(\\infty)$")
    colorAx = ax.get_color()
    #plt.axhline(y = cBest+cError,linestyle = '-',alpha=0.3,color=colorAx)
    #plt.axhline(y = cBest-cError,linestyle = '-',alpha=0.3,color=colorAx)
    plt.axhspan(cBest-cError,cBest+cError,linestyle = '-',alpha=0.3,color=colorAx)
    
    # make nice legend
    handles, labels = plt.gca().get_legend_handles_labels()    
    by_label = OrderedDict(zip(labels, handles))    
    legend=plt.legend(by_label.values(), by_label.keys(),numpoints=1,frameon=False,ncol=3,loc='upper center', bbox_to_anchor=(0.5, -0.25))
    legend.get_title().set_fontsize('21') 
    plt.setp(plt.gca().get_legend().get_texts(), fontsize='26')
    if debug:                   
        plt.show()
    else:
        
        plt.savefig('Fit_'+str(size)+'_'+str(p)+'_'+str(q)+'.pdf',bbox_inches = 'tight')
        plt.title('$p='+str(p)+'$ $q='+str(q)+'$')
        plt.savefig('Fit_'+str(size)+'_'+str(p)+'_'+str(q)+'.png',bbox_inches = 'tight')
        plt.close()
                      
                        
    # return infintite value               
    return cBest,cError,bBest,bError,aBest,aError


def showDataSet(dataSet,p,q,size,debug):
    # prepare numpy arra
    data = np.array(dataSet,dtype=float)
       
    # sort it by size
    data = data[data[:,0].argsort()]
    
    # extract x and y points
    x = data[:,0]
    y = data[:,1]
    err = data[:,2]
    
    # rescale for numerical stability
    xScale = x / np.max(x)
    
    ax = plt.subplot(111)
    ax.spines["top"].set_visible(False) 
    ax.spines["right"].set_visible(False)  
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    
    #plot 
    plt.errorbar(xScale,y, yerr=err, fmt='o')
    plt.xlabel('$N_{Scale}$')    
    plt.ylabel('$\\phi$')
    
    if debug:
        plt.show()
    
                   
# ========= main =======================
# *
# *
# ======================================

#  ========= parse the arguments ========= 
if argc != 3:    
    print errormessage
    sys.exit()
 
#get file with results
resultFile =  str(sys.argv[1])+".pickle"

#debugging
debug = False
if str(sys.argv[2]) == '-d':
    debug = True
elif str(sys.argv[2]) != '-n':
    print errormessage
    sys.exit()

# size set
sizeSet = set()

# qSet 
qSet = set()

#pSet
pSet = set()

# dict for exponents
expoDict = dict()
expoErrorDict = dict()

# dict for scale
scaleDict = dict()
scaleErrorDict = dict()

#  =========  load results from file ========= 
filePeaks = open(resultFile, "rb")
results = pickle.load(filePeaks)
filePeaks.close()

# store values
for key in results:
    p,q,size = key
    sizeSet.add(size)
    pSet.add(p)
    qSet.add(q)

#  ========= finite size scaling    ========= 

# prepare the peak list
infinitePeaksHigh = []
infinitePeaksLow = []
infinitePeaksMiddle = []

for q in qSet:
    for p in pSet:
        # prepare the peak lists
        lowerPeaks = []
        singlePeaks = []
        higherPeaks = []   
            
        fileOutput = open(str(p)+'_'+str(q)+'.dat','w')   
        
        # get data points 
        for size in sizeSet:
            key  = p,q,size
            if key in results:
                peaks = results[key][0]
                errors = results[key][1]     
                                   
                if len(peaks) == 1:                
                    singlePeaks.append((size,peaks[0],errors[0]))
                    fileOutput.write(str(size)+' '+str(peaks[0])+' '+str(errors[0])+'\n')
                else:
                    higherPeaks.append((size,np.max(peaks),errors[np.argmax(peaks)]))                
                    lowerPeaks.append((size,np.min(peaks),errors[np.argmin(peaks)]))
                    fileOutput.write(str(size)+' '+str(np.max(peaks))+' '+str(errors[np.argmax(peaks)])+' '+str(np.min(peaks))+' '+str(errors[np.argmin(peaks)])+'\n')
                                      
        doFit = 'y'
        # do the fits for infinite size
        if len(lowerPeaks) > 1:
           
           showDataSet(lowerPeaks,p,q,size,debug)
           
           if debug:
               doFit = raw_input('Do fit? [y/n]')
           if doFit == 'y':
               peak,error,expo,expoError,aBest,aError = findInfinitePeak(lowerPeaks,debug)
               infinitePeaksLow.append([q,peak,error])
               expoDict[q] = expo
               expoErrorDict[q] = expoError
               scaleDict[q] = aBest
               scaleErrorDict[q] = aError
        
        if len(higherPeaks) > 1:
            
           showDataSet(higherPeaks,p,q,size,debug)
           if debug:
               doFit = raw_input('Do fit? [y/n]')
           if doFit == 'y':
               peak,error,expo,expoError,aBest,aError = findInfinitePeak(higherPeaks,debug)
               infinitePeaksHigh.append([q,peak,error])
               expoDict[q] = expo
               expoErrorDict[q] = expoError
               scaleDict[q] = aBest
               scaleErrorDict[q] = aError
               
        if len(singlePeaks) > 1:
            
           showDataSet(singlePeaks,p,q,size,debug)
           if debug:
               doFit = raw_input('Do fit? [y/n]')
           if doFit == 'y':
               peak,error,expo,expoError,aBest,aError = findInfinitePeak(singlePeaks,debug)
               infinitePeaksMiddle.append([q,peak,error])   
               expoDict[q] = expo
               expoErrorDict[q] = expoError
               scaleDict[q] = aBest
               scaleErrorDict[q] = aError
               
        fileOutput.close()
           
# ========= Save peaks
fileInfPeaksHigh = open(str(sys.argv[1])+"InfHigh.pickle",'wb')
pickle.dump(infinitePeaksHigh,fileInfPeaksHigh)
fileInfPeaksHigh.close()
      
fileInfPeaksLow = open(str(sys.argv[1])+"InfLow.pickle",'wb')
pickle.dump(infinitePeaksLow,fileInfPeaksLow)
fileInfPeaksLow.close()

fileInfPeaksMiddle = open(str(sys.argv[1])+"InfMiddle.pickle",'wb')
pickle.dump(infinitePeaksMiddle,fileInfPeaksMiddle)
fileInfPeaksMiddle.close()

# ===== ExpoPlot ====
qArray = np.fromiter(iter(expoDict.keys()), dtype=float)
expoArray = np.fromiter(iter(expoDict.values()), dtype=float)
errorArray = np.fromiter(iter(expoErrorDict.values()), dtype=float)

ax = plt.subplot(111)
ax.spines["top"].set_visible(False) 
ax.spines["right"].set_visible(False)  
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

plt.errorbar(qArray,expoArray,yerr=errorArray,fmt='o')
plt.xlabel("$q$")
plt.ylabel("$b$")
plt.savefig('Expo_'+ str(sys.argv[1])+ '.pdf', bbox_inches='tight')
plt.close()

# ===== ScalePlot ====
qArray = np.fromiter(iter(scaleDict.keys()), dtype=float)
scaleArray = np.fromiter(iter(scaleDict.values()), dtype=float)
errorArray = np.fromiter(iter(scaleErrorDict.values()), dtype=float)

ax = plt.subplot(111)
ax.spines["top"].set_visible(False) 
ax.spines["right"].set_visible(False)  
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

plt.errorbar(qArray,scaleArray,yerr=errorArray,fmt='o')
plt.xlabel("$q$")
plt.ylabel("$a$")
plt.savefig('Scale_'+ str(sys.argv[1])+ '.pdf', bbox_inches='tight')
plt.close()

















