#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 15:05:44 2017

@author: konbob
"""

#import os
#os.environ["THEANO_FLAGS"] = "blas.ldflags=\"-L/usr/lib/ -lblas\""

import numpy as np
import sys
import matplotlib.pyplot as plt
import pymc3 as pm
from scipy.stats import uniform
import warnings
warnings.filterwarnings("ignore")
import math 

from multiprocessing import Pool

# Setup Latex for Plots
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'font.family': 'lmodern'})
plt.rcParams.update({'font.size': 26})
plt.rcParams.update({'backend': 'PDF'})
  
# ========= constants =
errormessage = "usage: python "+sys.argv[0] + "plotname datafiles"

# global variables 
argc = len(sys.argv)

# ========= functions
def Model(x,a,b,c):
    
    return a*(x**(-b))+c

def ComputeLowerBondUniform(mu,std):
    return  mu - std*0.5*math.sqrt(12.)
    
def ComputeUpperBondUniform(mu,std):    
    return  mu + std*0.5*math.sqrt(12.)        

def fit(x,y,lowerVec,upperVec):
    
    lA,lB,lC = lowerVec
    uA,uB,uC = upperVec
    
    model = pm.Model()
    
    with model:    
        # Priors for unknown model parameters
                 
        a = pm.Uniform('a',lower=lA,upper=uA)
        b = pm.Uniform('b',lower=lB,upper=uB) 
        c = pm.Uniform('c',lower=lC,upper=uC) 
        
        # Expected value of outcome
        mu = Model(x,a,b,c)
    
        # Likelihood (sampling distribution) of observations
        Like = pm.Normal('Like', mu=mu, sd=0.1*np.ones_like(y), observed=y)
        
        # do sampling        
        trace = pm.sample(1000,progressbar=False,init='ADVI',step = pm.NUTS(),njobs=1)
        
        # give summary
        summary = pm.df_summary(trace)
        
        return summary
    
def fitFlat(x,y):      
    
    model = pm.Model()
    
    with model:    
        # Priors for unknown model parameters
        a = pm.Flat('a')
        b = pm.Flat('b')
        c = pm.Flat('c')    
            
        # Expected value of outcome
        mu = Model(x,a,b,c)
    
        # Likelihood (sampling distribution) of observations
        Like = pm.Normal('Like', mu=mu, sd=0.01*np.ones_like(y), observed=y)
        
        # do sampling        
        trace = pm.sample(1000,progressbar=False,init='ADVI',step = pm.NUTS(),njobs=1)
        
        # give summary
        summary = pm.df_summary(trace)
        
        return summary
    
def calcD(summary,thetaZero):      
    # start with the varianc term
    vec = np.array([summary.loc['a']['sd'],summary.loc['b']['sd'],summary.loc['c']['sd']],dtype=float)
    vec = np.power(vec,2.0)
   
    # add means term
    meansBroad = np.array([summary.loc['a']['mean'],summary.loc['b']['mean'],summary.loc['c']['mean']],dtype=float)
    vec += np.power(meansBroad-thetaZero,2.0)    
        
    # return it
    return vec

def job(randomIndices):
        
    # reduced sample
    yRed = y[randomIndices]
    xRed = x[randomIndices]                 
          
    # broad run                          
    summaryBroad = fitFlat(xRed,yRed)
    dBroad = calcD(summaryBroad,thetaZero)               
    
    # fine run                  
    summaryFine = fit(xRed,yRed,vecMu,vecStd)
    dFine = calcD(summaryFine,thetaZero)
        
    return [dBroad,dFine]

def jobInitialiser(_vecMu,_vecStd,_thetaZero,_x,_y):
    global vecMu,vecStd,constBroad,thetaZero,x,y
    vecMu = _vecMu
    vecStd = _vecStd    
    thetaZero = _thetaZero
    x = _x
    y = _y

def computeEffSampleSizes(x,y,lowerVec,upperVec):
        
    if False:
        # show the pdfs
        for index in np.arange(len(lowerVec)):
            low = lowerVec[index]
            up = upperVec[index]
            PlotX = np.linspace(-1.5*low,1.5*up)
            pdfFine = uniform.pdf(PlotX,loc=low,scale =up-low )
            pdfBroad = np.ones_like(PlotX)/(3.0*(up-low))
            #plt.ylim([0.,0.5])
            
            plt.plot(PlotX,pdfFine,'-')
            plt.plot(PlotX,pdfBroad,'--')
            plt.title('index='+str(index))
            
            plt.show()
    
    # === get theta zero
    # run with broad prior
    summaryBroad = fitFlat(x,y)
    thetaZero = np.array([summaryBroad.loc['a']['mean'],summaryBroad.loc['b']['mean'],summaryBroad.loc['c']['mean']],dtype=float)    
    
    # === prepare runs with different sample sizes    
    # calculate max subsample size
    kMax = max(2*int(np.ceil(np.sqrt(len(y)))),6)
    
    # number of bootstraps
    bootstrapMax = 30
    
    # Storage for posterior uncertrainties
    Ubroad = np.zeros((kMax,3,bootstrapMax))
    Ufine = np.zeros((kMax,3,bootstrapMax))
    
    # === Run with different sample sizes
    for sampleSize in np.arange(1,kMax+1):   
          
        # prepare tuple of input
        argsJob = ()
                
        for indexPara in np.arange(bootstrapMax):
            # ====== get bootstraped data points   
            
            # get random shuffeld indices
            randomIndices = np.random.permutation(np.arange(len(x)))[:sampleSize]
                         
            # prepare args for parallelization
            argsJob = argsJob +(randomIndices,) 

        # === parallel run
        try:
            maxWorkers = 8
            pool = Pool(maxWorkers,jobInitialiser,(lowerVec,upperVec,thetaZero,x,y))
            
            
            # the heavy lifting happens here
            resObj =  pool.map_async(job,argsJob)
            
            # free the resources
            pool.close()
            pool.join()
            
            res = resObj.get()
            
            # store the results
            for indexResult in np.arange(len(res)):
                Ubroad[sampleSize-1,:,indexResult] = res[indexResult][0]
                Ufine[sampleSize-1,:,indexResult] = res[indexResult][1]
        except KeyboardInterrupt:
            print "Caught KeyboardInterrupt, terminating workers"
            pool.terminate()
            pool.join()
            sys.exit()
                    
    # average the results
    UbroadRed = np.mean(Ubroad,axis=2)
    UfineRed = np.mean(Ufine,axis=2) 
    
    UbroadError = np.std(Ubroad,axis=2)/np.sqrt(bootstrapMax)
    UfineError = np.std(Ufine,axis=2)/np.sqrt(bootstrapMax)
    
    titleDict = {}
    titleDict[0] = '$a$'
    titleDict[1] = '$b$'
    titleDict[2] = '$c$'
    
    # plot results =====================       
       
    for index in np.arange(3):
        ax = plt.subplot(111)
        ax.spines["top"].set_visible(False) 
        ax.spines["right"].set_visible(False)  
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        
        sampleSizeArray = np.arange(1,kMax+1)
        plt.errorbar(sampleSizeArray,UfineRed[:,index],yerr=UfineError[:,index],fmt='o',label='Informative')
        plt.errorbar(sampleSizeArray,UbroadRed[:,index],yerr=UbroadError[:,index],fmt='o',label='Flat')
        
        plt.legend(numpoints=1,frameon=False,ncol=2,loc='upper center', bbox_to_anchor=(0.5, -0.25))
        plt.xlabel('Subsample size $k$')
        plt.ylabel("$U(k)$")
        plt.savefig('U'+str(index)+'.pdf', bbox_inches='tight')
        plt.close()
        
        # write data to file
        np.savetxt('U_fine_red-'+str(index)+'.txt',UfineRed[:,index])
        np.savetxt('U_fine-error-'+str(index)+'.txt',UfineError[:,index])
        np.savetxt('U_broad_red-'+str(index)+'.txt',UbroadRed[:,index])
        np.savetxt('U_broad_error-'+str(index)+'.txt',UbroadError[:,index])
        
# ========= main =======================
# *
# *
# ======================================

#  ========= parse the arguments ========= 
if argc > 1:    
    print(errormessage)
    sys.exit()  

# ==== prepare test data ====

aTrue = 0.02
bTrue = 0.5
cTrue = 0.25

x = np.linspace(0.01,1,num=8)
xFine = np.linspace(0.01,1,num=100)

errors = np.random.normal(loc=0.0,scale= 0.01, size = len(x))

y = Model(x,aTrue,bTrue,cTrue) + errors


# ==== prepare the  model ====

lowerVec = np.array([0,0,0])
upperVec = np.array([.1,1.0,1.0])

# === fit the model

summary = fit(x,y,lowerVec,upperVec)

print summary

aBest = summary.loc['a']['mean']
bBest = summary.loc['b']['mean']
cBest = summary.loc['c']['mean']

if True:  
    plt.plot(x,y,'x')
    plt.plot(xFine,Model(xFine,aBest,bBest,cBest))    
    plt.show()    
    

# get effective prior sample sizes
computeEffSampleSizes(x,y,lowerVec,upperVec)  