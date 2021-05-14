#!/usr/bin/env python

# imports
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stat
import sys

# Setup Latex for Plots
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'font.family': 'lmodern'})
plt.rcParams.update({'font.size': 23})
plt.rcParams.update({'backend': 'PDF'})

# ========= constants
errormessage = "usage: "+sys.argv[0]+" inputfile"

# global variables 
argc = len(sys.argv)

# ========= functions
#empirical distribution function
def edf(dataSet,value):
    
    return np.searchsorted(dataSet, value, side='right') / float(dataSet.size)
    
def cdfBoltzmann(beta,value):
    return stat.expon.cdf(value,scale=1./beta)
    
def pdfBoltzmann(beta,value):
    return stat.expon.pdf(value,scale=1./beta)
    
def cdfUni(value):
    return stat.uniform.cdf(value)
        
def pdfUni(value):
    return stat.uniform.pdf(value) 
    
# ========= main =======================
if argc == 3:
    
    # get the names
    filename = sys.argv[1]  
    plotname = sys.argv[2] 
    
    # load data from inputfile
    data  = np.loadtxt(filename)
    
    
    # variables for later use
    N = len(data)
    sumData = np.sum(data)
    
    # calculate the mean of param beta
    betaEst = N/sumData 
    
    # calculate the error of param estimator
    betaError = np.sqrt(N)/sumData
    
    print "estimated beta ", betaEst, " +- ", betaError 
    
    # ========== Kolmogorov smirnov test ================
    
    # get end points of intervall
    x_min = np.min(data)
    x_max = np.max(data)
        
    # get the x array    
    x = np.linspace(x_min,x_max,500)
    
    # get the y arrays
    # sort the data before
    data = np.sort(data)
    y1 = edf(data,x)
    y2 = cdfBoltzmann(betaEst,x)
    
    plt.plot(x,y1)
    plt.plot(x,y2)
    plt.show()
    
    # calcualte the differences    
    delta = np.abs(y1-y2)
    
    # get the test variable
    D = np.max(delta)
    
    # threshold for 1% level    (no parameter estimation)
    #threshold = 1.63 / np.sqrt(N)
    
    # for estimated beta the threshold must be even smaller (1 % level)
    threshold = 1.25 / np.sqrt(N)     
     
    # test hypothesis
    if D > threshold:
        print "The two samples are probably not from the same distribution, chances are only 1%."
        result = "failed"
    else:
        print "No contradiction to hypothesis that the two samples belong to the same distribution."
        result = "passed"
    
    # ==== plot the data ==========
    
    # the histogram of the data
    n, bins, patches = plt.hist(data, 100, normed=1, facecolor='green', alpha=0.5)
    
    # add a 'best fit' line   
    theory = pdfBoltzmann(betaEst,bins)
    plt.plot(bins, theory, 'r--')
    plt.xlabel('Fitness')
    plt.ylabel('Probability density')
    plt.plot(x,delta/threshold,'b--')
    plt.savefig(plotname)
    
else:
    print errormessage