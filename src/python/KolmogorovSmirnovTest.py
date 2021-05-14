#!/usr/bin/env python
"""
Script for performing a Kolmogoro Smirnov test.
It can compare to two samples or one sample and a distribution.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stat
import sys

# ========= constants
errormessage = "usage: python "+sys.argv[0]+" -flag plotname datafile1 (datafile2) \n flags: \n -s two sample mode \n -n normal distribution mode \n -e exponential distribution mode \n -u  uniform distribution mode"

# global variables 
argc = len(sys.argv)
normalDistrMode=twoSampleMode=False;

# ========= functions

#empirical distribution function
def edf(dataSet,value):
    return np.searchsorted(dataSet, value, side='right') / float(dataSet.size)
    
def cdfNormal(mu,sigma,value):
    return stat.norm.cdf(value,loc=mu,scale=sigma)        
    
def pdfNormal(mu,sigma,value):
    return stat.norm.pdf(value,loc=mu,scale=sigma)
    
def cdfUni(value):
    return stat.uniform.cdf(value)  
    
def pdfUni(value):
    return stat.uniform.pdf(value)  
    
def cdfBoltzmann(beta,value):
    return stat.expon.cdf(value,scale=1./beta)
    
def pdfBoltzmann(beta,value):
    return stat.expon.pdf(value,scale=1./beta)

# ========= main =======================
# parse the arguments
if argc == 4 or argc == 5:
    
    # find out the deszired mode
    flag = sys.argv[1]
    
    # get the name of the plot
    plotname = sys.argv[2]
    
    # prepare the flags
    twoSampleMode = False
    normalDistrMode = False
    exponentialDistrMode = False
    uniformDistrMode = False
    
    # set the flag
    if flag == "-s":
        twoSampleMode = True
    elif flag == "-n":
        normalDistrMode = True
    elif flag == "-e":
        exponentialDistrMode = True
    elif flag == "-u":
        uniformDistrMode = True
    else:
        print errormessage
        exit
    
    # read the data from  file 1 and sort it 
    filename1 = sys.argv[3]
    data1 = np.sort(np.loadtxt(filename1+".dat",unpack=True))
    n1 = float(data1.size)
    
    # for later use
    theory= None
    
    if twoSampleMode==True:    
        # read the data from file 2 
        filename2 = sys.argv[4]
        data2  =  np.sort(np.loadtxt(filename2+".dat",unpack=True))  
        n2 = float(data2.size)           
    
    # calculate the x values
    if twoSampleMode:
        x_min = min(np.min(data1),np.min(data2))
        x_max = max(np.max(data1),np.max(data2))
    if normalDistrMode or exponentialDistrMode or uniformDistrMode:
        x_min = np.min(data1)
        x_max = np.max(data1)
        
    # get the x array    
    x = np.linspace(x_min,x_max,500)

    # get the y arrays    
    y1 = edf(data1,x)
    
    if twoSampleMode == True:
        y2 = edf(data2,x)
        
    if normalDistrMode == True:
        # calculate the estimators
        mu = np.mean(data1)
        sigma = np.std(data1)
        
        # prepare data set from cdf
        y2 = cdfNormal(mu,sigma,x)
        
        # prepare theory 
        theory = pdfNormal(mu,sigma,x)
        
        
    if exponentialDistrMode == True:
        # calculate the estimators        
        sumData = np.sum(data1)
    
        # calculate the mean of param beta
        betaEst = n1/sumData 
    
        # calculate the error of param estimator
        betaError = np.sqrt(n1)/sumData
        
        # prepare the cdf
        y2 = cdfBoltzmann(betaEst,x)
        
        # prepare theory 
        theory = pdfBoltzmann(betaEst,x)
        
    if uniformDistrMode == True:
        # prepare the cdf
        y2 = cdfUni(x)
        
        # prepare theory 
        theory = pdfUni(x)
        
    # calcualte the differences    
    delta = np.abs(y1-y2)
    
    # get the test variable
    D = np.max(delta)
    
    # threshold for 1% level
    if twoSampleMode == True:
        threshold = 1.63 * np.sqrt((n1+n2)/(n1*n2))
    if normalDistrMode == True:
        threshold = 1.031 * np.sqrt(1.0/n1) 
    if exponentialDistrMode == True:
        threshold = 1.25 * np.sqrt(1.0/n1)
    if uniformDistrMode == True:
        threshold = 1.63 * np.sqrt(1.0/n1) 
    
    # test hypothessis
    if D > threshold:
        print "The two samples are probably not from the same distribution, chances are only 1%."
        result = "failed"
    else:
        print "No contradiction to hypothesis that the two samples belong to the same distribution."
        result = "passed"
    
    # Setup Latex for Plots
    plt.rcParams.update({'text.usetex': True})
    plt.rcParams.update({'font.family': 'lmodern'})
    plt.rcParams.update({'font.size': 13})
    plt.rcParams.update({'backend': 'PDF'})
    
    plt.axes().set_aspect('equal', 'datalim')   
    
    # ======= plot the results ===========
    # cdf and test statistic    
    ax = plt.subplot(2, 1, 1)
    
    plt.plot(x,y1,'b',label='data')
    
    if twoSampleMode == True:
        plt.plot(x,y2,'g--',label=filename2)
    if normalDistrMode == True:
        plt.plot(x,y2,'g--',label="fitted distribution")
    if exponentialDistrMode == True:
        plt.plot(x,y2,'g--',label="fitted distribution")
    if uniformDistrMode == True:
        plt.plot(x,y2,'g--',label="predetermined distribution")
        
    plt.plot(x,delta/threshold,'r--',label="test statistic")
        
    plt.ylim([0,1])
    plt.xlabel("data")
    plt.ylabel("cumulative probability")
    plt.title("Kolmogorov-Smirnov-Test of " + filename1 + ". Result: " + result)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),ncol=3)      
    
    # pdf and histogram   
    plt.subplot(2, 1, 2)
    n, bins, patches =   plt.hist(data1, 100, normed=1, facecolor='blue', alpha=0.5)  
    if not twoSampleMode:
        plt.plot(x, theory,'g--')
    else:
        plt.hist(data2, 100, normed=1, facecolor='blue', alpha=0.5) 
        
    plt.xlabel('data')
    plt.ylabel('Probability density')
     
    # save fig           
    plt.savefig("KS-Test_"+plotname+".png")
    plt.show() 
    
else:
    print errormessage
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    