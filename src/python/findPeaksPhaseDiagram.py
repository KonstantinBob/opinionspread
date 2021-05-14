#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 10:04:37 2017

@author: konbob
"""

import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy import stats
from scipy.signal import find_peaks_cwt
import pandas as pd
import cPickle as pickle
from collections import OrderedDict

# Setup Latex for Plots
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'font.family': 'lmodern'})
plt.rcParams.update({'font.size': 26})
plt.rcParams.update({'backend': 'PDF'})
params = {'text.latex.preamble' : [r'\usepackage{amsmath}',r'\usepackage{amssymb}']}
plt.rcParams.update(params)
  
# ========= constants =
errormessage = "usage: python "+sys.argv[0] + " peakfile  debug flag (-d for true, -n for false) simulationFiles"

# global variables 
argc = len(sys.argv)



# ========= functions
def baseline(x):
    return 8*x+1

def getRawMoment(dataset,order):
    # sort data set
    np.ndarray.sort(dataset)
    
    # prepare the sum'
    res = 0.
    
    # compute powers
    for element in dataset:       
        res+= np.power(element,order)
        
    # normalize and return
    return res/(len(dataset))

# ========= compute the Binder cumulant
def getBinderCumulant(dataSet):     
    
    # compute forth raw moment of shifted data
    mom4 = getRawMoment(dataSet,4)
    
    # compute the raw central moment  of shifted data
    mom2 = getRawMoment(dataSet,2)
    
    # return the cumulant -> non central cumulant
    return mom4/np.power(mom2,2.0)


def findPeaksByDensity(dataSet,debug,NBins,p,q,size,h,net):
    try:
        # do density estimation with large bandwith
        bandWidth = 2. / np.power(len(dataSet),0.2)
        density = stats.gaussian_kde(dataSet,bw_method = bandWidth)
        
        # get the support of the data set
        epsilon = 1e-2
        Nsupport = 80
        x = np.linspace(np.min(dataSet)-epsilon,np.max(dataSet)+epsilon,num=Nsupport)
        
        #get smooth representation
        smooth = density.pdf(x)
        
        # find maxima by wavlets
        maxima = find_peaks_cwt(smooth,0.5*Nsupport*np.arange(.15,.65,step=0.05))
                  
                
        # drop irrelevant maxima
        newMax = list()
        maxPMF = np.max(smooth)        
        for maxi in maxima:        
            if  density.pdf(x[maxi])> 5e-2 * maxPMF:
                newMax.append(maxi)
        if len(newMax) > 2:
            newMax.pop(1)
            
        #plot single hist
        ax = plt.subplot(111)
        ax.spines["top"].set_visible(False) 
        ax.spines["right"].set_visible(False)  
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        
        plt.hist(timeSeries,bins=NBins,normed=1, alpha=0.5,label="Data")
        plt.plot(x,smooth,label="KDE")        
                
        # calculate the errors
        errorWidth = 4*(np.max(dataSet)-np.min(dataSet))/Nsupport
        errors = errorWidth*np.ones_like(x[newMax])
        
        for maxi in newMax:
            plotPeakLine = plt.axvline(x=x[maxi],ls='--',color='g',label='Peak')
            xLeft = x[maxi]-errorWidth
            #plt.axvline(x=xLeft,ls='--',color =  plotPeakLine.get_color(), alpha=0.5)
            xRight = x[maxi]+errorWidth
            #plt.axvline(x=xRight,color =  plotPeakLine.get_color(), alpha=0.5)
            plt.axvspan(xLeft,xRight,color =  plotPeakLine.get_color(), alpha=0.3)
           
            
        plt.xlabel('$\\frac{N_{\\uparrow}}{N}$')
        plt.ylabel('pmf') 
        
        # make nice legend
        handles, labels = plt.gca().get_legend_handles_labels()    
        by_label = OrderedDict(zip(labels, handles))    
        legend=plt.legend(by_label.values(), by_label.keys(),numpoints=1,frameon=False,ncol=3,loc='upper center', bbox_to_anchor=(0.5, -0.25))
        legend.get_title().set_fontsize('21') 
        plt.setp(plt.gca().get_legend().get_texts(), fontsize='26')
        
        if debug:                   
            plt.show()
        else:
            plt.savefig(net+'_'+str(size)+'_'+str(p)+'_'+str(q)+'.pdf',bbox_inches='tight')
            plt.title('p = '+str(p)+' q = '+str(q)+' size = '+str(size))
            plt.savefig(net+'_'+str(size)+'_'+str(p)+'_'+str(q)+'.png',bbox_inches='tight')
            plt.close()
                          
        # get the values
        return x[newMax],errors
    except: 
        #plot single hist
        ax = plt.subplot(111)
        ax.spines["top"].set_visible(False) 
        ax.spines["right"].set_visible(False)  
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        
        plt.hist(timeSeries,bins=NBins,normed=1, alpha=0.5)              
                   
        plt.xlabel('$\\frac{N_{\\uparrow}}{N}$')
        plt.ylabel('pmf')    
        plt.title(net+ ' p = '+str(p)+' q = '+str(q)+' size = '+str(size))
        if debug:                   
            plt.show()
        else:
            plt.savefig(str(size)+'_'+str(p)+'_'+str(q)+'.png', bbox_inches='tight')
            plt.close() 
        raise Exception

# ========= main =======================
# *                                    *
# *                                    *
# ======================================

#  ========= parse the arguments ========= 
if argc < 3:    
    print errormessage
    sys.exit()
    
# store filename for storage   
filenameStore = str(sys.argv[1])+".pickle"

#prepare map for results
results = {}
cumulants = {}

#debugging
debug = False
if str(sys.argv[2]) == '-d':
    debug = True
elif str(sys.argv[2]) != '-n':
    print 'Could not parse ',str(sys.argv[2])
    print errormessage
    sys.exit()

#  =========  file loop ========= 
for filename in sys.argv[3:]:
           
    # get the values of q and e 
    values = filename.split("_")
    net = str(values[0])
    size = values[1]
    p = float(values[2])
    q = float(values[3])  
    
     # print information
    if debug: 
        print "processing: p ",p,"q ",q,"size ",size       
        
    # store Data    
    df=pd.read_csv(filename,dtype='float',comment='#')
    timeSeries = df.as_matrix()[:,0]
        
    # calculate the number of bins
    h = 3.49 * np.std(timeSeries)/np.cbrt(len(timeSeries))
    if h != 0:    
        NBins =int( np.ceil((np.max(timeSeries)-np.min(timeSeries))/h))
                                    
        # analysis
        try:    
            peaks,errors = findPeaksByDensity(timeSeries,debug,NBins,p,q,size,h,net)
            
            if len(peaks) == 0:                                       
                print 'No peak for p ',p,"q ",q,"size ",size,' - Using other method!'
                hist,bins = np.histogram(timeSeries,bins=NBins)
                binWidth = bins[1]-bins[0]
                peaks = np.min(timeSeries)+binWidth*np.argmax(hist)
                errors = 1e-1
                
            else:       
                # store results
                results[(p,q,size)] = [peaks,errors]
                cumulants[(p,q,size)] = getBinderCumulant(timeSeries)
        except Exception as e:
            print "Error: ",e
    else:
        unique = np.unique(timeSeries) 
        if len(unique) == 1:
            peaks,errors = unique,1e-1
        else:
            print 'Singular values for p ',p,"q ",q,"size ",size  
    
# ======= save the results
filePeaks = open( filenameStore, "wb" )
pickle.dump( results, filePeaks )
filePeaks.close()

fileCumulants = open( str(sys.argv[1])+'Cumulants', "wb" )
if True:
    pickle.dump(cumulants,fileCumulants)
if False:
    for cumu in cumulants:
        p,q,size = cumu
        value = cumulants[cumu]
        fileCumulants.write(str(p)+' '+str(q)+' '+str(size)+" "+str(value)+'\n')
fileCumulants.close()