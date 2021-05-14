#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 14:16:19 2017

@author: konbob
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 10:04:37 2017

@author: konbob
"""

import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.signal import find_peaks_cwt,argrelmax,hann,convolve
import pandas as pd
import cPickle as pickle
from collections import OrderedDict

# Setup Latex for Plots
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'font.family': 'lmodern'})
plt.rcParams.update({'font.size': 26})
plt.rcParams.update({'backend': 'PDF'})
  
# ========= constants =
errormessage = "usage: python "+sys.argv[0] + " peakfile  debug flag (-d for true, -n for false) simulationFiles"

# global variables 
argc = len(sys.argv)

# ========= functions
def findPeaksByDensity(dataSet,debug,p,q,size,net):
    # append zeros at beginning and end to find peaks at the boundaries
    n = np.ceil(0.02*size)
    empty = np.zeros(int(n))
    data = np.concatenate((empty,dataSet,empty))
    
    # calculate fraction up
    x = np.arange(0-n,size+1+n)/(float(size+2.0*n))
              
    # smooth the data
    win = hann(15)    
    filtered =convolve(data, win, mode='same') / sum(win)
        
    # find maxima by wavlets    
    maxima =  find_peaks_cwt(filtered,0.5*size*np.arange(.10,.20,step=0.02))
    
    # if that failes, use the argelmax
    if len(maxima) == 0:
        maxima = argrelmax(data)[0]
        print 'Used argrelmax for p',p,"q ",q,"size ",size  
                          
    # drop irrelevant maxima
    newMax = list()
    maxPMF = np.max(data)        
    for maxi in maxima:        
        if  data[maxi]> 1e-2 * maxPMF:
            newMax.append(maxi)
    if len(newMax) > 2:
        newMax.pop(1)
        
    #plot single hist
    ax = plt.subplot(111)
    ax.spines["top"].set_visible(False) 
    ax.spines["right"].set_visible(False)  
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    
    plt.plot(x,data,label="PDF") 
    plt.plot(x,filtered,ls='--',label="smooth")        
            
    # calculate the errors
    errorWidth = 5./size
    errors = errorWidth*np.ones_like(x[newMax])
    
    for maxi in newMax:
        plotPeakLine = plt.axvline(x=x[maxi],ls='--',color='g',label='Peak')
        xLeft = x[maxi]-errorWidth        
        xRight = x[maxi]+errorWidth        
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
    size = float(values[1])
    p = float(values[2])
    q = float(values[3])  
    
     # print information
    if debug: 
        print "processing: p ",p,"q ",q,"size ",size       
        
    # store Data    
    df = pd.read_csv(filename,dtype='float',comment='#',header=None)    
    timeSeries = df.as_matrix()[:,1]
                 
    peaks,errors = findPeaksByDensity(timeSeries,debug,p,q,size,net)           
    results[(p,q,size)] = [peaks,errors]
          
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