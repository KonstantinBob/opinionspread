#!/usr/bin/env python

# imports
import numpy as np
import matplotlib.pyplot as plt
import sys
from  scipy.special import zeta
from  scipy.optimize import brent,brentq 


# Setup Latex for Plots
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'font.family': 'lmodern'})
plt.rcParams.update({'font.size': 23})
plt.rcParams.update({'backend': 'PDF'})

# ========= constants
errormessage = "usage: "+sys.argv[0]+" inputfile  plotname xmin fitFlag"

# global variables 
argc = len(sys.argv)

# ========= functions
def powerlaw(x,expo,xmin,scale):
    return scale*np.power(x,-1.0*expo)/zeta(expo,xmin)

def negLogLike(gamma,sumLog,n,xmin):
    return  n * np.log(zeta(gamma,xmin)) + sumLog*gamma

def fitPowerLaw(fitData,xmin):
    # store only x > xmin
    data = fitData[np.where( fitData > xmin)]
           
    # store sum of log data
    sumLogX = np.sum(np.log(data))
    
    # store data length
    n = np.size(data)
    
    if False:    
        gammas = np.linspace(1.1,4.0)
        plt.plot(gammas,negLogLike(gammas,sumLogX,n,xmin))    
        plt.show()
        
    # minimize negloglikelihood
    gammaHat = brent(lambda x : negLogLike(x,sumLogX,n,xmin),brack=(1.1,4.0))
    minimumValue = negLogLike(gammaHat,sumLogX,n,xmin)
        
    # calculate error left
    boundLeft = 0.9 * gammaHat
    gammaLeft = brentq(lambda x : minimumValue + 0.5 - negLogLike(x,sumLogX,n,xmin),boundLeft,gammaHat)
    errorLeft = gammaHat-gammaLeft
    
    # calculate error right
    boundRight = 1.1 * gammaHat
    gammaRight = brentq(lambda x : minimumValue + 0.5 - negLogLike(x,sumLogX,n,xmin),boundRight,gammaHat)
    errorRight = gammaRight-gammaHat
    
    return gammaHat,errorLeft,errorRight
    

# ========= main =======================
if argc == 5:
    
    # get the arguments
    filename = sys.argv[1]
    plotname = sys.argv[2]
    xminArg = float(sys.argv[3])
    fitFlag = sys.argv[4]
        
    # check if fitting is to be done
    doFit = False  
    if fitFlag == "y":
        doFit = True 
    
    # get the legend
    file = open(filename)
    title = file.readline().replace('\n', '')
    xLabel = file.readline().replace('\n', '')
    yLabel = file.readline().replace('\n', '')
    plotname = file.readline().replace('\n', '')
    file.close() 
    
    # get the data
    x,y = np.loadtxt(filename,unpack=True,skiprows=4)
    scale = np.sum(y)
      
    # prepare fit data
    fitDataList = list()
    for index in np.arange(len(x)):
        for counter in np.arange(y[index]):
            fitDataList.append(x[index])
    
    fitData = np.array(fitDataList,dtype='float')
       
    # do fit    
    xmin =  x[np.argmax(y)]
    if xminArg != 0:
        xmin = xminArg
        
    exponent,errorLeft,errorRight = fitPowerLaw(fitData,xmin)
    
    print  exponent,errorLeft,errorRight,xmin
            
    # ====== perform the plot 
    plt.title(title)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
        
    plt.xlim([1,1e4])
    plt.ylim([1,1e4])
    #plt.ticklabel_format(style='sci',axis='x',scilimits=(0,0))
    
    fitX = np.arange(xmin,1e4) 
    
    plt.loglog(x,y,"x",label="Data", zorder = 1)
    if doFit:
        plt.loglog(fitX,powerlaw(fitX,exponent,xmin,scale),label="Fit",lw=2,zorder = 2)
        plt.axvline(x=xmin,ls='--',alpha=0.5,label='min. Degree')
    
    ax = plt.subplot(111)
    ax.spines["top"].set_visible(False) 
    ax.spines["right"].set_visible(False)  
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    legend=plt.legend(loc='best',numpoints=1,frameon=False)    
    plt.setp(plt.gca().get_legend().get_texts(), fontsize='18') #legend 'list' fontsize
    plt.tight_layout()
   
    plt.savefig(plotname+".pdf", bbox_inches='tight')
    plt.show()
    
else:
    print errormessage
