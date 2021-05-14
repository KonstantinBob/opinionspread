#!/usr/bin/env python
"""
Script for fitting to Boltzmann distribution.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.integrate as integrate
import sys
import math as m

# ========= constants
errormessage = "usage: python "+sys.argv[0] + " data file"

# global variables 
argc = len(sys.argv)


# ========= functions    
def pot(x,a,b,c,d,e):    
    return a*pow(x,4.)+b*pow(x,3.)+c*pow(x,2.)+d*x+e
    
def u(x):
    return pot(x,160./3.,-104.,188./3.,-12.,1.0) 
    
def phi(beta):
    return integrate.quad(lambda x:np.exp(-beta*u(x)) ,0.,1.)[0]
    
def psi(beta):
    return integrate.quad(lambda x: u(x)*np.exp(-beta*u(x)) ,0.,1.)[0]/integrate.quad(lambda x:np.exp(-beta*u(x)) ,0.,1.)[0]
    
def theory(betaEst,x):
    return np.exp(-betaEst*u(x))/phi(betaEst)
    


# ========= main =======================
# parse the arguments
if argc == 2:    
    
    # read the file 
    filename = sys.argv[1]
    coord = np.loadtxt(filename,unpack = True)
    
    # ===== data processing =======
    
    # calculate the sum and size of sample
    data = u(coord)
    sumData = np.sum(data)
    lengthData = np.size(data)
    
    # possible beta values
    X = np.linspace(0.,10.,num=50)   
    
    # neg log likelihood
    Y = np.log(np.vectorize(phi)(X))+X*sumData/lengthData
    
    # function with root at beta estimate
    Z = np.vectorize(psi)(X)-sumData/lengthData
    
    # ===== parameter estimation =======
    # minimum value of neg log likelihood
    Ymin = np.argmin(Y)
    
    # beta estimate
    betaEstimate = X[np.argmin(Y)]    
    print betaEstimate
       
    # ===== errors of estimate
    # left error        
    #betaErrorleft, succLeft = opt.brentq(lambda x:lengthData*np.log(phi(x))+x*sumData+0.5+Ymin ,0.,1.0001*betaEstimate)    
    #betaErrorleft, succLeft = opt.brentq(lambda x: phi(x)*sumData+0.5+Ymin ,0.,1.0001*betaEstimate)    
    # right error
    #betaErrorRight, succRight = opt.brentq(lambda x: lengthData*m.log(phi(x))+x*sumData+0.5+Ymin,0.9999*betaEstimate,10)    
    #print betaEstimate, " + " ,betaErrorRight, " - ", betaErrorleft
    
    
    # show the neg log likelihood plot
    plt.plot(X,Y)
    plt.show()   
        
    # show data and fit
    n, bins, patches = plt.hist(coord, 50, normed=1, facecolor='blue', alpha=0.5,label="data")
    plt.plot(bins,theory(betaEstimate,bins),'r',label="theory pdf")
    plt.plot(bins,u(bins),'g',label="potential")
    plt.legend(loc="best",numpoints=1)
    plt.title(filename+" -  beta = "+"{0:.3g}".format(betaEstimate))
    plt.savefig("boltzmann"+filename+".png")
    plt.show()
    
        
else:
    print errormessage
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    