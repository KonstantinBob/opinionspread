"""
Script for estimating mean and variance of a normal distributed sample.
Errors for the estimators will also be calculated.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

# ========= constants
errormessage = "usage: python "+sys.argv[0]+" inputfile  outputfile"

# global variables 
argc = len(sys.argv)

# ========= functions

# ========= main =======================
# parse the arguments
if argc == 3 :
    
    # find out the deszired mode
       
    # read the data from  inputfile and sort it 
    inputFilename = sys.argv[1]
    data = np.sort(np.loadtxt(inputFilename+".dat",unpack=True,dtype=np.float64))
    n = float(data.size)
            
    # get the name of the outputfile
    outputFilename = sys.argv[2]
    
    # calculate the estimators according to Blobel and Lohrmann
    sumEntries = 0
    sumSqauredEntries = 0
    estimate = np.median(data)
    
    for entry in data:
        sumEntries = sumEntries + (entry-estimate)
        sumSqauredEntries = sumSqauredEntries + (entry-estimate)**2.0
        
    meanEstimator = estimate + (sumEntries)/(n)
    varianceEstimator = (sumSqauredEntries-(sumEntries**2.0)/(n))/(n-1.0) 
    
    # caluclate the errors aka standard deviations of the estimators 
    # ---> ONLY VALID FOR NORMAL DISTRIBUTED DATA <----
    meanError = np.sqrt((varianceEstimator)/(n)) 
    varianceError = np.sqrt((varianceEstimator)/(2.0*(n-1)))
    
    print meanEstimator, meanError
    print varianceEstimator, varianceError     

else:
    print errormessage
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    