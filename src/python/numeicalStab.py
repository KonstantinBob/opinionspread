#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 29 11:08:31 2017

@author: konbob
"""

import numpy as np

sampleSize = 100000

width=5.

# generate numbers
negNumbers = np.random.uniform(-10-width,-10+width,sampleSize)

posNumbers = np.random.uniform(10-width,10+width,sampleSize)

nullNumbers = np.random.uniform(-width,width,sampleSize)

mN = np.mean(negNumbers)

mP= np.mean(posNumbers)

mZero = np.mean(nullNumbers)

divisionMean = (mN+mP+mZero)/3.0

print divisionMean

longArray = np.append(np.append(negNumbers, posNumbers),nullNumbers)
np.random.shuffle(longArray)

mean=0.

for x in longArray:
    mean+=x

mean/=len(longArray)

print mean, divisionMean/mean, np.mean(longArray)