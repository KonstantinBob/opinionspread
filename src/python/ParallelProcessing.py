#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 14:44:46 2017

@author: konbob
"""
from functools import partial
from multiprocessing import Pool

def f(y,x):
    
    return y*x*x

fp = partial(f,1) 

print fp(2)

if __name__ == '__main__':
    pool = Pool(processes=4)              # start 4 worker processes

    squares =  pool.map(fp, range(10))

    print squares