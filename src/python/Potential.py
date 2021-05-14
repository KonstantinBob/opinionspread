#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0.0,1.0,100)

def pot(x,a,b,c,d,e):
    print a,b,c,d,e
    return a*np.power(x,4.)+b*np.power(x,3.)+c*np.power(x,2.)+d*x+e
    
print x, pot(x,160./3.,-104.,188./3.,-12.,1.0)    
plt.plot(x,pot(x,160./3.,-104.,188./3.,-12.,1.0))
plt.show()    

