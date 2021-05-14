import numpy as np
import matplotlib.pyplot as plt
import sys

binArray = 20 #np.linspace(5.0,7.0,50.)

# read the data from raw file
filename1 = sys.argv[1]
data1 = np.loadtxt(filename1+".dat",unpack=True)

plt.hist(data1,bins=binArray,alpha=0.7,normed=True,label=filename1)

# read the data from histo file
filename2 = sys.argv[2]
data2  = np.loadtxt(filename2+".dat",unpack=True)   

plt.hist(data2,bins=binArray,alpha=0.7,normed=True,label=filename2)

plt.legend(loc='best')

plt.savefig("CompareEffDiaDistr.pdf")
plt.show()