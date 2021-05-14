import numpy as np
import matplotlib.pyplot as plt
import sys

# Setup Latex for Plots
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'font.family': 'lmodern'})
plt.rcParams.update({'font.size': 23})
plt.rcParams.update({'backend': 'PDF'})

filename = sys.argv[1]

# get the legend
file = open(filename)
title = file.readline().replace('\n', '')
xLabel = file.readline().replace('\n', '')
yLabel = file.readline().replace('\n', '')
plotname = file.readline().replace('\n', '')
file.close()

# get the data
x,y = np.loadtxt(filename,unpack=True,skiprows=4)

# do fit?
#plt.xlim([1e2,1e5])
#plt.ylim([1e2,1e6])

#perform the plot
plt.title(title)
plt.xlabel(xLabel)
plt.ylabel(yLabel)
plt.loglog(x,y,'ro')
plt.tight_layout()
plt.savefig(plotname)
plt.show()