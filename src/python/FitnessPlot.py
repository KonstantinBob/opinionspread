import numpy as np
import matplotlib.pyplot as plt
import sys

# Setup Latex for Plots
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'font.family': 'lmodern'})
plt.rcParams.update({'font.size': 23})
plt.rcParams.update({'backend': 'PDF'})

filename = sys.argv[1]

# get the data
data= np.loadtxt(filename)

# legend etc
plt.title("fitness over time")
plt.xlabel("cycles")
plt.ylabel("fitness")

num_plots = len(data[1,:])

# Have a look at the colormaps here and decide which one you'd like:
# http://matplotlib.org/1.2.1/examples/pylab_examples/show_colormaps.html
colormap = plt.cm.rainbow
plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, num_plots)])

#perform the plot
for i in np.arange(num_plots):
	plt.plot(data[:,i],'o-',alpha=0.8)
	
plt.savefig(filename+".pdf")

plt.show()