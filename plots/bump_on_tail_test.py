import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib  import cm, colors
from mpl_toolkits.mplot3d import Axes3D
#put my modules here
import sub_rd_OSHUN as srd

#for my mac...
import warnings
warnings.simplefilter(action="ignore",category= FutureWarning)
warnings.simplefilter(action="ignore",category= UserWarning)
#debug
import sys
''' plot the distribution function '''

x = np.linspace(-0.12, 0.12,200)
vt  = 0.024
vb  = 0.065
vtb = 0.006
#f = np.exp(-x**2/vt**2 *0.5)+ x**2 * np.exp(-(x-vb)**2/vtb**2 *0.5)*100
f = np.exp(-x**2/vt**2 *0.5)+ x**2 * np.exp(-(x-vb)**2/vtb**2 *0.5)*100 + x**2 * np.exp(-(x+vb)**2/vtb**2 *0.5)*100
#f = x**2 * np.exp(-(x-vb)**2/vtb**2 *0.5)*0.2

fig = plt.figure()
plt.plot(x,f)
plt.show()

