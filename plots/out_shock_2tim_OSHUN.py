import numpy as np
import scipy as sci
import scipy.special as sp
import cmath
from math import e
#put my modules here
import sub_rd_OSHUN as srd

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib  import cm, colors, gridspec
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
#for my mac...
import warnings
warnings.simplefilter(action="ignore",category= FutureWarning)
warnings.simplefilter(action="ignore",category= UserWarning)
#debug
import sys
import imp
'''
multiplot variables using the same axes
'''

import pickle

fileObject = open('save_t0.p','rb')
d0 = pickle.load(fileObject)
fileObject = open('save_t2.p','rb')
d2 = pickle.load(fileObject)

title1  = d0[0]
x_title = d0[1]
y_tit1  = d0[2]
t_outF  = d0[3]
num_x   = d0[4]
num_y   = d0[5]
x_axis  = d0[6]
y_axis1 = d0[7]
data1F  = d0[8]

title2  = d0[9]
y_tit2  = d0[10]
data2F  = d0[11]

title3  = d0[12]
y_tit3  = d0[13]
y_axis3 = d0[14]
data3F  = d0[15]

title4  = d0[16]
y_tit4  = d0[17]
y_axis4 = d0[18]
data4F  = d0[19]

title5  = d0[20]
x_tit5  = d0[21]
x_axis5 = d0[22]
y_axis5 = d0[23]
data5F  = d0[24]
bar5lF  = d0[25]
bar5hF  = d0[26]

t_outL  = d2[3]
data1L  = d2[8]
data2L  = d2[11]
data3L  = d2[15]
data4L  = d2[19]
data5L  = d2[24]
bar5lL  = d2[25]
bar5hL  = d2[26]


range_u = 1.0 + 0.3
range_b = 2.0 - range_u


f, ((ax1,ax2), (ax3,ax4), (ax5,ax6), (ax7,ax8) ) = plt.subplots(4,2, sharey='row',sharex='col')

ax1.set_title('t= %.2e'%t_outF)
ax1.set_ylabel(title1)
ax1.plot(x_axis,data1F[int(num_y/2),:])
y2 = max((max(data1F[int(num_y/2),:]),max(data1L[int(num_y/2),:])))
y1 = min((min(data1F[int(num_y/2),:]),min(data1L[int(num_y/2),:])))
ax1.set_xlim([-512, 512])
ax1.set_ylim([y1*range_b, y2*range_u])
ax1.xaxis.set_major_formatter(ticker.FuncFormatter(srd.fmt))
ax1.yaxis.set_major_formatter(ticker.FuncFormatter(srd.fmt))
ax1.locator_params(nbins=6,axis='x')
ax1.locator_params(nbins=7,axis='y')
#ax1.legend(loc=3)

ax2.set_title('t= %.2e'%t_outL)
ax2.plot(x_axis,data1L[int(num_y/2),:])
ax2.set_xlim([-512, 512])
ax2.xaxis.set_major_formatter(ticker.FuncFormatter(srd.fmt))
ax2.locator_params(nbins=6,axis='x')

#---
ax3.set_ylabel(title2)
ax3.plot(x_axis,data2F[int(num_y/2),:])
y2 = max((max(data2F[int(num_y/2),:]),max(data2L[int(num_y/2),:])))
y1 = min((min(data2F[int(num_y/2),:]),min(data2L[int(num_y/2),:])))
ax3.set_ylim([y1*range_b, y2*range_u])
ax3.yaxis.set_major_formatter(ticker.FuncFormatter(srd.fmt))
ax3.locator_params(nbins=4,axis='y')

ax4.plot(x_axis,data2L[int(num_y/2),:])

#---
ax5.set_ylabel(title3)
ax5.plot(x_axis,data3F[int(num_y/2),:])
y2 = max((max(data3F[int(num_y/2),:]),max(data3L[int(num_y/2),:])))
y1 = min((min(data3F[int(num_y/2),:]),min(data3L[int(num_y/2),:])))
# here y1<0
ax5.set_ylim([y1*range_u, y2*range_u])
ax5.yaxis.set_major_formatter(ticker.FuncFormatter(srd.fmt))
ax5.locator_params(nbins=6,axis='y')

ax6.plot(x_axis,data3L[int(num_y/2),:])

#---
ax7.set_ylabel(title4)
ax7.set_xlabel(x_title)
ax7.plot(x_axis,data4F[int(num_y/2),:])
y2 = max((max(data4F[int(num_y/2),:]),max(data4L[int(num_y/2),:])))
y1 = min((min(data4F[int(num_y/2),:]),min(data4L[int(num_y/2),:])))
ax7.set_ylim([y1*range_b, y2*range_u])
ax7.yaxis.set_major_formatter(ticker.FuncFormatter(srd.fmt))
ax7.locator_params(nbins=7,axis='y')

ax8.set_xlabel(x_title)
ax8.plot(x_axis,data4L[int(num_y/2),:])

'''
#f, ((ax9,ax10)) = plt.subplots(1,2, sharey='row')
f = plt.figure()
gs = gridspec.GridSpec(1, 2, width_ratios=[1,1])
#---
#bar5h = np.max((bar5hF, bar5hL))
#bar5l = np.min((bar5lF, bar5lL))
bar5h = bar5hF
bar5l = bar5lF
levels     = np.linspace(bar5l,bar5h,300)

ax9 = plt.subplot(gs[0])
first = ax9.contourf(y_axis5, x_axis5, data5F.transpose(),cmap='jet', levels=levels)
cbar  = f.colorbar(first, orientation='horizontal', format=mpl.ticker.FuncFormatter(srd.fmt))
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.update_ticks()
ax9.set_xlabel(x_title)
ax9.set_ylabel(x_tit5)
ax9.locator_params(nbins=7,axis='y')
#-
bar5h = bar5hL
bar5l = bar5lL
levels     = np.linspace(bar5l,bar5h,300)

ax10 = plt.subplot(gs[1])
second = ax10.contourf(y_axis5, x_axis5, data5L.transpose(),cmap='jet', levels=levels)
cbar  = f.colorbar(second, orientation='horizontal', format=mpl.ticker.FuncFormatter(srd.fmt))
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.update_ticks()
ax10.set_xlabel(x_title)
ax10.locator_params(nbins=7,axis='y')
'''
#---
plt.show()
