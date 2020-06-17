import numpy as np
import scipy as sci
import scipy.special as sp
import cmath
from math import e
#put my modules here
import sub_rd_OSHUN as srd

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import ConnectionPatch
from matplotlib.gridspec import GridSpec
from matplotlib  import cm, colors
from mpl_toolkits.mplot3d import Axes3D
#for my mac...
import warnings
warnings.simplefilter(action="ignore",category= FutureWarning)
warnings.simplefilter(action="ignore",category= UserWarning)
#debug
import sys
import imp
'''
co-plot 1 variable at different time in one plot
'''

var_tit = 'N'
'denormalize in def-Input'
denor   = 0
'1D: direction x:4, y:5'
'plot the middle line'
direc_id = 5

########################################################################
path_id = '../../imass_noEx_noe/OUTPUT/'
T_list = ('00000','00030','00060','00100')
t_len = len(T_list)

for T_id in T_list:
    outputs = srd.rd_file(T_id, var_tit, path_id,denor)
    if(T_id == T_list[0]):
        var = np.zeros(shape=(len(T_list),outputs[direc_id]))
        tim = np.zeros(len(T_list))
        num = 0
        if(direc_id == 4):
            axis = outputs[7]
        if(direc_id == 5):
            axis = outputs[6]

    if(direc_id == 4):
        var[num,:]= outputs[8][:,int(outputs[5]/2)]
    if(direc_id == 5):
        var[num,:]= outputs[8][int(outputs[4]/2),:]
    tim[num] = outputs[3]
    num +=1

tim[1]='35893.'
tim[2]='71785.'
tim[3]='119642.'
'''plot'''
'plot1'
f = plt.figure(figsize=(15,15))
gs1 = GridSpec(3, 3)
gs1.update(left=0.05, right=0.95, wspace=0.25)

ax1 = plt.subplot(gs1[:-1, :])
for i in range(0,num):
    plt.plot(axis,var[i,:],label=tim[i])
plt.xlabel('x',size=20)
plt.ylabel('Density(x)',size=20)
plt.legend(loc=1)
plt.text(-750.0,0.9,'(a)',size=20)
########################################################################
path_id = '../OUTPUT_n3/'
T_list = ('00000','00030','00060','00100')
t_len = len(T_list)

for T_id in T_list:
    outputs = srd.rd_file(T_id, var_tit, path_id,denor)
    if(T_id == T_list[0]):
        var = np.zeros(shape=(len(T_list),outputs[direc_id]))
        tim = np.zeros(len(T_list))
        num = 0
        if(direc_id == 4):
            axis = outputs[7]
        if(direc_id == 5):
            axis = outputs[6]

    if(direc_id == 4):
        var[num,:]= outputs[8][:,int(outputs[5]/2)]
    if(direc_id == 5):
        var[num,:]= outputs[8][int(outputs[4]/2),:]
    tim[num] = outputs[3]
    num +=1

'plot2'
#fig = plt.figure(figsize=(15,5.5))
#ax2 = plt.subplot(1,2,1)
ax2 = plt.subplot(gs1[-1, :-1])
ax2.ticklabel_format(axis='x',style='sci',useOffset=False)
for i in range(0,num):
    plt.plot(axis,var[i,:],label=tim[i])
plt.legend(loc=2)
plt.xlabel('x',size=20)
plt.ylabel('Density(x)',size=20)
plt.text(35000,0.87,'(b)',size=20)

'plot3'
#ax3 = plt.subplot(1,2,2)
ax3 = plt.subplot(gs1[-1, -1])
ax3.ticklabel_format(axis='y',style='sci',useOffset=False)
for i in range(0,num):
    plt.plot(axis,var[i,:],label=tim[i])
#plt.legend()
plt.xlabel('x',size=20)
plt.xlim(671.0,684.0)
plt.ylim(0.9926+0.00005,0.9926+0.00042)
plt.text(682.5,0.99297,'(c)',size=20)

coordsA="data"
coordsB="data"

xy1 = (677,0.9926+0.00023)
xy2 = (671,0.9926+0.00025)
con = ConnectionPatch(xyA=xy2, xyB=xy1, coordsA=coordsA, coordsB=coordsB,
                      axesA=ax3, axesB=ax2,
                      arrowstyle="->", shrinkB=5)
ax3.add_artist(con)

plt.draw()
plt.show()
