T_id    = '00001'
var_tit = 'fsp'
path_id = '../OUTPUT/'

import numpy as np
import scipy as sci
import scipy.special as sp
import cmath
from math import e
#put my modules here
import sub_rd_OSHUN as srd

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib  import cm, colors
from mpl_toolkits.mplot3d import Axes3D
#for my mac...
import warnings
warnings.simplefilter(action="ignore",category= FutureWarning)
warnings.simplefilter(action="ignore",category= UserWarning)
#debug
import sys
import imp
''' for o_fsp '''

'''start here'''
outputs = srd.rd_file(T_id, var_tit, path_id,0)
title   = outputs[0]
x_title = outputs[1]
y_title = outputs[2]
t_out   = outputs[3]
num_x   = outputs[4]
num_y   = outputs[5]
x_axis  = outputs[6]
y_axis  = outputs[7]
data0   = outputs[8]

outputs = srd.rd_file(T_id, var_tit, path_id,1)
data1   = outputs[8]

outputs = srd.rd_file(T_id, var_tit, path_id,2)
data2   = outputs[8]

outputs = srd.rd_file(T_id, var_tit, path_id,3)
data3   = outputs[8]

outputs = srd.rd_file(T_id, var_tit, path_id,4)
data4   = outputs[8]

outputs = srd.rd_file(T_id, var_tit, path_id,5)
data5   = outputs[8]

'''plots'''
fig = plt.figure(figsize=(18.,11))

plt.subplot(421)
plt.plot(x_axis,data0[int(num_y/2),:],label='f0(p=199)')
plt.legend()

plt.subplot(423)
plt.plot(x_axis,data1[int(num_y/2),:],label='f1(p=199)')
plt.legend()

plt.subplot(425)
plt.plot(x_axis,data2[int(num_y/2),:],label='f2(p=199)')
plt.legend()

plt.subplot(427)
plt.plot(x_axis,data0[int(num_y/2),:],label='lp 0')
plt.plot(x_axis,np.abs(data1[int(num_y/2),:]),label='lp 1')
plt.plot(x_axis,np.abs(data2[int(num_y/2),:]),label='lp 2')
plt.legend()
#
plt.subplot(422)
plt.plot(x_axis,data3[int(num_y/2),:],label='f0(p=9)')
plt.legend()

plt.subplot(424)
plt.plot(x_axis,data4[int(num_y/2),:],label='f1(p=9)')
plt.legend()

plt.subplot(426)
plt.plot(x_axis,data5[int(num_y/2),:],label='f2(p=9)')
plt.legend()

plt.subplot(428)
plt.plot(x_axis,data3[int(num_y/2),:],label='sp 0')
plt.plot(x_axis,np.abs(data4[int(num_y/2),:]),label='sp 1')
plt.plot(x_axis,np.abs(data5[int(num_y/2),:]),label='sp 2')
plt.legend()





plt.xlabel(x_title)
plt.show()
