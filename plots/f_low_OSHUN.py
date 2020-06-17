#T_id    = '00200'
T_id    = '00000'
var_tit = 'pmulti'
path_id = '../OUTPUT/'
#path_id = '../OUTPUT_sm_dt/'

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
'''
contourf on 2D plane at a given time
'''

denor   = 0
'denormalize in def-Input'

plog = -1
'plog = 1 logspace levels, =0, linspace levels, =2 pcolormesh'
print('---   [',var_tit,'] reading time:', T_id,' ---')

'''start here'''
outputs = srd.rd_file(T_id, var_tit, path_id,denor)

title   = outputs[0]
x_title = outputs[1]
y_title = outputs[2]
t_out   = outputs[3]
num_x   = outputs[4]
num_y   = outputs[5]
x_axis  = outputs[6]
y_axis  = outputs[7]
data    = outputs[8]

print()
print('--- Dimension @time= %3.2f: ---'%(t_out) )

print('X(%i): [%6.2f,%6.2f], mid x: (%4.1e --- %4.1e)'%(\
        num_x,x_axis[0],x_axis[num_x-1],\
        min(data[:,int(num_x/2)]),max(data[:,int(num_x/2)]) ) )
print('Y(%i): [%6.2f,%6.2f], mid y: (%4.1e --- %4.1e)'%(\
        num_y,y_axis[0],y_axis[num_y-1],\
        min(data[int(num_y/2),:]),max(data[int(num_y/2),:]) ) )

'''plot here'''

'test fitting'
'x^n'
'x^n *exp(-x**2/B)'


fig2,axarr = plt.subplots(2,5,figsize=(19,9))
#axarr[0,0].plot(np.linspace(0,299,300),data[0,:])
axarr[0,0].plot(data[0,:])
axarr[0,0].set_title('f00.r')

#axarr[0,1].plot(np.linspace(0,299,300),data[2,:])
axarr[0,1].plot(data[2,:])
axarr[0,1].set_title('f10.r')

#axarr[0,2].plot(np.linspace(0,299,300),data[6,:])
axarr[0,2].plot(data[6,:])
axarr[0,2].set_title('f20.r')

#axarr[1,0].plot(np.linspace(0,299,300),data[12,:])
axarr[0,3].plot(data[12,:])
axarr[0,3].set_title('f30.r')

#axarr[0,4].plot(data[18,:])
axarr[0,4].plot(data[18,:])
axarr[0,4].set_title('f40.r')

plt.legend(loc=3,fontsize='small')

plt.show(block=False)


