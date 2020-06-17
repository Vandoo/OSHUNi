T_id    = '00000'
#T_id    = '00001'
var_tit1 = 'N'
var_tit2 = 'Vsq'
var_tit3 = 'Vx'
var_tit4 = 'Ex'
var_tit5 = 'p1x1'
path_id = '../OUTPUT/'

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
from matplotlib  import cm, colors
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

'denormalize in def-Input'
denor   = 0

'''start here'''
print('---   [',var_tit1,'] reading time:', T_id,' ---')
outputs = srd.rd_file(T_id, var_tit1, path_id,denor)
title1  = outputs[0]
x_title = outputs[1]
y_tit1  = outputs[2]
t_out   = outputs[3]
num_x   = outputs[4]
num_y   = outputs[5]
x_axis  = outputs[6]
y_axis1 = outputs[7]
data1   = outputs[8]

print('---         T= %.2e'%t_out,'  c/wp   ---')

print('Y(%i): %s,         mid y: (%4.2e --- %4.2e)'%(num_y,title1[:-1],\
        min(data1[int(num_y/2),:]),max(data1[int(num_y/2),:]) ) )
print('                    Peak x= %.2e'%( x_axis[np.argmax(data1[int(num_y/2),:])] ) )

outputs = srd.rd_file(T_id, var_tit2, path_id,denor)
title2  = outputs[0]
y_tit2  = outputs[2]
y_axis2 = outputs[7]
#data2   = np.sqrt(outputs[8]/3.0)
data2   = outputs[8]
print('Y(%i): %s,      mid y: (%4.2e --- %4.2e)'%(num_y,title2[:-1],\
        min(data2[int(num_y/2),:]),max(data2[int(num_y/2),:]) ) )
outputs = srd.rd_file(T_id, var_tit3, path_id,denor)
title3  = outputs[0]
y_tit3  = outputs[2]
y_axis3 = outputs[7]
data3   = outputs[8]
print('Y(%i): %s,         mid y: (%4.2e --- %4.2e)'%(num_y,title3[:-1],\
        min(data3[int(num_y/2),:]),max(data3[int(num_y/2),:]) ) )

outputs = srd.rd_file(T_id, var_tit4, path_id,denor)
title4  = outputs[0]
y_tit4  = outputs[2]
y_axis4 = outputs[7]
data4   = outputs[8]
print('Y(%i): %s,    mid y: (%4.2e --- %4.2e)'%(num_y,title4[:-1],\
        min(data4[int(num_y/2),:]),max(data4[int(num_y/2),:]) ) )

outputs = srd.rd_file(T_id, var_tit5, path_id,denor)
title5  = outputs[0]
x_tit5  = outputs[1]
y_tit5  = outputs[2]
x_axis5 = outputs[6]
y_axis5 = outputs[7]
data5   = outputs[8]
print('Y(%i): %s,       mid y: (%4.2e --- %4.2e)'%(num_y,title5[:-1],\
        min(data5[int(num_y/2),:]),max(data5[int(num_y/2),:]) ) )
#print('Min:', np.min(data5),'   (p,x)=',np.unravel_index(data5.argmin(), data5.shape),' in the map')
tmp = np.unravel_index(data5.argmin(), data5.shape)
print('Min: %.2e'%(np.min(data5)),'    [p,x]  = [%.2e, %.2e]'%(x_axis5[tmp[1]], y_axis5[tmp[0]]))

'''plots'''
fig = plt.figure(figsize=(12.,11))
#fig = plt.figure(figsize=(6.,5.5), dpi=80)
#plt.figure(1)

plt.subplot(511)
plt.plot(x_axis,data1[int(num_y/2),:])
plt.title('t= %.2e'%t_out)
plt.ylabel(title1)

plt.subplot(512)
plt.plot(x_axis,data2[int(num_y/2),:])
#plt.plot(x_axis,data2[int(num_y/2),:],label='Vt/c')
#plt.plot(x_axis,np.abs(data3[int(num_y/2),:]),label='|Vx/c|')
y2 = max(data2[int(num_y/2),:])
y1 = min(data2[int(num_y/2),:])
axes = plt.gca()
axes.set_ylim([1.5*y1-0.5*y2, 1.5*y2-0.5*y1])
#plt.title(title2)
plt.ylabel(title2)
plt.legend()

ax = plt.subplot(513)
ax.plot(x_axis,data3[int(num_y/2),:])
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
ax.yaxis.set_major_formatter(ticker.FuncFormatter(srd.fmt))
#plt.title(title3)
plt.ylabel(title3)

plt.subplot(514)
U = data4[0,:]
U[0] = 0
for i in range(1,num_x-1):
    #U[i] = (data4[0,i]+data4[0,i-1])*0.5*(x_axis[i]-x_axis[i-1])+U[i-1]
    U[i] = data4[0,i]*(x_axis[i]-x_axis[i-1])+U[i-1]
U = U*2/data3[0,125]**2

plt.plot(x_axis,data4[int(num_y/2),:])
#plt.plot(x_axis,U)
plt.ylabel(title4)

plt.subplot(515)

(x1,x2)=data5.shape
for i in range(0,x1):
    for j in range(0,x2):
        if data5[i,j] <0:
            data5[i,j] = data5[i,j] #* 1.0e-4
bar5h = np.max(data5)
bar5l = np.min(data5)
if bar5l <=0:
    #bar5l = -np.max((np.abs(bar5l),bar5h))
    levels     = np.linspace(bar5l,bar5h,300)
else:
    levels     = np.logspace(np.log10(bar5l),np.log10(bar5h),300)
cs = plt.contourf(y_axis5, x_axis5, data5.transpose(),cmap='jet', levels=levels)
#cs = plt.contourf(y_axis5, x_axis5, data5.transpose(),cmap='jet', norm=matplotlib.colors.SymLogNorm(linthresh=1.0,linscale=1.0,vmin=bar5l,vmax=bar5h),levels=levels)

#cs = plt.contourf(y_axis5, x_axis5, data5.transpose(),cmap='jet')
#cs = plt.contourf(data5)
cbar = plt.colorbar(cs,format=mpl.ticker.FuncFormatter(srd.fmt))

plt.xlabel(x_title)
plt.ylabel(x_tit5)

plt.show()
#fig.savefig('t1.eps',dpi=10)

import pickle
out_data = [title1, x_title, y_tit1, t_out, num_x, num_y, x_axis, y_axis1, data1, title2, y_tit2, data2, title3, y_tit3, y_axis3, data3, title4, y_tit4, y_axis4, data4, title5, x_tit5, x_axis5, y_axis5, data5, bar5l, bar5h]
pickle.dump(out_data,open("save_t2.p","wb"))
