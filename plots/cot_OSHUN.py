T_id    = '00000'
#var_tit = 'Ex'
#var_tit = 'Vx'
#var_tit = 'N'
#var_tit = 'Vsq'
var_tit = 'p1x1'
#var_tit = 'fsp'
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
'''
contourf on 2D plane at a given time
'''

denor   = 1
'denormalize in def-Input'

plog = 1
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
print('Min:', np.min(data),' at Np=',np.unravel_index(data.argmin(), data.shape),' in the map')

'''plot here'''
'''colorbar is using logspace'''
'''cm.jet is the default'''

fig = plt.figure()#figsize=(6.3,8))
plt.subplot(211)

if x_axis.size < 3 or y_axis.size <3:
    print('1D case')
    #need to work on this later
    #cs = plt.contourf(data,levels=levels,cmap=plt.cm.jet,norm=colors.LogNorm())
    cs = plt.contourf(data,cmap=plt.cm.jet)
else:
    #######
    #for POP.Tz.SecIII.A
    #cs = plt.contourf(x_axis, y_axis[int(num_y/2):], data[int(num_y/2):,:],levels=levels,cmap=plt.cm.jet,norm=colors.LogNorm())
    #######
    if plog == 1:
        #symmetric -/+ range
        barh = np.max(data)
        barl = np.min(data)
        if barl <=0:
            barl = -np.max((np.abs(barl),barh))
            levels     = np.linspace((np.min(data)),(np.max(data)),300)
            levels_bar = np.linspace((np.min(data)),(np.max(data)),2)
            levels     = np.linspace(0.0,(np.max(data)),300)
            levels_bar = np.linspace(0.0,(np.max(data)),2)
        else:
            levels     = np.logspace(np.log10(np.min(data)),np.log10(np.max(data)),300)
            levels_bar = np.logspace(np.log10(np.min(data)),np.log10(np.max(data)),2)
        ######
        #for POP.Tz.SecIV.B
        #cs = plt.contourf(x_axis[25:84], y_axis[int(num_y/2):], data[int(num_y/2):,25:84],cmap='jet', norm=colors.SymLogNorm(linthresh=1.0,linscale=-1.0,vmin=barl,vmax=barh),levels=levels)
        ######
        #cs = plt.contourf(x_axis, y_axis, data,cmap='jet', norm=colors.SymLogNorm(linthresh=1.0,linscale=1.0,vmin=barl,vmax=barh),levels=levels)
        ######
        #for p1x1
        cs = plt.contourf(y_axis, x_axis, data.transpose(),cmap='jet', norm=colors.SymLogNorm(linthresh=1.0,linscale=1.0,vmin=barl,vmax=barh),levels=levels)
        ######
        cbar = plt.colorbar(cs,format=mpl.ticker.FuncFormatter(srd.fmt))
        #cbar.set_ticks(levels_bar)
        #cbar.ax.set_ylabel(title)
    elif plog == 0:
        levels     = np.linspace((np.min(data)),(np.max(data)),50)
        levels_bar = np.linspace((np.min(data)),(np.max(data)),2)
        ######
        #for POP.Tz.SecIV.B
        #cs = plt.contourf(x_axis[25:84], y_axis[int(num_y/2):], data[int(num_y/2):,25:84],levels=levels,cmap=plt.cm.jet)
        ######
        cs = plt.contourf(x_axis, y_axis, data,levels=levels,cmap=plt.cm.jet)
        cbar = plt.colorbar(cs,format=mpl.ticker.FuncFormatter(srd.fmt))
        #cbar.set_ticks(levels_bar)
        #cbar.ax.set_ylabel(title)
    elif plog == 2:
        X,Y = np.mgrid[np.min(y_axis):np.max(y_axis):complex(0,num_y),np.min(x_axis):np.max(x_axis):complex(0,num_x)]
        cbar=plt.pcolormesh(X,Y,data,norm=colors.SymLogNorm(linthresh=0.0001,linscale=1.0,vmin=np.min(data),vmax=np.max(data)),cmap='jet')
        fig.colorbar(cbar)
    else:
        #for all
        cs = plt.contourf(x_axis, y_axis, data)
        cbar = plt.colorbar(cs,format=mpl.ticker.FuncFormatter(srd.fmt))
        #cbar.ax.set_ylabel(title)
    plt.xlabel(x_title)
    plt.ylabel(y_title)

#cbar = plt.colorbar(cs,format='%2.1e')
#cs2 = plt.contour(cs,levels=cs.levels[::50],hold='on')
#cbar.add_lines(cs2)

plt.title(title)
#plt.ticklabel_format(style='sci', scilimits=(0,0) )
plt.subplot(212)
plt.plot(x_axis,data[int(num_y/2),:])
#plt.plot(x_axis,data[0,:],'o')

plt.show()



