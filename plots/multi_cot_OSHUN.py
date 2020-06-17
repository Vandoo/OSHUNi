T_min    = 21
T_max    = 2
var_tit1 = 'N'
var_tit2 = 'p1x1'
var_tit3 = 'Vsq'
var_tit4 = 'Vx'
var_tit5 = 'Ex'
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
for movie
'''

'denormalize in def-Input'
denor   = 0

'''short compare'''
def cptwo(a,b,d):
    if d == 1:
        if a>=b:
                return a
        else:
                return b
    elif d == 2:
        if a>=b:
                return b
        else:
                return a
    else:
        print('Not numbers compared here!')
        sys.exit(0)

'''start here, find data range'''
print()
print('Now checking contour scales...')
i = T_min
if i == T_min:
    if i<10:
        T_id = '0000'+str(i)
    elif i<100 and i>9:
        T_id = '000'+str(i)
    elif i<1000 and i>99:
        T_id = '00'+str(i)
    elif i<10000 and i>999:
        T_id = '0'+str(i)
    else:
        print('T > 9999?')
        sys.exit(0)

    #read file #1
    outputs = srd.rd_file(T_id, var_tit1, path_id, denor)
    barh = np.max(outputs[8])
    barl = np.min(outputs[8])
	
    if i == T_min:
        bar1h = barh
        bar1l = barl
    else:
        bar1h = cptwo(bar1h,barh,1)
        bar1l = cptwo(bar1l,barl,2)
		
    #file #2
    outputs = srd.rd_file(T_id, var_tit2, path_id, denor)
    barh = np.max(outputs[8])
    barl = np.min(outputs[8])
    if i == T_min:
        bar2h = barh
        bar2l = barl
    else:
        bar2h = cptwo(bar2h,barh,1)
        bar2l = cptwo(bar2l,barl,2)
    #file #3
    outputs = srd.rd_file(T_id, var_tit3, path_id, denor)
    barh = np.max(outputs[8])
    barl = np.min(outputs[8])
    if i == T_min:
        bar3h = barh
        bar3l = barl
    else:
        bar3h = cptwo(bar3h,barh,1)
        bar3l = cptwo(bar3l,barl,2)
    #file #4
    outputs = srd.rd_file(T_id, var_tit4, path_id, denor)
    barh = np.max(outputs[8])
    barl = np.min(outputs[8])
    if i == T_min:
        bar4h = barh
        bar4l = barl
    else:
        bar4h = cptwo(bar4h,barh,1)
        bar4l = cptwo(bar4l,barl,2)
    #file #5
    outputs = srd.rd_file(T_id, var_tit5, path_id, denor)
    #outputs1 = srd.rd_file(T_id, var_tit3, path_id, denor)
    #outputs[8] = outputs[8]/sqrt(outputs1[8]/3)
    barh = np.max(outputs[8])
    barl = np.min(outputs[8])
    if i == T_min:
        bar5h = barh
        bar5l = barl
    else:
        bar5h = cptwo(bar5h,barh,1)
        bar5l = cptwo(bar5l,barl,2)
	

'''real plots'''
print()
print('New loop:')
i = T_min
if i == T_min:
    if i<10:
        T_id = '0000'+str(i)
    elif i<100 and i>9:
        T_id = '000'+str(i)
    elif i<1000 and i>99:
        T_id = '00'+str(i)
    elif i<10000 and i>999:
        T_id = '0'+str(i)
    else:
        print('T > 9999?')
        sys.exit(0)
    print('Now plotting...',T_id)
    #read file #1
    outputs = srd.rd_file(T_id, var_tit1, path_id, denor)
    if i == T_min:
        x_title = outputs[1]
        y_title = outputs[2]
        num_x   = outputs[4]
        num_y   = outputs[5]
        x_axis  = outputs[6]
        y_axis  = outputs[7]
    title1  = outputs[0]
    t_out   = outputs[3]
    data1   = outputs[8]
    #file #2
    outputs = srd.rd_file(T_id, var_tit2, path_id, denor)
    title2  = outputs[0]
    data2   = outputs[8]
    #for p1x1
    if i == T_min:
        xp_axis  = outputs[6]
        yp_axis  = outputs[7]
    #file #3
    outputs = srd.rd_file(T_id, var_tit3, path_id, denor)
    title3  = outputs[0]
    data3   = outputs[8]
    #file #4
    outputs = srd.rd_file(T_id, var_tit4, path_id, denor)
    title4  = outputs[0]
    data4   = outputs[8]
    #file #5
    outputs = srd.rd_file(T_id, var_tit5, path_id, denor)
    #outputs1 = srd.rd_file(T_id, var_tit3, path_id, denor)
    #outputs[8] = outputs[8]/sqrt(outputs1[8]/3)
    title5  = outputs[0]
    data5   = outputs[8]
    #plots
    fig = plt.figure(figsize=(12.,18.))
    
    plt.subplot(511)
    if bar1l <=0:
        bar1l = -np.max((np.abs(bar1l),bar1h))
        levels     = np.linspace(bar1l,bar1h,300)
        levels_bar = np.linspace(bar1l,bar1h,2)
    else:
        levels     = np.logspace(np.log10(bar1l),np.log10(bar1h),300)
        levels_bar = np.logspace(np.log10(bar1l),np.log10(bar1h),2)
    cs = plt.contourf(x_axis, y_axis, data1,cmap='jet', norm=matplotlib.colors.SymLogNorm(linthresh=1.0,linscale=1.0,vmin=bar1l,vmax=bar1h),levels=levels)
    cbar = plt.colorbar(cs,format=mpl.ticker.FuncFormatter(srd.fmt))
    plt.title('T = %.2e'%(t_out))
    plt.text(x_axis[int(num_x*0.86)],y_axis[int(num_y*0.74)],title1,fontsize=12)
    #plt.xlabel(x_title)
    plt.ylabel(y_title)

    plt.subplot(512)
    if bar2l <=0:
        bar2l = -np.max((np.abs(bar2l),bar2h))
        levels     = np.linspace(bar2l,bar2h,300)
    else:
        levels     = np.logspace(np.log10(bar2l),np.log10(bar2h),300)
    #cs = plt.contourf(x_axis, y_axis, data2,cmap='jet', norm=matplotlib.colors.SymLogNorm(linthresh=1.0,linscale=1.0,vmin=bar2l,vmax=bar2h),levels=levels)
    #cbar = plt.colorbar(cs,format=mpl.ticker.FuncFormatter(srd.fmt))
    #plt.text(x_axis[int(num_x*0.86)],y_axis[int(num_y*0.74)],title2,fontsize=12)

    #For p1x1
    levels     = np.linspace(0,bar2h,300)
    #cs = plt.contourf(yp_axis, xp_axis, data2.transpose(),cmap='jet', norm=matplotlib.colors.SymLogNorm(linthresh=1.0,linscale=1.0,vmin=bar2l,vmax=bar2h),levels=levels)
    cs = plt.contourf(yp_axis, xp_axis, data2.transpose(),cmap='jet', levels=levels)
    #cs = plt.contourf(yp_axis, xp_axis, data2.transpose(),cmap='jet')
    cbar = plt.colorbar(cs,format=mpl.ticker.FuncFormatter(srd.fmt))

    #plt.xlabel(x_title)
    plt.ylabel(y_title)

    plt.subplot(513)
    if bar3l <=0:
        bar3l = -np.max((np.abs(bar3l),bar3h))
        levels     = np.linspace(bar3l,bar3h,300)
    else:
        levels     = np.logspace(np.log10(bar3l),np.log10(bar3h),300)
    cs = plt.contourf(x_axis, y_axis, data3,cmap='jet', norm=matplotlib.colors.SymLogNorm(linthresh=1.0,linscale=1.0,vmin=bar3l,vmax=bar3h),levels=levels)
    cbar = plt.colorbar(cs,format=mpl.ticker.FuncFormatter(srd.fmt))
    plt.text(x_axis[int(num_x*0.86)],y_axis[int(num_y*0.74)],title3,fontsize=12)
    #plt.xlabel(x_title)
    plt.ylabel(y_title)

    plt.subplot(514)
    if bar4l <=0:
        bar4l = -np.max((np.abs(bar4l),bar4h))
        levels     = np.linspace(bar4l,bar4h,300)
    else:
        levels     = np.logspace(np.log10(bar4l),np.log10(bar4h),300)
    cs = plt.contourf(x_axis, y_axis, data4,cmap='jet', norm=matplotlib.colors.SymLogNorm(linthresh=1.0,linscale=1.0,vmin=bar4l,vmax=bar4h),levels=levels)
    cbar = plt.colorbar(cs,format=mpl.ticker.FuncFormatter(srd.fmt))
    plt.text(x_axis[int(num_x*0.86)],y_axis[int(num_y*0.74)],title4,fontsize=12)
    #plt.xlabel(x_title)
    plt.ylabel(y_title)

    plt.subplot(515)
    if bar5l <=0:
        bar5l = -np.max((np.abs(bar5l),bar5h))
        levels     = np.linspace(bar5l,bar5h,300)
    else:
        levels     = np.logspace(np.log10(bar5l),np.log10(bar5h),300)
    cs = plt.contourf(x_axis, y_axis, data5,cmap='jet', norm=matplotlib.colors.SymLogNorm(linthresh=1.0,linscale=1.0,vmin=bar5l,vmax=bar5h),levels=levels)
    cbar = plt.colorbar(cs,format=mpl.ticker.FuncFormatter(srd.fmt))
    plt.text(x_axis[int(num_x*0.86)],y_axis[int(num_y*0.74)],title5,fontsize=12)
    plt.xlabel(x_title)
    plt.ylabel(y_title)

    #savefig('Anime/t_%04d.eps'%(i))
    #plt.close(fig)
plt.show()
'''
for step in range(0,10):
    plt.plot(x_axis,data1[int(num_y/2),:])
    plt.title('test %d'%(step))
    savefig('anime/t_%03d.eps'%(step))
    '''
