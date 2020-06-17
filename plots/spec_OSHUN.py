T_id    = '00021'
var_tit = 'Ex'
path_id = '../OUTPUT/'


import numpy as np
import scipy as sci
import scipy.special as sp
from scipy import ndimage
import cmath
from math import e
#put my modules here
import sub_rd_OSHUN as srd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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
'function'
def spec_perp(data, dim=2):
    sz_x, sz_y = data.shape
    if dim == 1 :
        nl = np.max([sz_x,sz_y])
        if sz_x > sz_y :
            f_space = np.fft.fft(data[:,0])
        else:
            f_space = np.fft.fft(data[0,:])
        
        space_ps  = np.abs(f_space)
        space_ps *= space_ps
        nm = int(nl/2)
        fd = np.ndarray.copy(space_ps[0:nm])
        fd[1:nm] = fd[1:nm] + space_ps[:nm:-1]
        p  = np.ndarray.copy(fd[0:nm])

    elif dim == 2:
        nx = int(sz_x/2)
        ny = int(sz_y/2)

        f_space   = np.fft.fft2(data,s=[sz_x,sz_y])
        space_ps  = np.abs(f_space)
        space_ps *= space_ps

        fd = np.ndarray.copy(space_ps[0:nx,0:ny])
        fd[1:nx,0:ny] = fd[1:nx,0:ny] + space_ps[:nx:-1,0:ny]
        fd[0:nx,1:ny] = fd[0:nx,1:ny] + space_ps[0:nx,:ny:-1]
        fd[1:nx,1:ny] = fd[1:nx,1:ny] + space_ps[:nx:-1,:ny:-1] 

        fs = np.ndarray.copy(fd[0:nx, 0:ny])
        nm = np.min([nx,ny])
        kd = np.linspace(1,nm,nm)
        p  = np.zeros(nm)
        p[0] = (fs[0,1]+fs[1,0])*np.pi*0.25

        for i in range(1,nm):
            k   = kd[i]
            npk = 3*k+1
            dt  = 0.5*np.pi/(npk-1)
            x   = k*np.cos(np.arange(npk)*dt)
            y   = k*np.sin(np.arange(npk)*dt)
            coo = [y,x]
            pk  =ndimage.map_coordinates(fs,coo,order=3,mode='nearest')
            p[i]= 0.5*np.pi*k*np.sum(pk)/npk

    else:
        print('Dimension should be <= 2D!')
        sys.exit(0)

    return p

###################################################################'read data'
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
print('max point: %4.1e, min point: %4.1e'%(np.max(data),np.min(data)))
print('ABS max point: %4.1e, min point: %4.1e'%(np.max(np.abs(data)),np.min(np.abs(data))))

if x_axis.size < 3 or y_axis.size <3:
    print('1D case')
    spec = spec_perp(data,dim=1)
else:
    spec = spec_perp(data,dim=2)


'''plot here'''
fig = plt.figure(figsize=(16.8,6.5))
gs = gridspec.GridSpec(1,2,width_ratios=[1,1.2])
ax1= plt.subplot(gs[:,0])
plt.contourf(data)
ax1= plt.subplot(gs[:,1])
plt.plot(spec)
plt.xlabel('$k _perp$')
plt.ylabel('E($k _perp$)')
plt.title(title)
plt.show()



