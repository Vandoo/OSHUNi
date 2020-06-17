T_id    = '00001'
vp_tit  = 'V'
path_id  = '../OUTPUT/'
print(vp_tit,'is plotted.')


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
vector visualization on a plane at a given time
'''


ptype   = 'sc'
'vector fields on 2D plane. V, J, E, B?'
'qc=quive color, sc=stream color'
'qs=quive single,sl=stream line width'


print('---   Reading time:', T_id,' ---')
'''start here'''
if vp_tit =='V' or vp_tit == 'v':
    var_tit1 = 'Vx'
    #var_tit2 = 'Vy'
    var_tit2 = 'Vz'
elif vp_tit =='E' or vp_tit == 'e':
    var_tit1 = 'Ex'
    var_tit2 = 'Ey'
elif vp_tit =='J' or vp_tit == 'j':
    var_tit1 = 'Jx'
    var_tit2 = 'Jy'
elif vp_tit =='B' or vp_tit == 'b':
    var_tit1 = 'Bx'
    var_tit2 = 'By'
else:
    print('Dont know the field!')

#x dir
outputs = srd.rd_file(T_id, var_tit1, path_id)
title1  = outputs[0]
x_title = outputs[1]
y_title = outputs[2]
t_out   = outputs[3]
num_x   = outputs[4]
num_y   = outputs[5]
x_axis  = outputs[6]
y_axis  = outputs[7]
data1   = outputs[8]
print('X(%i): [%6.2f,%6.2f], mid x: (%4.1e --- %4.1e)'%(\
        num_x,x_axis[0],x_axis[num_x-1],\
        min(data1[:,int(num_x/2)]),max(data1[:,int(num_x/2)]) ) )
print('Y(%i): [%6.2f,%6.2f], mid y: (%4.1e --- %4.1e)'%(\
        num_y,y_axis[0],y_axis[num_y-1],\
        min(data1[int(num_y/2),:]),max(data1[int(num_y/2),:]) ) )

#y dir
outputs = srd.rd_file(T_id, var_tit2, path_id)
title2  = outputs[0]
data2   = outputs[8]
print('X(%i): [%6.2f,%6.2f], mid x: (%4.1e --- %4.1e)'%(\
        num_x,x_axis[0],x_axis[num_x-1],\
        min(data2[:,int(num_x/2)]),max(data2[:,int(num_x/2)]) ) )
print('Y(%i): [%6.2f,%6.2f], mid y: (%4.1e --- %4.1e)'%(\
        num_y,y_axis[0],y_axis[num_y-1],\
        min(data2[int(num_y/2),:]),max(data2[int(num_y/2),:]) ) )

'''plot here'''
amp = np.sqrt(data1**2+data2**2)
d1n = data1/amp
d2n = data2/amp

fig = plt.figure()#figsize=(6.3,8))
if(ptype == 'qc'):
    plt.quiver(x_axis, y_axis, d1n, d2n,
               data1,
               cmap=cm.seismic,
               headlength=7)
    plt.colorbar()
    plt.title('Quive Plot: '+title1[0:2]+'/'+title2+' Dynamic Colours')
elif(ptype == 'qs'):
    plt.quiver(x_axis, y_axis, d1n, d2n,
            color='Teal',
            headlength=7)
    plt.title('Quive Plot: '+title1[0:2]+'/'+title2+' Single Colour')
elif(ptype == 'sc'):
    plt.streamplot(x_axis, y_axis,data1, data2,
            color=amp,
            cmap=cm.cool,
            linewidth=1,
            arrowstyle='->',
            arrowsize=0.5)
    plt.colorbar()
    plt.title('Stream Plot: '+title1[0:2]+'/'+title2+' Dynamic Colours')
elif(ptype == 'sl'):
    lw = amp/amp.max()
    plt.streamplot(x_axis, y_axis,data1, data2,
            density=[0.5,1],
            color='DarkRed',
            linewidth=lw)
    plt.title('Stream Plot: '+title1[0:2]+'/'+title2+' Dynamic ling Width')

plt.xlabel(x_title)
plt.ylabel(y_title)
plt.show()
