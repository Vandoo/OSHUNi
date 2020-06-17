T_max   = 100
#var_tit1 = 'Vx'
var_tit2 = 'Ex'

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib  import cm, colors
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
#put my modules here
import sub_rd_OSHUN as srd

#for my mac...
import warnings
warnings.simplefilter(action="ignore",category= FutureWarning)
warnings.simplefilter(action="ignore",category= UserWarning)
#debug
import sys
'''
for grad B drift, Bx(y)
'''

path_id = '../OUTPUT/'
ppos    = 1
tim     = np.zeros(T_max)
ymt     = np.zeros(T_max)


for i in range(0,T_max):
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
    #read files
    #outputs1 = srd.rd_file(T_id, var_tit1, path_id)
    outputs2 = srd.rd_file(T_id, var_tit2, path_id)
    if i == 0:
        num_y   = outputs2[5]
        y_axis  = outputs2[7]
        dy      = y_axis[1]-y_axis[0]
        av2     = np.zeros((T_max,num_y),dtype=np.float64)
    tim[i]  = outputs2[3]
    av2[i,:]= outputs2[8][:,0]
    ymt[i]  = y_axis[np.argmax(av2[i,:])]
    if i == 1:
        if ymt[i] > ymt[i-1]:
            dr = 1
        else:
            dr = -1
    if i>1:
        if dr == 1 and ymt[i]<ymt[i-1]:
            ymt[i] = np.max(y_axis)+dy + ymt[i] - np.min(y_axis)
        if dr == -1 and ymt[i]>ymt[i-1]:
            ymt[i] = np.min(y_axis)-dy + ymt[i] - np.max(y_axis)
    #if i>0:
     #   av2[i,:] =av2[i-1,:]+outputs2[8][0,:]*(tim[i]-tim[i-1])


#y_axis = np.linspace(-np.pi,np.pi,num_y)
#vm = np.max(np.abs(av2))
#levels = np.linspace(-vm,vm,50)
m,b = np.polyfit(tim,ymt,1)
print('V_',var_tit2,'=',m)

fig  = plt.figure()
#cs   = plt.contourf(y_axis,tim,av2,cmap=cm.jet,levels=levels)
#cs   = plt.contourf(y_axis,tim,av2,cmap=cm.jet)
#cbar = plt.colorbar(cs,format=mpl.ticker.FuncFormatter(srd.fmt))
plt.plot(tim,ymt)
plt.plot(tim,b+m*tim)
plt.xlabel('t')
plt.ylabel('y')
plt.title(var_tit2)

plt.show()


'''
b = np.linspace(0.001,np.pi*0.999,50)
c = np.linspace(np.pi*1.001,np.pi*1.999,50)
e = np.zeros(100)
e[0:50] = b
e[50:100]= c
a = np.cos(e)
f = a/(1-a**2)
plt.plot(e,f)
plt.show()
'''
