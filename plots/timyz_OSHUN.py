import numpy as np
from math import e
import matplotlib.pyplot as plt
#put my modules here
import sub_rd_OSHUN as srd

#for my mac...
import warnings
warnings.simplefilter(action="ignore",category= FutureWarning)
warnings.simplefilter(action="ignore",category= UserWarning)
#debug
import sys

T_max   = 20
var_tit1 = 'Vy'
var_tit2 = 'Vz'
path_id = '../OUTPUT/'
ppos    = 1

data1 = np.zeros(T_max)
data2 = np.zeros(T_max)
tim  = np.zeros(T_max)
posx = np.zeros(T_max+1)
posy = np.zeros(T_max+1)
#data postion
pos_x = 0
pos_y = 0
v_drift = 0

for i in range(0,T_max):
    if i<10:
        T_id = '0000'+str(i)
    elif i<100 and i>9:
        T_id = '000'+str(i)
    elif i<1000 and i>99:
        T_id = '00'+str(i)
    else:
        print('T > 9999?')
        sys.exit(0)
    #read files
    outputs1 = srd.rd_file(T_id, var_tit1, path_id)
    outputs2 = srd.rd_file(T_id, var_tit2, path_id)
    data1[i] = outputs1[8][pos_x,pos_y]
    data2[i] = outputs2[8][pos_x,pos_y]
    tim[i]   = outputs1[3]
    if i == 1:
        dt = tim[1]-tim[0]
        posx[1] = posx[0]+dt*data1[0]
        posy[1] = posy[0]+dt*data2[0]
    if i>0:
        posx[i+1] = posx[i]+dt*data1[i]
        posy[i+1] = posy[i]+dt*(data2[i]+v_drift)
'''
'''
vmin = np.min([np.min(data1),np.min(data2)])
vmax = np.max([np.max(data1),np.max(data2)])
#plt.plot(tim1,data1,label=var_tit1,linestyle='-',marker='+')
plt.plot(tim,data1,label=var_tit1)
plt.plot(tim,data2,label=var_tit2)
plt.ylim([vmin,vmax])
plt.legend()
plt.show()

if ppos == 1:
    fig = plt.figure(2)
    plt.plot(posx,posy)
    plt.plot(posx[0],posy[0],marker='o',color='r')
    plt.xlabel('y')
    plt.ylabel('z')
    plt.show()
