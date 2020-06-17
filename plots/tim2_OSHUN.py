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
'''
test particle tracing on 2D plane.
'''

T_max   = 200
var_tit1 = 'Vx'
var_tit2 = 'Vy'
path_id = '../OUTPUT/'
ppos    = 1

data1 = np.zeros(T_max)
data2 = np.zeros(T_max)
tim   = np.zeros(T_max)
posx  = np.zeros(T_max+1)
posy  = np.zeros(T_max+1)
#initial postion
posx[0] = 300
posy[0] = 0
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
    if i == 0:
        num_x   = outputs1[4]
        num_y   = outputs1[5]
        x_axis  = outputs1[6]
        y_axis  = outputs1[7]
        idx = x_axis[1]-x_axis[0]
        idy = y_axis[1]-y_axis[0]
    #visual position in the domain
    pos_x = posx[i]
    if pos_x > x_axis[num_x-1]:
            pos_x = pos_x -int( (pos_x-x_axis[0])/(x_axis[num_x-1]-x_axis[0]))*(x_axis[num_x-1]-x_axis[0])
    elif pos_x < x_axis[0]:
            pos_x = pos_x -int( (pos_x-x_axis[num_x-1])/(x_axis[num_x-1]-x_axis[0]))*(x_axis[num_x-1]-x_axis[0])
            #pos_x = pos_x + x_axis[num_x-1] - x_axis[0]
    if pos_x < x_axis[0] or pos_x > x_axis[num_x-1] :
        print(i,'time, x too large!')
        sys.exit(0)

    pos_y = posy[i]
    if pos_y > y_axis[num_y-1]:
            pos_y = pos_y -int( (pos_y-y_axis[0])/(y_axis[num_y-1]-y_axis[0]))*(y_axis[num_y-1]-y_axis[0])
            #pos_y = pos_y - y_axis[num_y-1] + y_axis[0]
    elif pos_y < y_axis[0]:
            pos_y = pos_y -int( (pos_y-y_axis[num_y-1])/(y_axis[num_y-1]-y_axis[0]))*(y_axis[num_y-1]-y_axis[0])
            #pos_y = pos_y + y_axis[num_y-1] - y_axis[0]
    if pos_y < y_axis[0] or pos_y > y_axis[num_y-1] :
        print(i,'time, y too large!')
        sys.exit(0)
    #position in the 2D box
    dx = (pos_x - x_axis[0])/idx
    dy = (pos_y - y_axis[0])/idy
    x1 = int(dx)
    dx = dx-x1
    x2 = x1+1
    if x2 > num_x-1:
        x2 = x2-num_x
    y1 = int(dy)
    y2 = y1+1
    if y2 > num_y-1:
        y2 = y2-num_y
    dy = dy-y1
    #interpotation
    data1[i] = outputs1[8][y1,x1]*(1-dx)*(1-dy) + outputs1[8][y1,x2]*dx*(1-dy)\
                + outputs1[8][y2,x1]*(1-dx)*dy + outputs1[8][y2,x2]*dx*dy
    #print(i,': ',data1[i], outputs1[8][x1,y1],outputs1[8][x2,y1],outputs1[8][x1,y2],outputs1[8][x2,y2],dx,dy)
    #print(i,': ',data1[i], x1,y1,dx,dy)
    #print(i,': ',data1[i] )

    data2[i] = outputs2[8][y1,x1]*(1-dx)*(1-dy) + outputs2[8][y1,x2]*dx*(1-dy)\
                + outputs2[8][y2,x1]*(1-dx)*dy + outputs2[8][y2,x2]*dx*dy
    tim[i]  = outputs1[3]
    if i == 1:
        dt = tim[1]-tim[0]
        posx[1] = posx[0]+dt*data1[0]
        posy[1] = posy[0]+dt*data2[0]
    if i>0:
        posx[i+1] = posx[i]+dt*data1[i]
        posy[i+1] = posy[i]+dt*data2[i]
        #posy[i+1] = posy[i]+dt*(data2[i]+v_drift)

vmin = np.min([np.min(data1),np.min(data2)])
vmax = np.max([np.max(data1),np.max(data2)])
#plt.plot(tim1,data1,label=var_tit1,linestyle='-',marker='+')
plt.plot(tim,data1,label=var_tit1)
plt.plot(tim,data2,label=var_tit2)
#plt.ylim([vmin,vmax])
plt.legend()
plt.show()

if ppos == 1:
    fig = plt.figure(2)
    plt.plot(posx,posy)
    plt.plot(posx[0],posy[0],marker='o',color='r')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()
