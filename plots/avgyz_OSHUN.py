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

T_max   = 300
var_tit1 = 'Vy'
var_tit2 = 'Vz'
path_id = '../OUTPUT/'
ppos    = 1

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
        #av1 = np.zeros(num_x)
        #av2 = np.zeros(num_x)
    elif i == 1:
        av1 =  outputs1[8][:,0]
        av2 =  outputs2[8][:,0]
    else:
        av1 +=  outputs1[8][:,0]
        av2 +=  outputs2[8][:,0]

av1 = av1/T_max
av2 = av2/T_max

fig = plt.figure()
plt.plot(av1,label='Vy')
plt.plot(av2,label='Vz')
plt.legend()
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
