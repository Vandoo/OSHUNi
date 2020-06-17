var_tit = 'Bx'
T_list = ('00001','00201','00401','00601','00801','01001','01201','01401')
#T_list = ('00001','00101','00201','00301','00401','00501','00601','00701')
#T_list = ('00025','00050','00075','00100','00125','00150','00200')
#T_list = ('00000','00015','00030','00045','00060','00075','00090','00105','00120')
#T_list = ('00000','00012','00024','00036','00048','00060','00080','00100')
#T_list = ('00000','00012','00024','00036','00048','00060')
#T_list = ('00000','00006','00012','00018','00024','00030')
#T_list = ('00000','00004','00008','00012','00016','00018')
#T_list = ('00000','00002','00004','00006','00008')
#T_list = ('00000','00001','00002','00003','00004')
#T_list = ('00000','00001')
path_id = '../OUTPUT/'


print('[',var_tit,'] is plotted.')

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
co-plot 1 variable at different time in one plot
'''

'denormalize in def-Input'
denor   = 0
'1D: direction x:4, y:5'
'plot the middle line'
direc_id = 5

t_len = len(T_list)

for T_id in T_list:
    outputs = srd.rd_file(T_id, var_tit, path_id,denor)
    if(T_id == T_list[0]):
        var = np.zeros(shape=(len(T_list),outputs[direc_id]))
        tim = np.zeros(len(T_list))
        num = 0
        if(direc_id == 4):
            axis = outputs[6]
        if(direc_id == 5):
            axis = outputs[7]

    if(direc_id == 4):
        var[num,:]= outputs[8][int(outputs[5]/2),:]
    if(direc_id == 5):
        var[num,:]= outputs[8][:,int(outputs[4]/2)]
    tim[num] = outputs[3]
    num +=1

'plot'
for i in range(0,num):
    #plt.plot(axis,var[i,:],label=T_list[i])
    plt.plot(axis,var[i,:],label=tim[i])#,marker='o')
    print('|',var_tit,'| ---',T_list[i],'--- :',np.mean(np.abs(var[i,:])))
    print('Peak at [',np.argmax(var[i,:]),']:', var[i,np.argmax(var[i,:])])

#plt.xlim(0.0,30000)
plt.title(var_tit)
plt.legend()
plt.show()
