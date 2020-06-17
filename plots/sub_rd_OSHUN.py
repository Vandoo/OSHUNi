import numpy as np
#for debug
#from imp import reload
import sys

########################################################################
#file path
def file_name(var,denor=0):
    '''Vsq, T, Gam, P, N(=x1x2),
       Vx/y/z, Px, Qx/Kx, Qy/Ky, Ex/y/z, Bx/y/z, Jx/y/z
       p1x1, p1x1x2, fsp'''
    if var == 'Vsq' or var =='V':
        title_id = 'MOM/Vsq/Vsq'
    elif var == 'T':
        title_id = 'MOM/T_eV/T'
    elif var == 'P':
        title_id = 'MOM/P_Mbar/P'
    elif var == 'Gam':
        title_id = 'MOM/Gam/Gam'
    elif var == 'N' or var == 'x1x2':
        title_id = 'MOM/N/x1x2'
    elif var == 'Vx':
        title_id = 'MOM/Vx/Vx'
    elif var == 'Vy':
        title_id = 'MOM/Vy/Vy'
    elif var == 'Vz':
        title_id = 'MOM/Vz/Vz'
    elif var == 'Px':
        title_id = 'MOM/Px/GVx'
    elif var == 'Py':
        title_id = 'MOM/Py/GVy'
    elif var == 'Pz':
        title_id = 'MOM/Pz/GVz'
    elif var == 'Qx':
        title_id = 'MOM/Qx/Qx'
    elif var == 'Kx':
        title_id = 'MOM/Qx/Kx'
    elif var == 'Qy':
        title_id = 'MOM/Qx/Qy'
    elif var == 'Ky':
        title_id = 'MOM/Qy/Ky'
    elif var == 'Ex':
        if(denor == 0):
            title_id = 'FLD/EX/Ex'
        else:
            title_id = 'FLD/EX/Ex_d'
    elif var == 'Ey':
        if(denor == 0):
            title_id = 'FLD/EY/Ey'
        else:
            title_id = 'FLD/EY/Ey_d'
    elif var == 'Ez':
        if(denor == 0):
            title_id = 'FLD/EZ/Ez'
        else:
            title_id = 'FLD/EZ/Ez_d'
    elif var == 'Bx':
        if(denor == 0):
            title_id = 'FLD/BX/Bx'
        else:
            title_id = 'FLD/BX/Bx_d'
    elif var == 'By':
        if(denor == 0):
            title_id = 'FLD/BY/By'
        else:
            title_id = 'FLD/BY/By_d'
    elif var == 'Bz':
        if(denor == 0):
            title_id = 'FLD/BZ/Bz'
        else:
            title_id = 'FLD/BZ/Bz_d'
    elif var == 'Jx':
        title_id = 'FLD/JX/Jx'
    elif var == 'Jy':
        title_id = 'FLD/JY/Jy'
    elif var == 'Jz':
        title_id = 'FLD/JZ/Jz'
    elif var == 'p1x1':
        title_id = 'DISTR/P1X1/p1x1'
    elif var == 'fsp':
        if(denor == 0):
            title_id = 'DISTR/Fs/Fl0'
        if(denor == 1):
            title_id = 'DISTR/Fs/Fl1'
        if(denor == 2):
            title_id = 'DISTR/Fs/Fl2'
        if(denor == 3):
            title_id = 'DISTR/Fs/Fs0'
        if(denor == 4):
            title_id = 'DISTR/Fs/Fs1'
        if(denor == 5):
            title_id = 'DISTR/Fs/Fs2'
    elif var == 'p1x1x2':
        title_id = 'DISTR/P1X1X2/p1x1x2'
    elif var == 'pmulti':
        title_id = 'DISTR/Pmulti/pm'
    else:
        print('      Not Plotting! ', title_id)
        sys.exit(0)

    #print('          Plotting:',var)
    return title_id
########################################################################
#read file
########################################################################
def rd_file(T_id,var_tit,path_id,denor=0):
	'''
	a[0,1,2]: title,xtit,tit]
	a[3]: t_out
	a[4,5]:num x,y
	a[6,7]:axis x,y
	a[8]: data
	'''
	title_id = file_name(var_tit,denor)
	dname = path_id+title_id+"_"+T_id+".txt"
	try:
		dfile = open(dname,'r')
	except IOError:
		if (var_tit[0] == 'E' or var_tit[0] == 'B') and denor == 0:
			#title_id = 'FLD/EX/Ex_d'
			#print('          denormalized E/B')
			title_id = file_name(var_tit,1)
			dname = path_id+title_id+"_"+T_id+".txt"
			dfile = open(dname,'r')
		else:
			print('Could not open the file!', dname)
			sys.exit(0)
	except:
		sys.exit(0)

	title   = dfile.readline()
	x_title = dfile.readline()
	y_title = dfile.readline()
	t_out   = np.float64(dfile.readline())
	num_x   = int(dfile.readline())
	num_y   = int(dfile.readline())
	
	x_axis  = np.zeros(num_x)
	y_axis  = np.zeros(num_y)
	for i in range(0, num_x):
		x_axis[i] = np.float64(dfile.readline())
	for i in range(0, num_y):
		y_axis[i] = np.float64(dfile.readline())

	'''Notice: not Fortran style array for contourf'''
	data = np.zeros((num_y,num_x),dtype=np.float64)
	for j in range(0, num_y):
		for i in range(0, num_x):
			#print(i,j)
			data[j,i] = np.float64(dfile.readline())
	
	return [title,x_title,y_title, t_out, num_x, num_y,
			x_axis, y_axis,data]

########################################################################
#for T_id
########################################################################
def T_index(i):
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

    return T_id
########################################################################
########################################################################
#for colorbar ticklabel
########################################################################
def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${}\times10^{{{}}}$'.format(a, b)
def fmt3(x, pos):
    a, b = '{:.3e}'.format(x).split('e')
    b = int(b)
    return r'${}\times10^{{{}}}$'.format(a, b)
