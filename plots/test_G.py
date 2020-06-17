
import numpy as np
import matplotlib.pyplot as plt

#test G function
x = np.linspace(0.0001,0.0021,10)
a2 = 2.0* 0.004421**2
f = x**2 * np.exp(-x**2/a2)
dx = x[1]-x[0]

fd = 2*x*np.exp(-x**2/a2) -2*x**3/a2*np.exp(-x**2/a2)
fdn = np.gradient(f,dx)

lf = x
lf = lf/lf[0]*fd[0]

#plt.plot(x,fd)
#plt.plot(x,fdn)
#plt.plot(x,lf)


f3 = 2*x**3/a2*np.exp(-x**2/a2)
f3n = x**3
f3n = f3n/f3n[0]*f3[0]

#plt.plot(f3)
#plt.plot(x,f3n)

x = np.linspace(1.0e-4,601.0e-4,300)
dx = x[1]-x[0]
f1 = x * np.exp(-x**2/(4.421e-3*4)**2)

### f2
f2i = np.gradient(f1/x,dx)*x

f2t = np.gradient(f1,dx)- f1/x

f2 = np.gradient(f1,dx)
#f2[1] = (f1[1]*0.25 - f1[2])/2/dx
f2 = f2 -f1/x
f2[0] = 0
#f2[1] = (f1[1]*0.25 - f1[2])/2/dx

### f3
f3i = np.gradient(f2i/x**2,dx)*x**2

f3t = np.gradient(f2t,dx) - 2/x*f2t
#f3t = np.gradient(f2t,dx)/(1+(dx/x)**2/3)- 2/x*f2t

#f3 = np.gradient(f2,dx)
f3 = np.gradient(f2,dx)
f3 = f3 - 2/x*f2
f3[0] = 0
f3[1] = 0

#f3[0]= -(f2[0]*0.25 - f2[1])/dx/2
#f3[1]= -(f2[1]*0.25 - f2[2])/dx/2

### f4
#f33 = x**3 * np.exp(-x**2/a2)
f4i = np.gradient(f3i/x**3,dx)*x**3

#f4t = -2/a2*x*f33
f4t = np.gradient(f3t,dx)- 3/x*f3
f4 = np.gradient(f3,dx)/(1+(dx/x)**2/3)- 3/x*f3

### test H function ###
'''
f4 = x**4 * np.exp(-x**2/a2)*1.0e-9
f4[0:3]=0.0
f3i = (2*4+1 - 2*x**2/a2)*x**3*np.exp(-x**2/a2)*1.0e-9
f3 = np.gradient(f4*x**5,dx)/x**5
f3t = 5/x*f4+np.gradient(f4,dx)
f3t[2] = f3t[4]*(x[2]/x[4])**3
f3t[3] = f3t[4]*(x[3]/x[4])**3
'''

### test p^l v.s. p^l*exp(-p^2/B) ###

pt = 0.004421
fe0 = np.exp(-x**2/2/pt**2) / (np.sqrt(2.0*np.pi)*2*np.pi*pt**3)
fl0 = (fe0[0] - fe0[1]*x[0]**2/x[1]**2)/(1-x[0]**2/x[1]**2) +(fe0[1]-fe0[0])*(x/x[1])**2
#plt.plot(fe0,label='exp')
#plt.plot(fl0,label='2nd')


### f0->f4 ###
B = 2*pt**2
f1 = np.gradient(fe0,dx)
f1[0] = 2*x[0]*(fe0[1]-fe0[0])/(x[1]**2-x[0]**2)
invB = -(np.log(fe0[0])-np.log(fe0[1]))/(x[1]**2-x[0]**2)
#f1[0] = -2*x[0]*invB*fe0[0]
print((-2*x[0]/B) * fe0[0],f1[0],-2*x[0]*invB*fe0[0])
#f1 = (-2*x/B) * fe0

f2 = np.gradient(f1/x,dx)*x
z2 = np.gradient(f1,dx)-f1/x
f2[0:1] = 0
z2[0:1] = 0

f3 = np.gradient(f2/x**2,dx)*x**2
z3 = np.gradient(z2,dx)-z2/x*2.0
f3[0:2] = 0
z3[0:2] = 0
#f3[2] = -f2[3]*(1-(3/7)**2)/2/dx
f4 = np.gradient(f3/x**3,dx)*x**3
z4 = np.gradient(z3,dx)-z3/x*3.0
f4[0:3] = 0
z4[0:3] = 0
#f4[3] = f3[4]*(1-(5/9)**3)/2/dx
f5 = np.gradient(f4/x**4,dx)*x**4
f5[0:4] = 0

#f4i = (4/B**2 - 8*x/B**3 + 4/B + 8*x**2/B**2 - 3*8*x**2/B**3 + 16*x**4/B**4)*fe0
f2i = (-2*x/B)**2 * fe0
f3i = (-2*x/B)**3 * fe0
f4i = (-2*x/B)**4 * fe0
f5i = (-2*x/B)**5 * fe0

plt.plot(f4i[:],label='ideal')
plt.plot(z4[:],label='Tz')
plt.plot(f4[:],label='num')

plt.legend(loc=4,fontsize='small')

plt.show(block=False)

