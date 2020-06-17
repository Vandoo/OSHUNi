from numpy import *

un  = 'cgs'
'''ion or electron'''
spe = 'i'
#other parameters
Zp  = 1
#np  = 1.0e16#4.0e20 #cm^-3
#Tev = 1.0e0
np  = 1.0e13#4.0e20 #cm^-3
Tev = 1.0e0
#Tev = 5.0e-1
B   = 5.0e-01


print('----------------------------------------------------------')
if un == 'cgs':
    me  = 9.1094e-28
    mi  = 1.6726e-24
    ce  = 4.8032e-10 #e
    cs  = 2.9979e10  #light speed
elif un == 'si':
    me  = 9.1094e-31
    mi  = 1.6726e-27
    ce  = 1.6022e-19  #e
    cs  = 2.9979e8    #light speed
else:
    print('Unkown units')

if spe == 'e':
    print('Electron version:')
    lnGee = 23.5-0.5*log(np)+1.25*log(Tev)-sqrt(1.e-5+0.0625*(Tev-2)**2)
    lnG = 24-0.5*log(np)+log(Tev)
    mm = me
    Ze = 1
elif spe == 'i':
    print('Ion version:')
    lnG = 23.-0.5*log(2.0*np)+1.5*log(Tev)-3*log(Zp)
    mm = mi
else:
    print('Unkown species')

if lnG <2:
    lnG = 2

rp  = ce**2/mm/cs**2
print('r_p*: %e *Z'%(rp))

wp  = sqrt(4*pi*ce**2/mm)
print('wp* : %e *Z*sqrt(np)'%(wp))

#Export
print('-----')

beta = 1.0e7
epsilon0 = 8.8542e-12
alpha = 1.0e2
Esw  = 1.0/sqrt(4*pi*beta*epsilon0/alpha**3)
print('esu/cm -> V/m: %e'%(Esw))

E_tvpm = mm*cs**2/ce * Esw *1.0e-12
print('E_TV/m*: %e /Z/cwp'%(E_tvpm))

B_mg = mm*cs**2/ce *1.0e-6
print('B_MG*  : %e /Z/cwp'%(B_mg))

J_sa = cs*ce
print('J_sa*  : %e *np*Z'%(J_sa))
#For calculate the simulation parameters:
print('-------------------------------')
print('Simulation parameters based on:')
print('Density: %.1ecm^-3, Z=%i, T=%.1eeV, B=%.1eG'%(np,Zp,Tev,B))
print('-------------------------------')

print('-----')
if spe == 'e':
    w0  = wp*Ze*sqrt(np)
elif spe == 'i':
    w0  = wp*Zp*sqrt(np)
else:
    print('Unkown species')

t0  = 1.0/w0
x0  = cs/w0
cwp = x0
vc  = sqrt(Tev*1.6e-12/mm)/cs
print('x0(cwp)   : %e cm  = %e um'%(x0,x0*1.e4))
print('t0(1/wp)  : %e sec = %e fs'%(t0,t0*1.e15))
print('v/c       : %e '%(vc))
print()


kb  = 1.3807e-16
Tk  = Tev/8.621738e-5
t_zz = 0.75*sqrt(mm/pi)/np/Zp/(ce)**4/lnG*(Tk*kb)**1.5
print('ln(Gamma) : %.2e '%(lnG))
if spe == 'e':
    print('Tau_ei    : %.2e sec = %.2e fs, Tau_e = Tau_ei*Zp'%(t_zz,t_zz*1.e15))
elif spe == 'i':
    print('Tau_ii    : %.2e sec = %.2e fs'%(t_zz,t_zz*1.e15))
lambda_mfp = t_zz*vc*cs
print('lambda_mfp: %.2e cm  = %.2e um'%(lambda_mfp,lambda_mfp*1.e4))
print()


lambda_d = sqrt(kb*Tk/4/pi/np/ce**2)
N_d      = 4/3*pi*np*lambda_d**3
print('lambda_d  : %.2e cm  = %.2e um'%(lambda_d,lambda_d*1.e4))
print('N_d       : %.2e '%(N_d))
wpTau = 9*N_d/sqrt(2/pi)/Zp/lnG
print('wp*Tau_ei : %.2e '%(wpTau))
print()


w_g = B*ce/mm/cs
T_g = 2*pi/w_g
r_g = mm*vc/B/ce*cs**2
va  = B/sqrt(4*pi*np*mi)
print('If B = %.1e G:'%(B))
print('T_gyro(2pi/w_g) = %.2e sec = %.2e fs'%(T_g,T_g*1.e15))
print('R_gyro          = %.2e cm  = %.2e um'%(r_g,r_g*1.e4))
print('VA              = %.2e cm/s'%(va))

print('----------------------------------------------------------')
const1 = (2.0*me)**1.5*Zp*ce/3.0/sqrt(mi)
const2 = 16*sqrt(pi)*ce**2/3/me/(cs**2/me/1.6e-12)**1.5
print('(2*me)^1.5 *Ze/3/(mi*ni)^0.5 = %.4e /sqrt(ni) or %.4e'%(const1, const2))
print('----------------------------------------------------------')
