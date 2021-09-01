#!/usr/bin/env python3
'''
This script plot the double well free energy and its derivatives
@Author: Yang Bai
@Date  : 2021.08.31
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def F(x,ca,cb,factor):
    return factor*(x-ca)**2*(x-cb)**2

def dF(x,ca,cb,factor):
    return factor*2*(x-ca)*(x-cb)*(2*x-ca-cb)

def d2F(x,ca,cb,factor):
    return 2*(ca**2+4*ca*cb-6*ca*x+cb**2-6*cb*x+6*x**2)

ca=0.2;cb=0.8;factor=1.0

c=np.linspace(0.0,1.0,100)
y0=np.linspace(0.0,0.0,100)
f=F(c,ca,cb,factor)
df=dF(c,ca,cb,factor)

##################################
### find the roots of f and mu
##################################
root1=fsolve(dF,[0.0,1.0],args=(ca,cb,factor))
root2=fsolve(d2F,[0.0,1.0],args=(ca,cb,factor))
print('roots for dF=0 : ',root1)
print('roots for d2F=0: ',root2)

fig, ax1 = plt.subplots()
ax1.plot(c,f,'k',label=r'$\psi$',linewidth=2.0)
ax1.set_ylabel('free energy [-]',fontsize=15)
ax1.set_xlabel('concentration [-]',fontsize=15)

### for spinodal points
x1=root2[0];y1=F(x1,ca,cb,factor)
txt='%g'%(x1)
ax1.text(x1+0.025,y1-0.001,txt,fontsize=14)
x2=root2[-1];y2=F(x2,ca,cb,factor)
txt='%g'%(x2)
ax1.text(x2+0.025,y2-0.001,txt,fontsize=14)

ax1.plot(x1,y1,'r*',markersize=15)
ax1.plot(x2,y2,'r*',markersize=15)

ax1.tick_params(axis='both', which='major', labelsize=13)

ax1.legend(loc=1,fontsize=13)

ax2 = ax1.twinx()
ax2.plot(c,df,'b',label=r'$\mu$',linewidth=2.0)
ax2.plot(c,y0,'k--')
ax2.legend(loc=2,fontsize=13)

x=np.linspace(root2[0],root2[0],50)
y=np.linspace(-0.3,0.3,50)
ax2.plot(x,y,'r--',linewidth=2.0)

x=np.linspace(root2[-1],root2[-1],50)
ax2.plot(x,y,'r--',linewidth=2.0)

ax2.set_ylabel('chemical potential [-]',fontsize=15)

ax2.set_xlim([0.0,1.0])

plt.show()
