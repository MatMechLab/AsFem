import matplotlib.pyplot as plt

import numpy as np

tol=1.0e-2;nx=500
c=np.linspace(tol,1.0-tol,nx)



def F(x):
    return ((x-0.2)**2)*((x-0.8)**2)

def dF(x):
    return 4*(x**3-1.5*x**2+0.66*x-0.08)

def d2F(x):
    return 12*(x*x-x+0.22)


calpha=0.2
cbeta=0.8
cmid=0.5
calphas=(5-np.sqrt(3))/10
cbetas=(5+np.sqrt(3))/10

print('calphas=%g,cbetas=%g'%(calphas,cbetas))

f=F(c[:])
df=dF(c[:])
d2f=d2F(c[:])

fig, ax1 = plt.subplots()
ax1.plot(c,f,'k',label='F')
ax1.plot(calpha,F(calpha),'k*')
ax1.plot(cbeta,F(cbeta),'k*')
ax1.plot(cmid,F(cmid),'k*')
ax1.plot(calphas,F(calphas),'ro')
ax1.plot(cbetas,F(cbetas),'ro')
# for spinodal decomposition area
line1x=np.linspace(calphas,calphas,10)
line1y=np.linspace(-0.001,0.025,10)
ax1.plot(line1x,line1y,'r--')
# for spinodal decomposition area
line2x=np.linspace(cbetas,cbetas,10)
line2y=np.linspace(-0.001,0.025,10)
ax1.plot(line2x,line2y,'r--')
ax1.legend()
ax1.set_ylabel('Free energy',fontsize=15)
ax1.text(0.35,0.022,'spinodal',color='red',fontsize=13)
ax1.text(0.35,0.020,'decomposition',color='red',fontsize=13)

ax1.text(0.1,0.022,r'$\alpha$-phase',color='black',fontsize=13)
ax1.text(0.72,0.022,r'$\beta$-phase',color='black',fontsize=13)

ax1.text(0.18,0.002,r'$c_{\alpha}$',fontsize=16)
ax1.text(0.78,0.002,r'$c_{\beta}$',fontsize=16)

ax1.text(calphas+0.02,0.0025,r'$c_{\alpha}^{s}$',color='red',fontsize=16)
ax1.text(cbetas-0.07,0.0025,r'$c_{\beta}^{s}$',color='red',fontsize=16)


ax1.text(0.1,0.015,'stable',color='black',fontsize=15)
ax1.text(0.4,0.015,'unstable',color='red',fontsize=15)
ax1.text(0.7,0.015,'stable',color='black',fontsize=15)


ax2 = ax1.twinx()
ax2.plot(c,df,'b-',label='$\mu$')
ax2.plot(calpha,dF(calpha),'k*')
ax2.plot(cbeta,dF(cbeta),'k*')
ax2.plot(cmid,dF(cmid),'k*')

ax2.plot(calphas,dF(calphas),'r*')
ax2.plot(cbetas,dF(cbetas),'r*')

ax2.set_ylabel('Chemical potential',fontsize=15)
ax2.legend()
x=np.linspace(0.0,1.0,20)
y=np.linspace(0.0,0.0,20)
ax2.plot(x,y,'b--')

plt.savefig('free.jpg',dpi=500,bbox_inches='tight')
