import numpy as np
import matplotlib.pyplot as plt


def Y(alpha,Yieldstres,Kinf,K0,delta):
	return Yieldstres+(Kinf-K0)*(1.0-np.exp(-delta*alpha))


Sigma0=0.8
Kinf=3.5;K0=1.0;delta=35.0
x=np.linspace(0.0,0.4,200)
y=Y(x,Sigma0,Kinf,K0,delta)

plt.figure(1)
plt.plot(x,y)
plt.xlabel('effective plastic strain',fontsize=13)
plt.ylabel('Yield stress')

plt.show()