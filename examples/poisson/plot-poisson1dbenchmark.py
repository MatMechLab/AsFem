'''
This script plot the benchmark results for 1d poisson problem
Author: Yang Bai
Date  : 2022.08.26
'''
import numbers
import numpy as np
import matplotlib.pyplot as plt
import re

def splitnum(str):
    numbers=[]
    substr=str.split('"')
    for i in substr:
        if i.isdigit():
            numbers.append(int(i))
    return numbers

def exactsolution(x):
    a=2.0;c0=a/12.0-0.5;c1=0.0
    return 0.5*x*x-a*x**4/12.0+c0*x+c1

def exactsolutionderiv(x):
    a=2.0;c0=a/12.0-0.5;c1=0.0
    return x-a*x**3/3+c0

filename='poisson-1d-benchmark.vtu'
inp=open(filename,'r')
line=inp.readline()
nodesnum=0;elmtnum=0
x=np.zeros(0)
phi=np.zeros(0);proj_phi=np.zeros(0)
gradphi=np.zeros(0)
while '</VTKFile>' not in line:
    if '<Piece NumberOfPoints=' in line:
        nodesnum,elmtnum=splitnum(line)
        x=np.zeros(nodesnum)
        phi=np.zeros(nodesnum)
        proj_phi=np.zeros(nodesnum)
        gradphi=np.zeros(nodesnum)
    elif '<DataArray type="Float64" Name="nodes"' in line:
        for i in range(nodesnum):
            line=inp.readline()
            linevalue=line.split()
            x[i]=float(linevalue[0])
    elif '<DataArray type="Float64" Name="phi"' in line:
        for i in range(nodesnum):
            line=inp.readline()
            phi[i]=float(line)
    elif '<DataArray type="Float64" Name="exactsolution"' in line:
        for i in range(nodesnum):
            line=inp.readline()
            proj_phi[i]=float(line)
    elif '<DataArray type="Float64" Name="gradu"' in line:
        for i in range(nodesnum):
            line=inp.readline()
            linevalue=line.split()
            gradphi[i]=float(linevalue[0])

    line=inp.readline()
inp.close()

sum=0.0
for i in range(nodesnum):
    sum+=(phi[i]-exactsolution(x[i]))**2
sum=np.sqrt(sum)/nodesnum
print('average error is: %14.6e'%(sum))

x0=np.linspace(0.0,1.0,5000)
phi_exact=exactsolution(x0)
phi_exact_deriv=exactsolutionderiv(x0)

plt.figure(1)
plt.plot(x,phi,'k*',markevery=2,label='AsFem solution')
plt.plot(x,proj_phi,'r*',markevery=2,label='AsFem projected solution')
plt.plot(x0,phi_exact,label='Exact solution')
plt.ylabel('solution',fontsize=15)
plt.xlabel('x',fontsize=15)
plt.xlim([0.0,1.0])
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.legend(fontsize=13)
plt.savefig('poisson1d-benchmark-phi.jpg',dpi=200,bbox_inches='tight')


plt.figure(2)
plt.plot(x,gradphi,'k*',markevery=2,label='AsFem solution derivative')
plt.plot(x0,phi_exact_deriv,label='Exact solution derivative')
plt.xlim([0.0,1.0])
plt.ylabel('solution derivative',fontsize=15)
plt.xlabel('x',fontsize=15)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.legend(fontsize=13)
plt.savefig('poisson1d-benchmark-gradphi.jpg',dpi=200,bbox_inches='tight')

plt.show()
