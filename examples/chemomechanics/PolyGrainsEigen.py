'''
This script is used for the input file generation of 
polycrystalline structe. It will ouput the [mates] and [elmts]
block for AsFem.
@author: Yang Bai
@Date  : 2021.11.07
@Usage : PolyGrains.py -grains 5
'''
import numpy as np
import sys

D=1.5e1;Omega=0.08
E=210.0;nu=0.25

Gc=2.7e-3;L=0.01;viscosity=1.0e-6
UseHist=0

nGrains=50

i=0
for name in sys.argv:
    if '-grains' in name:
        if i+1<len(sys.argv):
            nGrains=int(sys.argv[i+1])
            if nGrains<1:
                nGrains=1
    i+=1
###################################
### for random euler angles 
###################################
E0=160+np.random.normal(50,25,nGrains)
Omega=0.001+np.random.normal(0.04,0.0095,nGrains)
###################################
filename='polygrains.txt'
inp=open(filename,'w')

##########################
### for [elmts] block 
##########################
inp.write('[elmts]\n')
for i in range(nGrains):
    str ='  [myelmt%d]\n'%(i+1)
    str+='    type=diffusionfracture\n'
    str+='    dofs=c d ux uy\n'
    str+='    mate=mymate%d\n'%(i+1)
    str+='    domain=%d\n'%(i+1)
    str+='  [end]\n'
    inp.write(str)
inp.write('[end]\n\n\n')

#################################
### for [mates] block
#################################
inp.write('[mates]\n')
str=''
for i in range(nGrains):
    #E=E0[i]
    str+='  [mymate%d]\n'%(i+1)
    str+='    type=diffusionfracturemate\n'
    str+='    params=%11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n'%(D,Omega[i],Gc,L,viscosity,E,nu,UseHist)
    str+='  [end]\n'
str+='[end]\n'
inp.write(str)
