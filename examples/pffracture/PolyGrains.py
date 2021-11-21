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

Lambda1=121.15;mu1=80.77
Lambda2=Lambda1/2;mu2=mu1/2
Lambda3=Lambda1/2;mu3=mu1/2
K1=Lambda1+2*mu1/3.0
K2=Lambda2+2*mu2/3.0;K3=K2

# elastic constants
C11=K1+4*mu1/3;C12=K1-2*mu2/3;C13=K1-2*mu3/3
C22=K2+4*mu2/3;C23=K2-2*mu3/3
C33=K3+4*mu3/3
C44=mu3;C55=mu2;C66=mu1

Gc=2.7e-3;L=0.012;viscosity=1.0e-6

nGrains=4
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
Theta1=np.random.rand(nGrains)*180.0
Theta2=np.random.rand(nGrains)*90.0
Theta3=np.random.rand(nGrains)*60.0
###################################
filename='polygrains.txt'
inp=open(filename,'w')

##########################
### for [elmts] block 
##########################
inp.write('[elmts]\n')
for i in range(nGrains):
    str ='  [myelmt%d]\n'%(i+1)
    str+='    type=miehefrac\n'
    str+='    dofs=d ux uy\n'
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
    theta1=Theta1[i];theta2=0.0;theta3=0.0
    str+='  [mymate%d]\n'%(i+1)
    str+='    type=stressdecompositionmate\n'
    str+='    params=%11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e '%(C11,C12,C13,C22,C23,C33,C44,C55,C66)
    str+='%11.4e %11.4e %11.4e '%(Gc,L,viscosity)
    str+='%11.4e %11.4e %11.4e\n'%(theta1,theta2,theta3)
    str+='  [end]\n'
str+='[end]\n'
inp.write(str)
