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

LambdaVec=121.15;muVec=80.77




Gc=2.7e-3;L=0.012;viscosity=1.0e-6

nGrains=250
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
#LambdaVec=np.random.rand(nGrains)*121.15
muVec=np.random.rand(nGrains)*80.77

LambdaVec=np.random.normal(121.15, 30, nGrains)
muVec    =np.random.normal(80.77, 25, nGrains)
###################################
filename='polygrainsE.txt'
inp=open(filename,'w')

##########################
### for [elmts] block 
##########################
inp.write('[elmts]\n')
for i in range(nGrains):
    str ='  [myelmt%d]\n'%(i+1)
    str+='    type=miehefrac\n'
    str+='    dofs=d ux uy uz\n'
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
    k=LambdaVec[i];mu=muVec[i]
    str+='  [mymate%d]\n'%(i+1)
    str+='    type=miehefracmate\n'
    str+='    params=%11.4e %11.4e %11.4e %11.4e %11.4e\n'%(k,mu,Gc,L,viscosity)
    str+='  [end]\n'
str+='[end]\n'
inp.write(str)
