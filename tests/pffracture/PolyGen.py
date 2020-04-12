#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 17:13:18 2020

@author: Yang Bai
@purpose: for the input file generation of polycrystalline simulation
""" 
import numpy as np

inputfilename='tensile50.i'
meshfile='rect50.msh'
grains=50
E1=210.0;nu1=0.3
E2=90.0;nu2=0.2
Gc=2.7e-3
L=4.0e-3
viscosity=1.0e-5

Alpha=65.0;dAlpha=25.0
Beta=75.0;dBeta=15.0
Gamma=90.0;dGamma=5.0

Theta1=np.random.normal(Alpha, dAlpha, grains)
Theta2=np.random.normal(Beta,dBeta,grains)
Theta3=np.random.normal(Gamma,dGamma,grains)

inp=open(inputfilename,'w')


# for mesh block
str='[mesh]\n'
str+='\ttype=gmsh\n'
str+='\tfile=%s\n'%(meshfile)
str+='[end]\n\n'
inp.write(str)

# for dof block
str='[dofs]\n'
str+='\tname=ux uy d\n'
str+='[end]\n\n'
inp.write(str)

# for projection block
str='[projection]\n'
str+='\tname=vonMises Hydro sigxx sigyy sigxy\n'
str+='[end]\n\n'
inp.write(str)

# for elmts
str='[elmts]\n'
inp.write(str)
for i in range(grains):
    str='\t[block%d]\n'%(i+1)
    str+='\t\ttype=miehefrac\n'
    str+='\t\tdofs=ux uy d\n'
    str+='\t\tmate=anisofrac%g\n'%(i+1)
    str+='\t\tdomain=%d\n'%(i+1)
    str+='\t[end]\n'
    inp.write(str)
str='[end]\n\n'
inp.write(str)

# for elmts
str='[mates]\n'
inp.write(str)
for i in range(grains):
    theta1=Theta1[i]
    theta2=Theta2[i]
    theta3=Theta3[i]
    str='\t[anisofrac%d]\n'%(i+1)
    str+='\t\ttype=anisopffrac\n'
    str+='\t\tparams=%g %g %g %g %g %g %g %g %g %g\n'%(E1,E2,nu1,nu2,Gc,L,viscosity,theta1,theta2,theta3)
    str+='\t[end]\n'
    inp.write(str)
str='[end]\n\n'
inp.write(str)


str='[bcs]\n'

str+='\t[fixUx]\n'
str+='\t\ttype=dirichlet\n'
str+='\t\tdof=ux\n'
str+='\t\tboundary=left\n'
str+='\t\tvalue=0.0\n'
str+='\t[end]\n'

str+='\t[fixUy]\n'
str+='\t\ttype=dirichlet\n'
str+='\t\tdof=uy\n'
str+='\t\tboundary=bottom\n'
str+='\t\tvalue=0.0\n'
str+='\t[end]\n'

str+='\t[loadUy]\n'
str+='\t\ttype=dirichlet\n'
str+='\t\tdof=uy\n'
str+='\t\tboundary=top\n'
str+='\t\tvalue=1.0*t\n'
str+='\t[end]\n'

str+='[end]\n\n'
inp.write(str)

str ='[timestepping]\n'
str+='\ttype=be\n'
str+='\tdt=1.0e-5\n'
str+='\tdtmax=4.0e-4\n'
str+='\tdtmin=5.0e-9\n'
str+='\tendtime=1.0e2\n'
str+='\topts=4\n'
str+='\tadaptive=true\n'
str+='[end]\n\n'
inp.write(str)

str='[nonlinearsolver]\n'
str+='\ttype=newtonls\n'
str+='\tr_abs_tol=5.0e-7\n'
str+='\tr_rel_tol=1.5e-8\n'
str+='\tmaxiters=25\n'
str+='[end]\n\n'
inp.write(str)


str='[job]\n'
str+='\ttype=transient\n'
str+='\tdebug=dep\n'
str+='\tprojection=true\n'
str+='[end]\n'
inp.write(str)




inp.close()