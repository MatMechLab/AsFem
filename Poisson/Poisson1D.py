#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  6 14:03:25 2018

QQ group: 797998860 
@email: walkandthinker@gmail.com
@author: FlyFox
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

from Mesh import *
from GaussPoint import *
from ShapeFuns import *

def Welcome():
    print('********************************************')
    print('*** Welcome to use this simple script!   ***')
    print('*** Simple fem code for laplace equation ***')
    print('*** Author: FlyFox                       ***')
    print('*** Email: walkandthinker@gmail.com      ***')
    print('*** QQ group: 797998860                  ***')
    print('********************************************')


F=1.0;A=0.5 # for equation looks like: A*div(grad(phi))=F
PhiL=0.5    # for left side boundary condition phi(x=xmin)=PhiL
#######################################################
### For analytical solution of 1D linear laplace equation
#######################################################
def AnalyticalSolution(NodeCoords,xmin,xmax):
    nNodes=np.size(NodeCoords)
    U0=np.zeros(nNodes)
    B=-(1.0*F/A)*xmax
    C=PhiL-(0.5*F/A)*xmin*xmin+F*xmin*xmax/A

    for i in range(nNodes):
        U0[i]=(0.5*F/A)*NodeCoords[i]**2+B*NodeCoords[i]+C
    
    return U0

#######################################################
### You can modify the user element for your own model
#######################################################
def elmt(elCoords,u):
    # Here F and A is just two constants for different solution of laplace equation
    nnodes=np.size(elCoords)
    Lint=nnodes
    gs=Int1D(Lint)

    k=np.zeros((nnodes,nnodes)) # for local k matrix
    rhs=np.zeros(nnodes)        # for local rhs

    for gp in range(Lint):
        xi=gs[gp,1]
        shp,xsj=Shp1D(elCoords,xi)
        JxW=xsj*gs[gp,0] # |J|*weight

        phi=0.0;gradphi=0.0
        for i in range(nnodes):
            phi    +=u[i]*shp[i,0]
            gradphi+=u[i]*shp[i,1]  


        for I in range(nnodes):
            rhs[I]+=-F*shp[I,0]*JxW-A*gradphi*shp[I,1]*JxW
            for J in range(nnodes):
                k[I,J]+=A*shp[J,1]*shp[I,1]*JxW
    
    return k,rhs

def FormKR(mesh,U):
    nElmts=mesh.nElmts
    nNodesPerElmt=mesh.nNodesPerElmts
    nDofs=mesh.nNodes
    AMATRX=np.zeros((nDofs,nDofs)) # [K]{u}=F-->here AMATRX is K
    RHS=np.zeros(nDofs)

    for e in range(nElmts):
        elConn=mesh.Conn[e,:]
        print(elConn)
        elCoords=mesh.NodeCoords[elConn-1]
        elU=U[elConn-1]

        k,rhs=elmt(elCoords,elU)

        # Assemble k,rhs to system matrix
        for i in range(nNodesPerElmt):
            RHS[elConn[i]-1]+=rhs[i]
            for j in range(nNodesPerElmt):
                AMATRX[elConn[i]-1,elConn[j]-1]+=k[i,j]
    

    # Apply dirichlet boundary condition for the left side node
    AMATRX[0,0]=1.0e12
    RHS[0]=0.0


    print('*** System equation: Ax=F generated!     ***')
    print('***--------------------------------------***')
    return AMATRX,RHS


if __name__=="__main__":
    Welcome()
    xmin=0.0;xmax=1.0 # define solution domain
    ne=10;P=3         # define element number and order

    mesh=Mesh1D(ne,xmin,xmax,P)
    mesh.CreateMesh()


    
    nDofs=mesh.nNodes
    U=np.zeros(nDofs)

    print('*** Node number=%5d'%(mesh.nNodes))
    print('*** Node number per element=%5d'%(mesh.nNodesPerElmts))
    print('*** Element number=%5d'%(mesh.nElmts))
    print('*** Dofs number=%5d'%(nDofs))
    
    # start newton-raphson iteration
    iters=0;MaxIters=50
    atol=1.0e-12;rtol=1.0e-9 # absolute error and relative error
    IsConvergent=False
    while iters<MaxIters and (not IsConvergent):
        U[0]=PhiL # for left side dirichlet boundary condition
        AMATRX,RHS=FormKR(mesh,U)
        dU=np.linalg.solve(AMATRX,RHS)
        U[:]+=dU[:] # update, U=U+deltU

        if iters==0:
            R0_norm=np.linalg.norm(RHS)
            dU0_norm=np.linalg.norm(dU)
        
        R_norm=np.linalg.norm(RHS)
        dU_norm=np.linalg.norm(dU)

        if (R_norm<rtol*R0_norm or R_norm<atol) or (dU_norm<rtol*dU0_norm or dU_norm<atol):
            IsConvergent=True
            print('Iteration=%2d, |R|=%.5e, |dU=|%.5e'%(iters,R_norm,dU_norm))
            break
        
        iters+=1


    
    Uanalytical=AnalyticalSolution(mesh.NodeCoords,xmin,xmax)

    error=np.linalg.norm(U-Uanalytical)
    print('|U-U_analytical|=%.5e'%(error))

    plt.plot(mesh.NodeCoords,U,'r*',label='FEM')
    plt.plot(mesh.NodeCoords,Uanalytical,label='Analytical')
    plt.legend()
    plt.savefig('compare.png',dpi=800)
    plt.show()