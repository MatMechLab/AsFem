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
from numba import jit # If you don't have numba, comment out it!

from Mesh import *
from GaussPoint import *
from ShapeFuns import *

def Welcome():
    print('********************************************')
    print('*** Welcome to use this simple script!   ***')
    print('*** Simple fem code for 2d CH model      ***')
    print('*** Author: FlyFox                       ***')
    print('*** Email: walkandthinker@gmail.com      ***')
    print('*** QQ group: 797998860                  ***')
    print('********************************************')

############################################################
### Governing equation for CH(mixed form):               ###
### dc/dt=div(M*grad(mu))                                ###
###    mu=df/dc-kappa*lap(c)                             ###
############################################################

############################################################
### define the free energy of your own model
############################################################
@jit
def FreeEnergyAndItsDerivative(c,choice):
    if choice==1:
        f=100*c*c*(1-c)*(1-c)
        dfdc=200*c*(c -1)*(2*c-1) # chemical potential
        d2fdc2=1200*c**2-1200*c+200
    elif choice==2:
        f=(c**4-2*c**2)/4.0
        dfdc=c**3-c
        d2fdc2=3*c**2-1.0
    
    return f,dfdc,d2fdc2
    # write your own free energy here
    #else:

############################################################
### For projection
############################################################
@jit
def Projection(nNodes,IX,xsj,shp,val,Proj):
    for i in range(nNodes):
        k=IX[i]
        xs=shp[i,0]*xsj
        Proj[k,0]+=xs
        for j in range(9):
            Proj[k,1+j]+=val[j]*xs

############################################################
### user element for cahn hillard equation
############################################################
@jit
def elmt(IX,elCoords,elU,Proj,ctan):
    nNodes=np.size(elCoords,0)
    nDofs=2*nNodes

    ngp=2
    if nNodes>=8:
        ngp=3
    
    gs=Int2D(ngp)
    Lint=ngp*ngp

    K=np.zeros((nDofs,nDofs))
    RHS=np.zeros(nDofs)
    val=np.zeros(9)

    kappa=2.0e-2 # for the thickness of interface
    M=1.0        # for the mobility term


    for gpInd in range(Lint):
        xi =gs[gpInd,1]
        eta=gs[gpInd,2]
        shp,xsj=Shp2D(elCoords,xi,eta)
        JxW=gs[gpInd,0]*xsj

        conc=0.0;cdot=0.0;gradc=np.zeros(3)
        mu=0.0;gradmu=np.zeros(3)
        for i in range(nNodes):
            conc    +=shp[i,0]*elU[2*i][0]
            cdot    +=shp[i,0]*elU[2*i][1]
            gradc[1]+=shp[i,1]*elU[2*i][0]
            gradc[2]+=shp[i,2]*elU[2*i][0]

            mu       +=shp[i,0]*elU[2*i+1][0]
            gradmu[1]+=shp[i,1]*elU[2*i+1][0]
            gradmu[2]+=shp[i,2]*elU[2*i+1][0]

        f,dfdc,d2fdc2=FreeEnergyAndItsDerivative(conc,1)
        
        
        #############################
        ### For projection value
        val[1-1]=f
        val[2-1]=dfdc
        val[3-1]=d2fdc2
        # Do projection from gauss point to nodal point
        Projection(nNodes,IX,xsj,shp,val,Proj)


        for iInd in range(nNodes):
            # For residual
            # For R_c
            RHS[2*iInd  ]+=cdot*shp[iInd,0]*JxW\
                          +M*(gradmu[1]*shp[iInd,1]+gradmu[2]*shp[iInd,2])*JxW
            # For R_mu
            RHS[2*iInd+1]+=mu*shp[iInd,0]*JxW-dfdc*shp[iInd,0]*JxW\
                          -kappa*(gradc[1]*shp[iInd,1]+gradc[2]*shp[iInd,2])*JxW
            # For jacobian              
            for jInd in range(nNodes):
                # Kc,cdot
                K[2*iInd  ,2*jInd  ]+=-shp[jInd,0]*shp[iInd,0]*JxW*ctan[1]
                # Kc,mu
                K[2*iInd  ,2*jInd+1]+=-M*(shp[jInd,1]*shp[iInd,1]+shp[jInd,2]*shp[iInd,2])*JxW*ctan[0]

                # Kmu,c
                K[2*iInd+1,2*jInd  ]+=d2fdc2*shp[jInd,0]*shp[iInd,0]*JxW*ctan[0]\
                                     +kappa*(shp[jInd,1]*shp[iInd,1]+shp[jInd,2]*shp[iInd,2])*JxW*ctan[0]
                # Kmu,mu
                K[2*iInd+1,2*jInd+1]+=-shp[jInd,0]*shp[iInd,0]*JxW*ctan[0]
        
    return K,RHS

#########################################################
@jit
def FormKR(mesh,U,V,AMATRIX,RHS,Proj,ctan):
    nElmts=mesh.nElmts
    nNodesPerElmt=mesh.nNodesPerElmts
    nDofs=mesh.nNodes*2
    AMATRIX[:,:]=0.0 # [K]{u}=F-->here AMATRX is K
    RHS[:]=0.0
    Proj[:,:]=0.0

    for e in range(nElmts):
        elConn=mesh.Conn[e,:]-1
        elCoords=mesh.NodeCoords[elConn]
        elU=np.zeros((2*nNodesPerElmt,2))
        for i in range(nNodesPerElmt):
            elU[2*i  ][0]=U[2*elConn[i]  ]
            elU[2*i+1][0]=U[2*elConn[i]+1]
            elU[2*i  ][1]=V[2*elConn[i]  ]
            elU[2*i+1][1]=V[2*elConn[i]+1]

        k,rhs=elmt(elConn,elCoords,elU,Proj,ctan)
        

        # Assemble k,rhs to system matrix
        for i in range(nNodesPerElmt):
            iInd=elConn[i]
            RHS[2*iInd  ]+=rhs[2*i  ]
            RHS[2*iInd+1]+=rhs[2*i+1]
            for j in range(nNodesPerElmt):
                jInd=elConn[j]
                AMATRIX[2*iInd  ,2*jInd  ]+=k[2*i  ,2*j  ]
                AMATRIX[2*iInd  ,2*jInd+1]+=k[2*i  ,2*j+1]
                AMATRIX[2*iInd+1,2*jInd  ]+=k[2*i+1,2*j  ]
                AMATRIX[2*iInd+1,2*jInd+1]+=k[2*i+1,2*j+1]
    
    # For projection value
    for i in range(mesh.nNodes):
        for j in range(9):
            Proj[i,1+j]/=Proj[i,0]
#######################################################
def ApplyDispBC(sidename,dofname,mesh,U,value):
    if sidename=='left':
        # For left edge
        BCConn=mesh.LeftBCConn
    elif sidename=='right':
        BCConn=mesh.RightBCConn
    elif sidename=='bottom':
        BCConn=mesh.BottomBCConn
    elif sidename=='top':
        BCConn=mesh.TopBCConn
    else:
        sys.exit('Side name=%s is invalid!!!'%(sidename))
    
    if dofname=='ux':
        component=1
    elif dofname=='uy':
        component=2
    else:
        sys.exit('dof name=%s in invalid!!!'%(dofname))
    
    ne=np.size(BCConn,0)
    nNodesPerBCElmt=np.size(BCConn,1)
    for e in range(ne):
        elConn=BCConn[e,:]-1
        for i in range(nNodesPerBCElmt):
            iInd=2*elConn[i]+component-1
            U[iInd]=value
########################################################
### Apply constrain condition
########################################################
def ApplyConstrainBC(sidename,dofname,mesh,AMATRIX,RHS):
    if sidename=='left':
        # For left edge
        BCConn=mesh.LeftBCConn
    elif sidename=='right':
        BCConn=mesh.RightBCConn
    elif sidename=='bottom':
        BCConn=mesh.BottomBCConn
    elif sidename=='top':
        BCConn=mesh.TopBCConn
    else:
        sys.exit('Side name=%s is invalid!!!'%(sidename))
    
    if dofname=='ux':
        component=1
    elif dofname=='uy':
        component=2
    else:
        sys.exit('dof name=%s in invalid!!!'%(dofname))
    
    ne=np.size(BCConn,0)
    nNodesPerBCElmt=np.size(BCConn,1)
    penalty=1.0e16
    for e in range(ne):
        elConn=BCConn[e,:]-1
        for i in range(nNodesPerBCElmt):
            iInd=2*elConn[i]+component-1
            AMATRIX[iInd,iInd]=penalty
            RHS[iInd]=0.0
#########################################################
@jit
def RandomIC(val,noise,mesh,U):
    np.random.seed()
    valmin=val-noise
    valmax=val+noise
    for i in range(mesh.nNodes):
        val=valmin+(valmax-valmin)*np.random.rand()
        U[2*i  ]=val
        U[2*i+1]=0.0
##########################################################
def PlotDisp(mesh,U,Flag):
    x=mesh.NodeCoords[:,0]
    y=mesh.NodeCoords[:,1]
    c =U[0::2]
    mu=U[1::2]

    plt.figure(1)
    plt.title('$c$',fontsize=16)
    plt.tricontourf(x,y,c,60)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.colorbar()
    plt.clim(np.min(c),np.max(c))
    plt.savefig('c.png',dpi=800,bbox_inches='tight')

    plt.figure(2)
    plt.title('$\mu$',fontsize=16)
    plt.tricontourf(x,y,mu,60)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.colorbar()
    plt.clim(np.min(mu),np.max(mu))
    plt.savefig('mu.png',dpi=800,bbox_inches='tight')

    if Flag==True:
        plt.show()
    

if __name__=='__main__':
    Welcome()
    nx=60;ny=60
    W=1.0;H=1.0
    mesh=Mesh2D(nx,ny,0.0,W,0.0,H,'quad4')
    mesh.CreateMesh()
    mesh.SplitBCMesh()
    
    nDofs=2*mesh.nNodes
    U=np.zeros(nDofs)
    U0=np.zeros(nDofs)
    V=np.zeros(nDofs)

    print('********************************************')
    print('*** number of node   = %6d            ***'%(mesh.nNodes))
    print('*** number of element= %6d            ***'%(mesh.nElmts))
    print('*** number of dofs   = %6d            ***'%(nDofs))
    print('*** nodes per element= %6d            ***'%(mesh.nNodesPerElmts))
    print('********************************************')
    print('*** Simulation start...                  ***')

    # Apply initial condition
    RandomIC(0.63,0.02,mesh,U0)

    #PlotDisp(mesh,U0,True)

    AMATRIX=np.zeros((nDofs,nDofs))
    RHS=np.zeros(nDofs)
    Proj=np.zeros((mesh.nNodes,10))

    nSteps=10

    dt=5.0e-6 # For phase field with random distribution, you can't give a large dt!!!
    ctan=np.zeros(2)
    ctan[0]=1.0;ctan[1]=1.0/dt
    atol=1.0e-12;rtol=1.0e-9 # absolute error and relative error
    MaxIters=50
    for step in range(nSteps):
        V[:]=0.0;U[:]=U0[:]
        iters=0;IsConvergent=False
        while iters<MaxIters and (not IsConvergent):
            FormKR(mesh,U,V,AMATRIX,RHS,Proj,ctan)

            dU=np.linalg.solve(AMATRIX,RHS)
            U[:]+=dU[:]
            V[:]=(U[:]-U0[:])*ctan[1]

            R_norm=np.linalg.norm(RHS)
            dU_norm=np.linalg.norm(dU)
            if iters==0:
                R0_norm=R_norm
                dU0_norm=dU_norm

            if R_norm<rtol*R0_norm or R_norm<atol:
                IsConvergent=True
                break
            iters+=1

        print('Step=%5d===> iters=%2d,|R0|=%12.5e,|R|=%12.5e,|dU0|=%12.5e,|dU|=%12.5e'%(step+1,iters,R0_norm,R_norm,dU0_norm,dU_norm))
        if IsConvergent:
            U0[:]=U[:]

    
    PlotDisp(mesh,U,True)




