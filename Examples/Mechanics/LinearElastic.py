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
    print('*** Simple fem code for 2d linear elastic***')
    print('*** Author: FlyFox                       ***')
    print('*** Email: walkandthinker@gmail.com      ***')
    print('*** QQ group: 797998860                  ***')
    print('********************************************')

############################################################
### Projection 
############################################################
def Projection(nNodes,IX,xsj,shp,val,Proj):
    for i in range(nNodes):
        k=IX[i]
        xs=shp[i,0]*xsj
        Proj[k,0]+=xs
        for j in range(9):
            Proj[k,1+j]+=val[j]*xs
############################################################
### user element for linear elastic problem
############################################################
def elmt(IX,elCoords,elU,Proj):
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

    E=10.0e6;nu=0.3 # Young's modulus and poisson ratio
    D=np.zeros((3,3)) # Constitutive law

    term1=E/(1-nu**2);term2=E*nu/(1-nu**2)
    D[0,0]=term1;D[0,1]=term2
    D[1,0]=term2;D[1,1]=term1
    D[2,2]=0.5*E/(1+nu)

    for gpInd in range(Lint):
        xi =gs[gpInd,1]
        eta=gs[gpInd,2]
        shp,xsj=Shp2D(elCoords,xi,eta)
        JxW=gs[gpInd,0]*xsj

        B=np.zeros((3,2*nNodes))
        for i in range(nNodes):
            B[0,2*i  ]=shp[i,1]
            B[0,2*i+1]=0.0

            B[1,2*i  ]=0.0
            B[1,2*i+1]=shp[i,2]

            B[2,2*i  ]=shp[i,2]
            B[2,2*i+1]=shp[i,1]
        Bt=B.transpose()
        S=np.dot(D,B)
        strain=np.zeros(3)
        stress=np.zeros(3)
        
        for i in range(nDofs):
            strain[0]+=B[0,i]*elU[i]
            strain[1]+=B[1,i]*elU[i]
            strain[2]+=B[2,i]*elU[i]

            stress[0]+=S[0,i]*elU[i]
            stress[1]+=S[1,i]*elU[i]
            stress[2]+=S[2,i]*elU[i]
        #############################
        ### For projection value
        # For stress
        val[1-1]=stress[1-1]
        val[2-1]=stress[2-1]
        val[3-1]=stress[3-1]
        # For strain
        val[4-1]=strain[1-1]
        val[5-1]=strain[2-1]
        val[6-1]=strain[3-1]
        # For von Mises
        val[7-1]=np.sqrt(stress[0]**2+stress[1]**2+3*stress[2]**2-stress[0]*stress[1])
        val[8-1]=(stress[0]+stress[1])/2.0


        # Do projection from gauss point to nodal point
        Projection(nNodes,IX,xsj,shp,val,Proj)


        C=np.dot(Bt,S)
        for iInd in range(2*nNodes):
            # For residual
            for k in range(3):
                RHS[iInd]+=-Bt[iInd,k]*stress[k]*JxW
            for jInd in range(2*nNodes):
                K[iInd,jInd]+=C[iInd,jInd]*JxW
        
    return K,RHS
#######################################################
def FormKR(mesh,U,AMATRIX,RHS,Proj):
    nElmts=mesh.nElmts
    nNodesPerElmt=mesh.nNodesPerElmts
    nDofs=mesh.nNodes*2
    AMATRIX[:,:]=0.0 # [K]{u}=F-->here AMATRX is K
    RHS[:]=0.0
    Proj[:,:]=0.0

    for e in range(nElmts):
        elConn=mesh.Conn[e,:]-1
        elCoords=mesh.NodeCoords[elConn]
        elU=np.zeros(2*nNodesPerElmt)
        for i in range(nNodesPerElmt):
            elU[2*i  ]=U[2*elConn[i]  ]
            elU[2*i+1]=U[2*elConn[i]+1]

        k,rhs=elmt(elConn,elCoords,elU,Proj)
        

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
##########################################################
def PlotDisp(mesh,U,Flag):
    x=mesh.NodeCoords[:,0]
    y=mesh.NodeCoords[:,1]
    ux=U[0::2]
    uy=U[1::2]

    plt.figure(1)
    plt.title('$u_{x}$',fontsize=16)
    plt.tricontourf(x,y,ux,60)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.colorbar()
    plt.clim(np.min(ux),np.max(ux))

    plt.figure(2)
    plt.title('$u_{y}$',fontsize=16)
    plt.tricontourf(x,y,uy,60)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.colorbar()
    plt.clim(np.min(uy),np.max(uy))

    if Flag==True:
        plt.show()

def PlotStressStrain(mesh,Proj,Flag):
    x=mesh.NodeCoords[:,0]
    y=mesh.NodeCoords[:,1]
    sxx=Proj[:,1]
    sxy=Proj[:,2]
    syy=Proj[:,3]

    stxx=Proj[:,4]
    stxy=Proj[:,5]
    styy=Proj[:,6]

    vonMises=Proj[:,7]


    plt.figure(3)
    plt.title('$\sigma_{xx}$',fontsize=16)
    plt.tricontourf(x,y,sxx,60)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.colorbar()
    plt.clim(np.min(sxx),np.max(sxx))

    plt.figure(4)
    plt.title('$\sigma_{xy}$',fontsize=16)
    plt.tricontourf(x,y,sxy,60)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.colorbar()
    plt.clim(np.min(sxy),np.max(sxy))

    plt.figure(5)
    plt.title('$\sigma_{yy}$',fontsize=16)
    plt.tricontourf(x,y,syy,60)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.colorbar()
    plt.clim(np.min(syy),np.max(syy))

    ### For strain
    plt.figure(6)
    plt.title('$\epsilon_{xx}$',fontsize=16)
    plt.tricontourf(x,y,stxx,60)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.colorbar()
    plt.clim(np.min(stxx),np.max(stxx))

    plt.figure(7)
    plt.title('$\epsilon_{xy}$',fontsize=16)
    plt.tricontourf(x,y,stxy,60)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.colorbar()
    plt.clim(np.min(stxy),np.max(stxy))

    plt.figure(8)
    plt.title('$\epsilon_{yy}$',fontsize=16)
    plt.tricontourf(x,y,styy,60)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.colorbar()
    plt.clim(np.min(styy),np.max(styy))

    plt.figure(9)
    plt.title('$vonMises$',fontsize=16)
    plt.tricontourf(x,y,vonMises,60)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.colorbar()
    plt.clim(np.min(vonMises),np.max(vonMises))

    if Flag==True:
        plt.show()
    

if __name__=='__main__':
    Welcome()
    nx=40;ny=20
    W=10.0;H=2.0
    mesh=Mesh2D(nx,ny,0.0,W,0.0,H,'quad4')
    mesh.CreateMesh()
    mesh.SplitBCMesh()
    
    nDofs=2*mesh.nNodes
    U=np.zeros(nDofs)

    ngp=2
    AMATRIX=np.zeros((nDofs,nDofs))
    RHS=np.zeros(nDofs)
    Proj=np.zeros((mesh.nNodes,10))



    iters=0;MaxIters=10;IsConvergent=False
    atol=1.0e-13;rtol=1.0e-10 # absolute error and relative error
    iters=0;IsConvergent=False
    while iters<MaxIters and (not IsConvergent):
        ApplyDispBC('right','ux',mesh,U,0.1)
        FormKR(mesh,U,AMATRIX,RHS,Proj)
        ## For constrain condition
        ApplyConstrainBC('left','ux',mesh,AMATRIX,RHS)
        ApplyConstrainBC('bottom','uy',mesh,AMATRIX,RHS)
        ApplyConstrainBC('right','ux',mesh,AMATRIX,RHS)

        dU=np.linalg.solve(AMATRIX,RHS)
        U[:]+=dU[:]
        

        if iters==0:
            R0_norm=np.linalg.norm(RHS)
            dU0_norm=np.linalg.norm(dU)
        R_norm=np.linalg.norm(RHS)
        dU_norm=np.linalg.norm(dU)

        print('iteration=%2d, |R0|=%.5e <--->|R|=%.5e, |dU0|=%.5e <--->|dU=|%.5e'%(iters,R0_norm,R_norm,dU0_norm,dU_norm))
        if (R_norm<rtol*R0_norm or R_norm<atol) or (dU_norm<rtol*dU0_norm or dU_norm<atol):
            IsConvergent=True
            break
        iters+=1
    
    PlotDisp(mesh,U,False)
    PlotStressStrain(mesh,Proj,True)




