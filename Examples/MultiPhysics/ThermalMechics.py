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
    print('*** Simple fem code for themo-mechanics  ***')
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
### user element for mechanics coupled with thermal conduct
############################################################
def elmt(IX,elCoords,elU,Proj,ctan):
    nNodes=np.size(elCoords,0)
    nDofs=3*nNodes # dofs: ux,uy,T

    ngp=2
    if nNodes>=8:
        ngp=3
    
    gs=Int2D(ngp)
    Lint=ngp*ngp

    K=np.zeros((nDofs,nDofs))
    RHS=np.zeros(nDofs)
    val=np.zeros(9)

    # all the parameters are normlized one
    E=120.0;nu=0.3   # Young's modulus and poisson ratio
    omega=0.08       # volume expansion coefficient
    conductivity=1.0 # thermal conductivity coefficient

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
        T=0.0;Tdot=0.0;gradT=np.zeros(3)
        for i in range(nNodes):
            B[0,2*i  ]=shp[i,1]
            B[0,2*i+1]=0.0

            B[1,2*i  ]=0.0
            B[1,2*i+1]=shp[i,2]

            B[2,2*i  ]=shp[i,2]
            B[2,2*i+1]=shp[i,1]
            
            T       +=shp[i,0]*elU[3*i+2][0]
            Tdot    +=shp[i,0]*elU[3*i+2][1]
            gradT[1]+=shp[i,1]*elU[3*i+2][0]
            gradT[2]+=shp[i,2]*elU[3*i+2][0]


        Bt=B.transpose()
        strain=np.zeros(3)# for total strain
        stress=np.zeros(3)# for strain free stress
        I=np.zeros(2)     # for eigen strain
        I[1-1]=1.0;I[2-1]=1.0
        

        for i in range(nNodes):
            # calulate mechanical strain
            strain[0]+=B[0,2*i]*elU[3*i][0]+B[0,2*i+1]*elU[3*i+1][0]
            strain[1]+=B[1,2*i]*elU[3*i][0]+B[1,2*i+1]*elU[3*i+1][0]
            strain[2]+=B[2,2*i]*elU[3*i][0]+B[2,2*i+1]*elU[3*i+1][0]
        
        eigen_strain=(omega*T/3.0)*I
        mechstrain=np.zeros(3)
        mechstrain[1-1]=strain[1-1]-eigen_strain[1-1]
        mechstrain[2-1]=strain[2-1]-eigen_strain[2-1]
        mechstrain[3-1]=strain[3-1]

        stress[1-1]=D[0,0]*mechstrain[0]+D[0,1]*mechstrain[1]+D[0,2]*mechstrain[2]
        stress[2-1]=D[1,0]*mechstrain[0]+D[1,1]*mechstrain[1]+D[1,2]*mechstrain[2]
        stress[3-1]=D[2,0]*mechstrain[0]+D[2,1]*mechstrain[1]+D[2,2]*mechstrain[2]
        
        dstressdT=np.zeros(3)
        dstressdT[1-1]=(D[0,0]+D[0,1])*(-omega/3.0)
        dstressdT[2-1]=(D[1,0]+D[1,1])*(-omega/3.0)
        dstressdT[3-1]=(D[2,0]+D[2,1])*(-omega/3.0)
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
        # For von Mises and hydrostatic stress
        val[7-1]=np.sqrt(stress[0]**2+stress[1]**2+3*stress[2]**2-stress[0]*stress[1])
        val[8-1]=(stress[0]+stress[1])/2.0


        # Do projection from gauss point to nodal point
        Projection(nNodes,IX,xsj,shp,val,Proj)


        C=np.dot(np.dot(Bt,D),B)
        for iInd in range(nNodes):
            # For residual
            for k in range(3):
                # R_ux
                RHS[3*iInd  ]+=-Bt[2*iInd  ,k]*stress[k]*JxW
                # R_uy
                RHS[3*iInd+1]+=-Bt[2*iInd+1,k]*stress[k]*JxW
            # R_t
            RHS[3*iInd+2]+=Tdot*shp[iInd,0]*JxW\
                          +conductivity*(gradT[1]*shp[iInd,1]+gradT[2]*shp[iInd,2])*JxW
            for jInd in range(nNodes):
                # Kux,ux
                K[3*iInd  ,3*jInd  ]+=C[2*iInd  ,2*jInd  ]*JxW
                # Kux,uy
                K[3*iInd  ,3*jInd+1]+=C[2*iInd  ,2*jInd+1]*JxW
                # Kuy,ux
                K[3*iInd+1,3*jInd  ]+=C[2*iInd+1,2*jInd  ]*JxW
                # Kuy,uy
                K[3*iInd+1,3*jInd+1]+=C[2*iInd+1,2*jInd+1]*JxW
                for k in range(3):
                    # Kux,c
                    K[3*iInd  ,3*jInd+2]+=Bt[2*iInd  ,k]*dstressdT[k]*shp[jInd,0]*JxW*ctan[0]
                    # Kuy,c
                    K[3*iInd+1,3*jInd+2]+=Bt[2*iInd+1,k]*dstressdT[k]*shp[jInd,0]*JxW*ctan[0]
                #######################################################
                #  for thermal conduct, here is just simple one way coupling
                #  temperature will influence stress, but stress has no effect on temperature    
                # Kt,tdot
                K[3*iInd+2,3*jInd+2]+=-shp[jInd,0]*shp[iInd,0]*JxW*ctan[1]
                # Kt,t
                K[3*iInd+2,3*jInd+2]+=-conductivity*(shp[jInd,1]*shp[iInd,1]+shp[jInd,2]*shp[iInd,2])*JxW*ctan[0]
        
    return K,RHS
#######################################################
def FormKR(mesh,U,V,AMATRIX,RHS,Proj,ctan):
    nElmts=mesh.nElmts
    nNodesPerElmt=mesh.nNodesPerElmts
    AMATRIX[:,:]=0.0 # [K]{u}=F-->here AMATRX is K
    RHS[:]=0.0
    Proj[:,:]=0.0

    for e in range(nElmts):
        elConn=mesh.Conn[e,:]-1
        elCoords=mesh.NodeCoords[elConn]
        elU=np.zeros((3*nNodesPerElmt,2))
        for i in range(nNodesPerElmt):
            elU[3*i  ][0]=U[3*elConn[i]  ]
            elU[3*i+1][0]=U[3*elConn[i]+1]
            elU[3*i+2][0]=U[3*elConn[i]+2]

            elU[3*i  ][1]=V[3*elConn[i]  ]
            elU[3*i+1][1]=V[3*elConn[i]+1]
            elU[3*i+2][1]=V[3*elConn[i]+2]

        k,rhs=elmt(elConn,elCoords,elU,Proj,ctan)
        

        # Assemble k,rhs to system matrix
        for i in range(nNodesPerElmt):
            iInd=elConn[i]
            RHS[3*iInd  ]+=rhs[3*i  ]
            RHS[3*iInd+1]+=rhs[3*i+1]
            RHS[3*iInd+2]+=rhs[3*i+2]
            for j in range(nNodesPerElmt):
                jInd=elConn[j]
                AMATRIX[3*iInd  ,3*jInd  ]+=k[3*i  ,3*j  ]
                AMATRIX[3*iInd  ,3*jInd+1]+=k[3*i  ,3*j+1]
                AMATRIX[3*iInd  ,3*jInd+2]+=k[3*i  ,3*j+2]

                AMATRIX[3*iInd+1,3*jInd  ]+=k[3*i+1,3*j  ]
                AMATRIX[3*iInd+1,3*jInd+1]+=k[3*i+1,3*j+1]
                AMATRIX[3*iInd+1,3*jInd+2]+=k[3*i+1,3*j+2]

                AMATRIX[3*iInd+2,3*jInd  ]+=k[3*i+2,3*j  ]
                AMATRIX[3*iInd+2,3*jInd+1]+=k[3*i+2,3*j+1]
                AMATRIX[3*iInd+2,3*jInd+2]+=k[3*i+2,3*j+2]
    
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
    elif dofname=='T':
        component=3
    else:
        sys.exit('dof name=%s in invalid!!!'%(dofname))
    
    ne=np.size(BCConn,0)
    nNodesPerBCElmt=np.size(BCConn,1)
    for e in range(ne):
        elConn=BCConn[e,:]-1
        for i in range(nNodesPerBCElmt):
            iInd=3*elConn[i]+component-1
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
    elif dofname=='T':
        component=3
    else:
        sys.exit('dof name=%s in invalid!!!'%(dofname))
    
    ne=np.size(BCConn,0)
    nNodesPerBCElmt=np.size(BCConn,1)
    penalty=1.0e16
    for e in range(ne):
        elConn=BCConn[e,:]-1
        for i in range(nNodesPerBCElmt):
            iInd=3*elConn[i]+component-1
            AMATRIX[iInd,iInd]=penalty
            RHS[iInd]=0.0
##########################################################
def PlotDisp(mesh,U,Flag):
    x=mesh.NodeCoords[:,0]
    y=mesh.NodeCoords[:,1]
    ux=U[0::3]
    uy=U[1::3]
    T =U[2::3]

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

    plt.figure(3)
    plt.title('Temperature',fontsize=16)
    plt.tricontourf(x,y,T,60)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.colorbar()
    plt.clim(np.min(T),np.max(T))

    plt.figure(4)
    plt.title('Temperature(deformed)',fontsize=16)
    factor=1.0
    xx=x+ux*factor
    yy=y+uy*factor
    plt.tricontourf(xx,yy,T,60)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.colorbar()
    plt.clim(np.min(T),np.max(T))
    plt.savefig('Temperature.png',dpi=800,bbox_inches='tight')

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
    hyStress=Proj[:,8]


    plt.figure(5)
    plt.title('$\sigma_{xx}$',fontsize=16)
    plt.tricontourf(x,y,sxx,60)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.colorbar()
    plt.clim(np.min(sxx),np.max(sxx))

    plt.figure(6)
    plt.title('$\sigma_{xy}$',fontsize=16)
    plt.tricontourf(x,y,sxy,60)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.colorbar()
    plt.clim(np.min(sxy),np.max(sxy))

    plt.figure(7)
    plt.title('$\sigma_{yy}$',fontsize=16)
    plt.tricontourf(x,y,syy,60)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.colorbar()
    plt.clim(np.min(syy),np.max(syy))

    ### For strain
    plt.figure(8)
    plt.title('$\epsilon_{xx}$',fontsize=16)
    plt.tricontourf(x,y,stxx,60)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.colorbar()
    plt.clim(np.min(stxx),np.max(stxx))

    plt.figure(9)
    plt.title('$\epsilon_{xy}$',fontsize=16)
    plt.tricontourf(x,y,stxy,60)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.colorbar()
    plt.clim(np.min(stxy),np.max(stxy))

    plt.figure(10)
    plt.title('$\epsilon_{yy}$',fontsize=16)
    plt.tricontourf(x,y,styy,60)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.colorbar()
    plt.clim(np.min(styy),np.max(styy))

    plt.figure(11)
    plt.title('$vonMises$',fontsize=16)
    plt.tricontourf(x,y,vonMises,60)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.colorbar()
    plt.clim(np.min(vonMises),np.max(vonMises))

    plt.figure(12)
    plt.title('$\sigma_{h}$',fontsize=16)
    plt.tricontourf(x,y,hyStress,60)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.colorbar()
    plt.clim(np.min(hyStress),np.max(hyStress))

    if Flag==True:
        plt.show()
    

if __name__=='__main__':
    Welcome()
    nx=50;ny=10
    W=10.0;H=2.0
    mesh=Mesh2D(nx,ny,0.0,W,0.0,H,'quad4')
    mesh.CreateMesh()
    mesh.SplitBCMesh()
    
    nDofs=3*mesh.nNodes

   
    AMATRIX=np.zeros((nDofs,nDofs))
    RHS=np.zeros(nDofs)
    Proj=np.zeros((mesh.nNodes,10))

    U0=np.zeros(nDofs)
    U=np.zeros(nDofs)
    V=np.zeros(nDofs)

    dt=5.0e-2;ctan=np.zeros(2)
    ctan[0]=1.0;ctan[1]=1.0/dt

    nSteps=20;MaxIters=50
    atol=1.0e-13;rtol=1.0e-10 # absolute error and relative error

    for step in range(nSteps):
        U[:]=U0[:]
        iters=0;IsConvergent=False
        while iters<MaxIters and (not IsConvergent):
            ApplyDispBC('left','ux',mesh,U,0.0)
            ApplyDispBC('left','uy',mesh,U,0.0)
            ApplyDispBC('bottom','T',mesh,U,2.0)
            ApplyDispBC('top','T',mesh,U,0.0)

            V[:]=(U[:]-U0[:])*ctan[1]
            FormKR(mesh,U,V,AMATRIX,RHS,Proj,ctan)
            ## For constrain condition
            ApplyConstrainBC('left','ux',mesh,AMATRIX,RHS)
            ApplyConstrainBC('left','uy',mesh,AMATRIX,RHS)
            ApplyConstrainBC('bottom','T',mesh,AMATRIX,RHS)
            ApplyConstrainBC('top','T',mesh,AMATRIX,RHS)

            dU=np.linalg.solve(AMATRIX,RHS)
            U[:]+=dU[:]

            R_norm=np.linalg.norm(RHS)
            dU_norm=np.linalg.norm(dU)
            if iters==0:
                R0_norm=R_norm
                dU0_norm=dU_norm

            if R_norm<rtol*R0_norm or R_norm<atol:
                IsConvergent=True
                break
            iters+=1
            #print('Step=%5d===> iters=%2d,|R0|=%12.5e,|R|=%12.5e,|dU0|=%12.5e,|dU|=%12.5e'%(step+1,iters,R0_norm,R_norm,dU0_norm,dU_norm))
        

        print('Step=%5d===> iters=%2d,|R0|=%12.5e,|R|=%12.5e,|dU0|=%12.5e,|dU|=%12.5e'%(step+1,iters,R0_norm,R_norm,dU0_norm,dU_norm))
        if IsConvergent:
            U0[:]=U[:]
    
    
    PlotDisp(mesh,U,False)
    PlotStressStrain(mesh,Proj,True)




