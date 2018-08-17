#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: FlyFox
@email: walkandthinker@gmail.com
@discuss: QQ group: 797998860
@Description: 1D heatconduction problem with Robin BC
"""
import numpy as np
import sys
import matplotlib.pyplot as plt

###############################
### For PyFem's own module
from Mesh import *
from ShapeFuns import *

def Welcome():
    print('********************************************')
    print('*** Welcome to use this simple script!   ***')
    print('*** Simple fem code for headconduct1D    ***')
    print('*** Author: FlyFox                       ***')
    print('*** Email: walkandthinker@gmail.com      ***')
    print('*** QQ group: 797998860                  ***')
    print('********************************************')


def Int1D(ngp):
    # generate 1D gausss integration point and weight
    if ngp<2 or ngp>5:
        sys.exit('Error: unsupported gauss point number(=%d)'%(ngp))
    gs=np.zeros((ngp,2))# 0->for w,1->xi
    if ngp == 2:
        gs[0,1] =-0.577350269189625764509148780502
        gs[1,1] = 0.577350269189625764509148780502
        gs[0,0] = 1.0
        gs[1,0] = 1.0
    elif ngp == 3:
        gs[0,0] = 0.555555555555555555555555555556
        gs[0,1] =-0.774596669241483377035853079956
        gs[1,0] = 0.888888888888888888888888888889
        gs[1,1] = 0.0
        gs[2,0] = 0.555555555555555555555555555556
        gs[2,1] = 0.774596669241483377035853079956
    return gs

def Int2D(ngp):
    gs=np.zeros((ngp*ngp,3))# 0->for w,1->for xi,2-> for eta
    k=0
    gs1d=Int1D(ngp)
    for i in range(ngp):
        for j in range(ngp):
            gs[k,0]=gs1d[i,0]*gs1d[j,0]
            gs[k,1]=gs1d[i,1]
            gs[k,2]=gs1d[j,1]
            k+=1
    return gs

#######################################################
### You can modify the user element for your own model
#######################################################
def elmt(elCoords,u,v,ctan):
    # Here F and A is just two constants for different solution of laplace equation
    nnodes=np.size(elCoords,0)
    Lint=nnodes
    ngp=2
    gs=Int2D(ngp)
    Lint=np.size(gs,0)
    conductivity=1.5 # heat conductivity coefficient

    k=np.zeros((nnodes,nnodes)) # for local k matrix
    rhs=np.zeros(nnodes)        # for local rhs
    gradT=np.zeros(2)

    for gp in range(Lint):
        xi=gs[gp,1]
        eta=gs[gp,2]
        shp,xsj=Shp2D(elCoords,xi,eta)
        JxW=xsj*gs[gp,0] # |J|*weight

        T=0.0;Tdot=0.0;gradT[0]=0.0;gradT[1]=0.0
        # for temperature, its gradient and dTdt on gauss point
        for i in range(nnodes):
            Tdot    +=v[i]*shp[i,0]
            T       +=u[i]*shp[i,0]
            gradT[0]+=u[i]*shp[i,1]
            gradT[1]+=u[i]*shp[i,2]  
        for I in range(nnodes):
            # for residual
            rhs[I]+=Tdot*shp[I,0]*JxW+conductivity*(gradT[0]*shp[I,1]+gradT[1]*shp[I,2])*JxW
            for J in range(nnodes):
                # for jacobain
                # For Kt,tdot part
                k[I,J]+=-shp[J,0]*shp[I,0]*ctan[1]*JxW
                # For Kt,t part
                k[I,J]+=-conductivity*(shp[J,1]*shp[I,1]+shp[J,2]*shp[I,2])*ctan[0]*JxW
    return k,rhs
#########################
def FormKR(NodeCoords,Conn,U,V,ctan):
    nElmts=np.size(Conn,0)
    nNodesPerElmt=np.size(Conn,1)
    nDofs=np.size(NodeCoords,0)
    AMATRX=np.zeros((nDofs,nDofs)) # [K]{u}=F-->here AMATRX is K
    RHS=np.zeros(nDofs)

    for e in range(nElmts):
        elConn=Conn[e,:]
        elCoords=NodeCoords[elConn-1]
        elU=U[elConn-1]
        elV=V[elConn-1]

        k,rhs=elmt(elCoords,elU,elV,ctan)

        # Assemble k,rhs to system matrix
        for i in range(nNodesPerElmt):
            RHS[elConn[i]-1]+=rhs[i]
            for j in range(nNodesPerElmt):
                AMATRX[elConn[i]-1,elConn[j]-1]+=k[i,j]
    return AMATRX,RHS
####################################################################


####################################################################
### For surface integration(here is Robin, you can add Neumann case)
def SurfaceElmt(elCoords,u,v,ctan):
    Tf=0.1
    nnodes=np.size(elCoords)
    Lint=nnodes
    gs=Int1D(Lint)
    k=np.zeros((nnodes,nnodes)) # for local k matrix
    rhs=np.zeros(nnodes)        # for local rhs
    for gp in range(Lint):
        xi=gs[gp,1]
        shp,xsj=Shp1D(elCoords,xi)
        JxW=xsj*gs[gp,0] # |J|*weight
        T=0.0
        for i in range(nnodes):
            T    +=u[i]*shp[i,0]
        for I in range(nnodes):
            rhs[I]+=-(T-Tf)*shp[I,0]*JxW
            for J in range(nnodes):
                k[I,J]+=shp[J,0]*shp[I,0]*JxW
    
    return k,rhs
#######################################
def ApplyRobinBC(sidename,mesh,U,V,ctan,AMATRX,RHS):
    if sidename=='left':
        # For left edge
        BCConn=mesh.LeftBCConn
        component=2-1 # use Y-coordinates
    elif sidename=='right':
        BCConn=mesh.RightBCConn
        component=2-1 # use Y-coordinates
    elif sidename=='bottom':
        BCConn=mesh.BottomBCConn
        component=1-1 # use X-coordinates
    elif sidename=='top':
        BCConn=mesh.TopBCConn
        component=1-1 # use X-coordinates
    else:
        sys.exit('Side name=%s is invalid!!!'%(sidename))
    
    # Now we do surface integration, in 2D case, it is simple line element case
    ne=np.size(BCConn,0)
    nNodesPerBCElmt=np.size(BCConn,1)
    for e in range(ne):
        elConn=BCConn[e,:]-1
        elCoords=mesh.NodeCoords[elConn,component]
        elU=U[elConn]
        elV=V[elConn]
        k,rhs=SurfaceElmt(elCoords,elU,elV,ctan)
        # Assemble surface integration to global K and R
        for i in range(nNodesPerBCElmt):
                RHS[elConn[i]]+=rhs[i]
                for j in range(nNodesPerBCElmt):
                    AMATRX[elConn[i],elConn[j]]+=k[i,j]
    return AMATRIX,RHS
#################################################################################
def PlotT(mesh,T):
    # plot T field(contour plot)
    # for e in range(12):
    #     elConn=mesh.Conn[e,:]-1
    #     np.append(elConn,elConn[0])
    #     elCoordsX=mesh.NodeCoords[elConn,0]
    #     elCoordsY=mesh.NodeCoords[elConn,1]
    #     elT=T[elConn]
    #     plt.tricontourf(elCoordsX,elCoordsY,elT,60)
    x=mesh.NodeCoords[:,0]
    y=mesh.NodeCoords[:,1]
    plt.tricontourf(x,y,T,60)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.colorbar()
    plt.savefig('HeatConduct2D.png',dpi=800,bbox_inches='tight')
    plt.show()
def PlotLine(mesh,T):
    x=[];t=[]
    for i in range(mesh.nNodes):
        if mesh.NodeCoords[i,1]==0.5:
            x.append(mesh.NodeCoords[i,0])
            t.append(T[i])
    plt.plot(x,t)
    plt.show()

if __name__=='__main__':
    Welcome()
    # Mesh generation
    W=10.0;H=1.0 # Width and height
    mesh2d=Mesh2D(40,4,0.0,W,0.0,H,'quad4')
    mesh2d.CreateMesh()
    mesh2d.SplitBCMesh()
    
    print('Node number=%6d'%(mesh2d.nNodes))
    print('Element number=%6d'%(mesh2d.nElmts))
    print('Node number of each element=%6d'%(mesh2d.nNodesPerElmts))


    

    nDofs=mesh2d.nNodes
    print('Dof number=%6d'%(nDofs))

    mesh2d.PlotMesh()

    T0=np.zeros(nDofs)
    T=np.zeros(nDofs)
    V=np.zeros(nDofs)
    dt=1.0e-2
    ctan=np.zeros(2)
    ctan[0]=1.0;ctan[1]=1.0/dt # for time integration, here is backward euler

    # Apply initial condition
    t1=10.0;t2=5.0
    for i in range(mesh2d.nNodes):
        if mesh2d.NodeCoords[i,0]<=0.5*W:
            T0[i]=t1
        else:
            T0[i]=t2
    

    PlotT(mesh2d,T0)

    nSteps=60
    MaxIters=50
    atol=1.0e-12;rtol=1.0e-9 # absolute error and relative error
    for step in range(nSteps):
        iters=0
        T[:]=T0[:]
        V[:]=0.0
        
        iters=0;IsConvergent=False
        while iters<=MaxIters and (not IsConvergent):
            AMATRIX,RHS=FormKR(mesh2d.NodeCoords,mesh2d.Conn,T,V,ctan)
            AMATRIX,RHS=ApplyRobinBC('right',mesh2d,T,V,ctan,AMATRIX,RHS)
            dT=np.linalg.solve(AMATRIX,RHS)
            T[:]+=dT[:]
            V[:]=(T[:]-T0[:])*ctan[1]

            if iters==0:
                R0_norm=np.linalg.norm(RHS)
                dU0_norm=np.linalg.norm(dT)
        
            R_norm=np.linalg.norm(RHS)
            dU_norm=np.linalg.norm(dT)

            if (R_norm<rtol*R0_norm or R_norm<atol) or (dU_norm<rtol*dU0_norm or dU_norm<atol):
                IsConvergent=True
                print('Step=%4d===>iteration=%2d, |R|=%.5e, |dU=|%.5e'%(step+1,iters,R_norm,dU_norm))
                break
        
            iters+=1
        if IsConvergent:
            T0[:]=T[:]
            if step%10==0:
                print('step=%2d'%(step))
                PlotLine(mesh2d,T0)
    PlotT(mesh2d,T0)

    