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


def Welcome():
    print('********************************************')
    print('*** Welcome to use this simple script!   ***')
    print('*** Simple fem code for headconduct1D    ***')
    print('*** Author: FlyFox                       ***')
    print('*** Email: walkandthinker@gmail.com      ***')
    print('*** QQ group: 797998860                  ***')
    print('********************************************')


conductivity=1.0 # heat conductivity coefficient
L=10.0 # Length of 1D domain
ne=100;P=2  # element number and mesh order    

#######################################################
### Create 1D lagrange mesh
#######################################################
def CreateMesh(xmin,xmax,ne,p):
    if xmin>=xmax:
        sys.exit('Invalid mesh information xmin=%.5f >= xmax=%.5f'%(xmin,xmax))
    if ne<=1:
        sys.exit('ne=%d is too less for fem calculation!'%(ne))
    if p<1 or p>3:
        sys.exit('Invalid element order(p=%d)!'%(p))
    
    nNodesPerElmt=p+1
    nElmts=ne
    nNodes=nElmts*p+1

    dx=(xmax-xmin)/(nNodes-1)

    NodeCoords=np.zeros(nNodes)
    Conn=np.zeros((nElmts,nNodesPerElmt),dtype=np.int)

    # For node's coordinates
    for i in range(nNodes):
        NodeCoords[i]=xmin+i*dx
    
    # For element's connectivity
    # 1D lagrange mesh should looks like: 1-------2          // for p=1
    #                                     1--2--3 or 1--3--2 // for p=2(here we use the first one)
    for e in range(1,nElmts+1):
        for j in range(1,nNodesPerElmt+1):
            Conn[e-1,j-1]=(e-1)*p+j
    
    return NodeCoords,Conn
    

    print('*** Mesh generation finished!            ***')
    print('***--------------------------------------***')

#######################################################
### Define 1D gauss integration rule
#######################################################
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

#######################################################
### Define 1D shape function
#######################################################
def Shp1D(elCoords,xi):
    nnodes=np.size(elCoords)
    shp=np.zeros((nnodes,2)) # 0->for N
                             # 1->for dN/dx
    # for 1D lagrange mesh, one is referred to:
    # https://en.wikipedia.org/wiki/Lagrange_polynomial                         
    if nnodes == 2:
        #Linear element
        shp[0,0] = 0.5*(1.0 - xi)
        shp[0,1] =-0.5

        shp[1,0] = 0.5*(xi + 1.0)
        shp[1,1] = 0.5
    elif nnodes == 3:
        # Quadratic line element(3 nodes)
        shp[0,0] = 0.5*xi*(xi - 1.0)
        shp[0,1] = 0.5*(2.0*xi - 1.0)

        shp[1,0] = -(xi + 1.0)*(xi - 1.0)
        shp[1,1] = -2.0*xi

        shp[2,0] = 0.5*xi*(xi + 1.0)
        shp[2,1] = 0.5*(2.0*xi + 1.0)

    # Now we calculate the |J| of current mesh
    dxdxi=0.0
    for i in range(nnodes):
        dxdxi+=elCoords[i]*shp[i,1]
    DetJac=np.abs(dxdxi)
    if DetJac<1.0e-13:
        sys.exit('Error: current mesh is singular!!!')
    # now the dN/dx is still in local coordinate(xi,eta),
    # we should return back to nature coordinate(x,y)                

    for i in range(nnodes):
        shp[i,1]/=DetJac
    
    return shp,DetJac

#######################################################
### You can modify the user element for your own model
#######################################################
def elmt(elCoords,u,v,ctan):
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

        T=0.0;gradT=0.0;Tdot=0.0 # for temperature, its gradient and dTdt on gauss point
        for i in range(nnodes):
            Tdot +=v[i]*shp[i,0]
            T    +=u[i]*shp[i,0]
            gradT+=u[i]*shp[i,1]  


        for I in range(nnodes):
            # for residual
            rhs[I]+=Tdot*shp[I,0]*JxW+conductivity*gradT*shp[I,1]*JxW
            for J in range(nnodes):
                # for jacobain
                # For Kt,tdot part
                k[I,J]+=-shp[J,0]*shp[I,0]*ctan[1]*JxW
                # For Kt,t part
                k[I,J]+=-conductivity*shp[J,1]*shp[I,1]*ctan[0]*JxW
    
    return k,rhs

def FormKR(NodeCoords,Conn,U,V,ctan):
    nElmts=np.size(Conn,0)
    nNodesPerElmt=np.size(Conn,1)
    nDofs=np.size(NodeCoords)
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
    
    ################################################
    ### Now we apply the robin bc
    ### since it is just 1D case, user can directly apply it on node
    ### for 2D and 3D cases, surface integration is required!
    ################################################
    T=U[-1];Tf=1.0 # for grad(T)*n=T-Tf
    RHS[-1]+=-(T-Tf)
    AMATRX[-1,-1]+=1



    #print('*** System equation: Ax=F generated!     ***')
    #print('***--------------------------------------***')
    return AMATRX,RHS


if __name__=='__main__':
    Welcome()

    temperature1=10.0;temperature2=2.0

    NodeCoords,Conn=CreateMesh(0.0,L,ne,P)
    
    nDofs=np.size(NodeCoords)
    T0=np.zeros(nDofs)
    T=np.zeros(nDofs)
    dT=np.zeros(nDofs)
    V=np.zeros(nDofs)
    ctan=np.zeros(2)

    for i in range(nDofs):
        if i<nDofs/2:
            T0[i]=temperature1
        else:
            T0[i]=temperature2


    dt=1.0e-2 # time step size
    ctan[0]=1.0;ctan[1]=1.0/dt # for time integration
    nstep=100

    MaxIters=50
    atol=1.0e-12;rtol=1.0e-9 # absolute error and relative error
    for step in range(nstep):
        T[:]=T0[:]
        V[:]=0.0
        
        iters=0;IsConvergent=False
        while iters<=MaxIters and (not IsConvergent):
            AMTRX,RHS=FormKR(NodeCoords,Conn,T,V,ctan)
            dT=np.linalg.solve(AMTRX,RHS)
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
                plt.plot(NodeCoords,T0)
    plt.savefig('HeatConduct1D.png',dpi=500,bbox_inches='tight')
    plt.show()
    
    

