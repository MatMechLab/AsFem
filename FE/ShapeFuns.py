#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: FlyFox
@email: walkandthinker@gmail.com
@discuss: QQ group: 797998860
"""
import numpy as np
import sys
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
    elif nnodes==4:
        #Third order line element
        shp[0,0]=-(3.0*xi+1.0)*(3.0*xi-1.0)*(    xi-1.0)/16.0
        shp[1,0]= (3.0*xi+3.0)*(3.0*xi-1.0)*(3.0*xi-3.0)/16.0
        shp[2,0]=-(3.0*xi+3.0)*(3.0*xi+1.0)*(3.0*xi-3.0)/16.0
        shp[3,0]= (    xi+1.0)*(3.0*xi+1.0)*(3.0*xi-1.0)/16.0
        shp[0,1]=-27.0*xi*xi/16.0+9.0*xi/8.0+ 1.0/16.0
        shp[1,1]= 81.0*xi*xi/16.0-9.0*xi/8.0-27.0/16.0
        shp[2,1]=-81.0*xi*xi/16.0-9.0*xi/8.0+27.0/16.0
        shp[3,1]= 27.0*xi*xi/16.0+9.0*xi/8.0- 1.0/16.0   

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
### Define 2D shape function
#######################################################
def Shp2D(elCoords,xi,eta):
    nnodes=np.size(elCoords,0)

    if nnodes<3 or nnodes>9:
        sys.exit('*** Error: unsupported mesh type(node=%d)'%(nnodes))

    shp=np.zeros((nnodes,3)) # 0->for N
                             # 1->for dN/dx
                             # 2->for dN/dy

    if nnodes==4:
        #2D-4Nodes rectangle element
        shp[0,0]=(1.0-xi)*(1.0-eta)/4.0
        shp[1,0]=(1.0+xi)*(1.0-eta)/4.0
        shp[2,0]=(1.0+xi)*(1.0+eta)/4.0
        shp[3,0]=(1.0-xi)*(1.0+eta)/4.0

        shp[0,1]= (eta-1.0)/4.0
        shp[0,2]= (xi -1.0)/4.0

        shp[1,1]= (1.0-eta)/4.0
        shp[1,2]=-(1.0+xi )/4.0

        shp[2,1]= (1.0+eta)/4.0
        shp[2,2]= (1.0+xi )/4.0

        shp[3,1]=-(1.0+eta)/4.0
        shp[3,2]= (1.0-xi )/4.0
    elif nnodes==8:
        #2D-8Nodes rectangle element
        shp[0,0]=(1.0-xi)*(1.0-eta)*(-xi-eta-1.0)/4.0
        shp[1,0]=(1.0+xi)*(1.0-eta)*( xi-eta-1.0)/4.0
        shp[2,0]=(1.0+xi)*(1.0+eta)*( xi+eta-1.0)/4.0
        shp[3,0]=(1.0-xi)*(1.0+eta)*(-xi+eta-1.0)/4.0
        shp[4,0]=(1.0-xi*xi)*(1.0-eta    )/2.0
        shp[5,0]=(1.0+xi   )*(1.0-eta*eta)/2.0
        shp[6,0]=(1.0-xi*xi)*(1.0+eta    )/2.0
        shp[7,0]=(1.0-xi   )*(1.0-eta*eta)/2.0

        #derivatives over xi and eta
        shp[0,1]=(1.0-eta)*(2.0*xi+eta)/4.0
        shp[0,2]=(1.0-xi )*(xi+2.0*eta)/4.0

        shp[1,1]=(1.0-eta)*(2.0*xi-eta)/4.0
        shp[1,2]=(1.0+xi )*(2.0*eta-xi)/4.0

        shp[2,1]=(1.0+eta)*(2.0*xi+eta)/4.0
        shp[2,2]=(1.0+xi )*(xi+2.0*eta)/4.0

        shp[3,1]=(1.0+eta)*(2.0*xi-eta)/4.0
        shp[3,2]=(1.0-xi )*(2.0*eta-xi)/4.0

        shp[4,1]=xi*(eta-1.0)
        shp[4,2]=(xi*xi-1.0)/2.0

        shp[5,1]=(1.0-eta*eta)/2.0
        shp[5,2]=-(1.0+xi)*eta

        shp[6,1]=-xi*(1.0+eta)
        shp[6,2]=(1.0-xi*xi)/2.0

        shp[7,1]=(eta*eta-1.0)/2.0
        shp[7,2]=(xi-1.0)*eta
    elif nnodes==9:
        #2D-9Nodes rectangle element
        shp[0,0]=(xi*xi-xi )*(eta*eta-eta)/4.0
        shp[1,0]=(xi*xi+xi )*(eta*eta-eta)/4.0
        shp[2,0]=(xi*xi+xi )*(eta*eta+eta)/4.0
        shp[3,0]=(xi*xi-xi )*(eta*eta+eta)/4.0
        shp[4,0]=(1.0-xi*xi)*(eta*eta-eta)/2.0
        shp[5,0]=(xi*xi+xi )*(1.0-eta*eta)/2.0
        shp[6,0]=(1.0-xi*xi)*(eta*eta+eta)/2.0
        shp[7,0]=(xi*xi-xi )*(1.0-eta*eta)/2.0
        shp[8,0]=(1.0-xi*xi)*(1.0-eta*eta)

        shp[0,1]=(2.0*xi-1.0)*(eta*eta-eta)/4.0
        shp[0,2]=(xi*xi-xi  )*(2.0*eta-1.0)/4.0

        shp[1,1]=(2.0*xi+1.0)*(eta*eta-eta)/4.0
        shp[1,2]=(xi*xi+xi  )*(2.0*eta-1.0)/4.0

        shp[2,1]=(2.0*xi+1.0)*(eta*eta+eta)/4.0
        shp[2,2]=(xi*xi+xi  )*(2.0*eta+1.0)/4.0

        shp[3,1]=(2.0*xi-1.0)*(eta*eta+eta)/4.0
        shp[3,2]=(xi*xi-xi  )*(2.0*eta+1.0)/4.0

        shp[4,1]=-xi*(eta*eta-eta)
        shp[4,2]=(1.0-xi*xi )*(2.0*eta-1.0)/2.0

        shp[5,1]=(2.0*xi+1.0)*(1.0-eta*eta)/2.0
        shp[5,2]=-(xi*xi+xi )*eta

        shp[6,1]=-xi*(eta*eta+eta)
        shp[6,2]=(1.0-xi*xi )*(2.0*eta+1.0)/2.0

        shp[7,1]=(2.0*xi-1.0)*(1.0-eta*eta)/2.0
        shp[7,2]=-(xi*xi-xi )*eta

        shp[8,1]=-2.0*xi*(1.0-eta*eta)
        shp[8,2]=-2.0*eta*(1.0-xi*xi)
    
    #compute jacob transform matrix
    dxdxi=0.0;dxdeta=0.0
    dydxi=0.0;dydeta=0.0
    for i in range(nnodes):
        dxdxi +=elCoords[i,0]*shp[i,1]
        dxdeta+=elCoords[i,0]*shp[i,2]

        dydxi +=elCoords[i,1]*shp[i,1]
        dydeta+=elCoords[i,1]*shp[i,2]
    
    Jac =np.zeros((2,2))
    Jac[0,0]=dxdxi ;Jac[0,1]=dydxi
    Jac[1,0]=dxdeta;Jac[1,1]=dydeta
    DetJac=Jac[0,0]*Jac[1,1]-Jac[0,1]*Jac[1,0]
    if np.abs(DetJac)<1.0e-12:
        sys.exit('*** Error: singular element!!!')
    
    XJac=np.inv(Jac)

    for i in range(nnodes):
        temp    =XJac[0,0]*shp[i,1]+XJac[0,1]*shp[i,2]
        shp[i,2]=XJac[1,0]*shp[i,1]+XJac[1,1]*shp[i,2]
        shp[i,1]=temp


