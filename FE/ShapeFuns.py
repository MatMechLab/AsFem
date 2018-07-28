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