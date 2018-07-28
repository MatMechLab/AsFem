#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: FlyFox
@email: walkandthinker@gmail.com
@discuss: QQ group: 797998860
"""

import numpy as np
import matplotlib.pyplot as plt
import sys 

class Mesh1D:
    def __init__(self,nx,xmin,xmax,p):
        if nx<1:
            sys.exit('*** Error: nx=%d is invalid !!!'%(nx))
        if xmin>=xmax:
            sys.exit('*** Error: xmin=%g >= xmax=%g !!!'%(xmin,xmax)) 

        if p<1 or p>3:
            sys.exit('*** Error: unsupported element order(p=%d) !!!'%(p))   
        self.nElmts=nx
        self.Xmin=xmin 
        self.Xmax=xmax
        self.P=p
    
    def CreateMesh(self):
        print('***--------------------------------------***')
        print('*** Start to create 1D mesh...           ***')
        self.nNodesPerElmts=self.P+1
        self.nNodes=self.nElmts*self.P+1
        self.NodeCoords=np.zeros(self.nNodes)
        self.Conn=np.zeros((self.nElmts,self.nNodesPerElmts),dtype=np.int)

        dx=(self.Xmax-self.Xmin)/(self.nNodes-1)
        for i in range(self.nNodes):
            self.NodeCoords[i]=self.Xmin+i*dx
        
        for e in range(1,self.nElmts+1):
            for j in range(1,self.nNodesPerElmts+1):
                self.Conn[e-1,j-1]=(e-1)*self.P+j

        print('*** Mesh1D generation finished!          ***')
        print('***--------------------------------------***')
    
    def GetBCConn(self,sidename):
        if sidename=='left':
            return self.Conn[0,0]
        elif sidename=='right':
            return self.Conn[self.nElmts-1,-1]
        else:
            sys.exit('*** Error: wrong side name')

    def GetIthElCoords(self,e):
        if e<1 or e>self.nElmts:
            sys.exit('*** Error: e=%d is out of range(nElmts=)!!!'%(e,self.nElmts))
        elConn=self.Conn[e-1,:]
        return self.NodeCoords[elConn]

    def GetIthElConn(self,e):
        if e<1 or e>self.nElmts:
            sys.exit('*** Error: e=%d is out of range(nElmts=)!!!'%(e,self.nElmts))
        return self.Conn[e-1,:]

    
    def PlotMesh(self):
        y=np.zeros(self.nNodesPerElmts)
        for e in range(self.nElmts):
            elConn=self.Conn[e,:]-1
            x=self.NodeCoords[elConn]
            plt.plot(x,y,'-+')
        
        plt.show()

#################################################
### For 2D Mesh
#################################################
class Mesh2D:
    def __init__(self,nx,ny,xmin,xmax,ymin,ymax,meshtype):
        if nx<1:
            sys.exit('*** Error: nx=%d is invalid!!!'%(nx))
        self.Nx=nx 

        if ny<1:
            sys.exit('*** Error: ny=%d is invalid!!!'%(ny))
        self.Ny=ny

        if xmin>=xmax:
            sys.exit('*** Error: xmin=%g >= xmax=%g !!!'%(xmin,xmax))
        self.Xmin=xmin
        self.Xmax=xmax

        if ymin>=ymax:
            sys.exit('*** Error: ymin=%g >= ymax=%g !!!'%(ymin,ymax))
        self.Ymin=ymin
        self.Ymax=ymax

        
        if (meshtype != 'quad4') and (meshtype !='quad8') and (meshtype !='quad9'):
            sys.exit('*** Error: unsupported mesh type !!!')
        self.MeshType=meshtype
    
    def CreateMesh(self):
        if self.MeshType=='quad4':
            dx=(self.Xmax-self.Xmin)/self.Nx
            dy=(self.Ymax-self.Ymin)/self.Ny

            self.nElmts=self.Nx*self.Ny
            self.nNodes=(self.Nx+1)*(self.Ny+1)
            self.nNodesPerElmts=4

            self.NodeCoords=np.zeros((self.nNodes,2))
            self.Conn=np.zeros((self.nElmts,self.nNodesPerElmts),dtype=np.int)

            for j in range(1,self.Ny+1+1):
                for i in range(1,self.Nx+1+1):
                    k=(j-1)*(self.Nx+1)+i
                    self.NodeCoords[k-1,0]=self.Xmin+(i-1)*dx
                    self.NodeCoords[k-1,1]=self.Ymin+(j-1)*dy
            
            for j in range(1,self.Ny+1):
                for i in range(1,self.Nx+1):
                    e=(j-1)*self.Nx+i
                    i1=(j-1)*(self.Nx+1)+i
                    i2=i1+1
                    i3=i2+self.Nx+1
                    i4=i3-1
                    self.Conn[e-1,0]=i1
                    self.Conn[e-1,1]=i2
                    self.Conn[e-1,2]=i3
                    self.Conn[e-1,3]=i4
        elif self.MeshType=='quad9':
            dx=(self.Xmax-self.Xmin)/(2.0*self.Nx)
            dy=(self.Ymax-self.Ymin)/(2.0*self.Ny)

            self.nElmts=self.Nx*self.Ny
            self.nNodes=(2*self.Nx+1)*(2*self.Ny+1)-self.nElmts
            self.nNodesPerElmts=8

            self.NodeCoords=np.zeros((self.nNodes,2))
            self.Conn=np.zeros((self.nElmts,self.nNodesPerElmts),dtype=np.int)

            for j in range(self.Ny+1):
                for i in range(2*self.Nx+1+1):
                    k=(j-1)*(2*self.Nx+1+self.Nx+1)+i
                    self.NodeCoords[k-1,0]=self.Xmin+(i-1)*dx 
                    self.NodeCoords[k-1,1]=self.Ymin+(j-1)*2*dy

                for i in range(self.Nx+1+1):
                    k=(j-1)*(2*self.Nx+1+self.Nx+1)+2*self.Nx+1+i
                    self.NodeCoords[k-1,0]=self.Xmin+(i-1)*2*dx
                    self.NodeCoords[k-1,1]=self.Ymin+(j-1)*2*dy+dy
                    
    
    def SplitBCMesh(self):


    def PlotMesh(self):
        plt.hold
        for e in range(self.nElmts):
            elConn=self.Conn[e,:]-1
            np.append(elConn,elConn[0])
            x=self.NodeCoords[elConn,0]
            y=self.NodeCoords[elConn,1]
            plt.plot(x,y,'-k+')
            i=np.array([elConn[-1],elConn[0]])
            x=self.NodeCoords[i,0]
            y=self.NodeCoords[i,1]
            plt.plot(x,y,'-k+')
        plt.show()
        




    

