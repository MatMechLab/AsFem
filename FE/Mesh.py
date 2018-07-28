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
        elif self.MeshType=='quad8':
            dx=(self.Xmax-self.Xmin)/(2.0*self.Nx)
            dy=(self.Ymax-self.Ymin)/(2.0/self.Ny)
            
            self.nElmts=self.Nx*self.Ny
            self.nNodes=(2*self.Nx+1)*(2*self.Ny+1)-self.nElmts
            self.nNodesPerElmts=8
            self.NodeCoords=np.zeros((self.nNodes,2))
            self.Conn=np.zeros((self.nElmts,self.nNodesPerElmts))

            for j in range(1,self.Ny+1):
                for i in range(1,2*self.Nx+1+1):
                    k=(j-1)*(2*self.Nx+1+self.Nx+1)+i
                    self.NodeCoords[k-1,0]=self.Xmin+(i-1)*dx
                    self.NodeCoords[k-1,1]=self.Ymin+(j-1)*2*dy
                for i in range(1,self.Nx+1+1):
                    k=(j-1)*(2*self.Nx+1+self.Nx+1)+2*self.Nx+1+i
                    self.NodeCoords[k-1,0]=self.Xmin+(i-1)*2*dx
                    self.NodeCoords[k-1,1]=self.Ymin+(j-1)*2*dy+dy 
            j=self.Ny+1
            for i in range(1,2*self.Nx+1+1):
                k=(j-1)*(2*self.Nx+1+self.Nx+1)+i
                self.NodeCoords[k-1,0]=self.Xmin+(i-1)*dx
                self.NodeCoords[k-1,1]=self.Ymin+(j-1)*2*dy
            for j in range(1,self.Ny+1):
                for i in range(1,self.Nx+1):
                    e=(j-1)*self.Nx+i
                    i1=(j-1)*(2*self.Nx+1+self.Nx+1)+2*i-1
                    i2=i1+2
                    i3=i2+(2*self.Nx+1+self.Nx+1)
                    i4=i3-2
                    i5=i1+1
                    i6=i2+(2*self.Nx+1)-i
                    i7=i3-1
                    i8=i1+(2*self.Nx+1)-(i-1)

                    self.Conn[e-1,1-1]=i1
                    self.Conn[e-1,2-1]=i2
                    self.Conn[e-1,3-1]=i3
                    self.Conn[e-1,4-1]=i4
                    self.Conn[e-1,5-1]=i5
                    self.Conn[e-1,6-1]=i6
                    self.Conn[e-1,7-1]=i7
                    self.Conn[e-1,8-1]=i8
        elif self.MeshType=='quad9':
            dx=(self.Xmax-self.Xmin)/(2.0*self.Nx)
            dy=(self.Ymax-self.Ymin)/(2.0*self.Ny)

            self.nElmts=self.Nx*self.Ny
            self.nNodes=(2*self.Nx+1)*(2*self.Ny+1)
            self.nNodesPerElmts=9
            self.NodeCoords=np.zeros((self.nNodes,2))
            self.Conn=np.zeros((self.nElmts,self.nNodesPerElmts),dtype=np.int)

            for j in range(1,2*self.Ny+1+1):
                for i in range(1,2*self.Nx+1+1):
                    k=(j-1)*(2*self.Nx+1)+i
                    self.NodeCoords[k-1,0]=self.Xmin+(i-1)*dx 
                    self.NodeCoords[k-1,1]=self.Ymin+(j-1)*dy
            
            for j in range(1,self.Ny+1):
                for i in range(1,self.Nx+1):
                    e=(j-1)*self.Nx+i
                    i1=(j-1)*2*(2*self.Nx+1)+2*i-1
                    i2=i1+2
                    i3=i2+2*(2*self.Nx+1)
                    i4=i3-2

                    i5=i1+1
                    i6=i2+(2*self.Nx+1)
                    i7=i3-1
                    i8=i1+(2*self.Nx+1)
                    i9=i8+1

                    self.Conn[e-1,1-1]=i1
                    self.Conn[e-1,2-1]=i2
                    self.Conn[e-1,3-1]=i3
                    self.Conn[e-1,4-1]=i4
                    self.Conn[e-1,5-1]=i5
                    self.Conn[e-1,6-1]=i6
                    self.Conn[e-1,7-1]=i7
                    self.Conn[e-1,8-1]=i8
                    self.Conn[e-1,9-1]=i9
        ###################################################
        #split boundary mesh
    
    def SplitBCMesh(self):
        if self.MeshType=='quad4':
            nBCNodes=2
        else:
            nBCNodes=3
        self.LeftBCConn=np.zeros((self.Ny,nBCNodes),dtype=int)
        self.RightBCConn=np.zeros((self.Ny,nBCNodes),dtype=int)
        self.BottomBCConn=np.zeros((self.Nx,nBCNodes),dtype=int)
        self.TopBCConn=np.zeros((self.Nx,nBCNodes),dtype=int)

        k1=0;k2=0
        for j in range(1,self.Ny+1):
            i1=1;i2=self.Nx
            e1=(j-1)*self.Nx+i1
            e2=(j-1)*self.Nx+i2
            if self.MeshType=='quad4':
                # 4---3
                # |   |
                # 1---2
                self.LeftBCConn[k1,0]=self.Conn[e1-1,4-1]
                self.LeftBCConn[k1,1]=self.Conn[e1-1,1-1]

                self.RightBCConn[k2,0]=self.Conn[e2-1,2-1]
                self.RightBCConn[k2,1]=self.Conn[e2-1,3-1]
            else:
                # 4--7--3
                # |     |
                # 8     6
                # |     |
                # 1--5--2
                self.LeftBCConn[k1,0]=self.Conn[e1-1,4-1]
                self.LeftBCConn[k1,1]=self.Conn[e1-1,8-1]
                self.LeftBCConn[k1,2]=self.Conn[e1-1,1-1]

                self.RightBCConn[k2,0]=self.Conn[e2-1,2-1]
                self.RightBCConn[k2,1]=self.Conn[e2-1,6-1]
                self.RightBCConn[k2,2]=self.Conn[e2-1,3]
            k1+=1
            k2+=1
        k1=0;k2=0
        for i in range(1,self.Nx+1):
            j1=1;j2=self.Ny
            e1=(j1-1)*self.Nx+i
            e2=(j2-1)*self.Nx+i
            if self.MeshType=='quad4':
                # 4---3
                # |   |
                # 1---2
                self.BottomBCConn[k1,1-1]=self.Conn[e1-1,1-1]
                self.BottomBCConn[k1,2-1]=self.Conn[e1-1,2-1]

                self.TopBCConn[k2,1-1]=self.Conn[e2-1,3-1]
                self.TopBCConn[k2,2-1]=self.Conn[e2-1,4-1]
            else:
                # 4--7--3
                # |     |
                # 8     6
                # |     |
                # 1--5--2
                self.BottomBCConn[k1,1-1]=self.Conn[e1-1,1-1]
                self.BottomBCConn[k1,2-1]=self.Conn[e1-1,5-1]
                self.BottomBCConn[k1,3-1]=self.Conn[e1-1,2-1]

                self.TopBCConn[k2,1-1]=self.Conn[e2-1,3-1]
                self.TopBCConn[k2,2-1]=self.Conn[e2-1,7-1]
                self.TopBCConn[k2,3-1]=self.Conn[e2-1,4-1]
            k1+=1
            k2+=1

    def PlotMesh(self):
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
        for e in range(self.Nx):
            elConn=self.BottomBCConn[e,:]-1
            x=self.NodeCoords[elConn,0]
            y=self.NodeCoords[elConn,1]
            plt.plot(x,y,'-r+')
            elConn=self.TopBCConn[e,:]-1
            x=self.NodeCoords[elConn,0]
            y=self.NodeCoords[elConn,1]
            plt.plot(x,y,'-g+')
        for e in range(self.Ny):
            elConn=self.LeftBCConn[e,:]-1
            x=self.NodeCoords[elConn,0]
            y=self.NodeCoords[elConn,1]
            plt.plot(x,y,'-b+')
            elConn=self.RightBCConn[e,:]-1
            x=self.NodeCoords[elConn,0]
            y=self.NodeCoords[elConn,1]
            plt.plot(x,y,'-y+')

        plt.show()
        




    

