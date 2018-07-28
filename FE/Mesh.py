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

    
    def PlotMesh(self):
        y=np.zeros(self.nNodesPerElmts)
        for e in range(self.nElmts):
            elConn=self.Conn[e,:]-1
            x=self.NodeCoords[elConn]
            plt.plot(x,y,'-+')
        
        plt.show()



    

