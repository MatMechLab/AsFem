//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.06.29
//+++ Purpose: Implement the mesh generation for normal lagrange
//+++          mesh, for exmaple:
//+++             1D-->edge2,edge3,edge4,edge5
//+++             2D-->quad4,quad8,quad9
//+++             3D-->hex8 ,hex20,hex27
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Mesh/LagrangeMesh.h"

bool LagrangeMesh::CreateLagrangeMesh(){
    if(GetBulkMeshDim()==1){
        return Create1DLagrangeMesh();
    }
    else if(GetBulkMeshDim()==2){
        return Create2DLagrangeMesh();
    }
    else if(GetBulkMeshDim()==3){
        return Create3DLagrangeMesh();
    }
    else{
        MessagePrinter::PrintErrorTxt("unsupported dim(>3) for mesh generation");
        return false;
    }
    return false;
}