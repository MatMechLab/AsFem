//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2022.05.06
//+++ Purpose: the mesh generator for AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Mesh/MeshGenerator.h"

MeshGenerator::MeshGenerator(){
    m_isMeshGenerated=false;
}

bool MeshGenerator::createMesh(const int &dim,const MeshType &meshtype,MeshData &meshdata){
    if(dim==1){
        m_isMeshGenerated=Lagrange1DMeshGenerator::generateMesh(meshtype,meshdata);
    }
    else if(dim==2){
        m_isMeshGenerated=Lagrange2DMeshGenerator::generateMesh(meshtype,meshdata);
    }
    else if(dim==3){
        m_isMeshGenerated=Lagrange3DMeshGenerator::generateMesh(meshtype,meshdata);
    }
    else{
        MessagePrinter::printErrorTxt("dim>3 is invalid for the mesh generation, error raised in MeshGenerator::generateMesh");
        MessagePrinter::exitAsFem();
    }

    return m_isMeshGenerated;

}