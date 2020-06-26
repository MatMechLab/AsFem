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
//+++ Date   : 2020.06.26
//+++ Purpose: implement the mesh import for msh(ver=2) file format
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Mesh/Gmsh2IO.h"

bool Gmsh2IO::ReadMeshFromFile(Mesh &mesh){
    if(mesh.GetDim()){
        return true;
    }
    return false;
}