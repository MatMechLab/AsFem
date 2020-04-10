//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "Mesh/Mesh.h"

bool Mesh::CreateMesh(){
    if(_IsBuiltInMesh){
        if(GetDim()==1){
            return Create1DMesh();
        }
        else if(GetDim()==2){
            return Create2DMesh();
        }
        else if(GetDim()==3){
            return Create3DMesh();
        }
        else{
            // PetscPrintf(PETSC_COMM_WORLD,"*************************************************************************\n");
            PetscPrintf(PETSC_COMM_WORLD,"*** Error: invalid dimension number for built-in mesh generation !!!  ***\n");
            // PetscPrintf(PETSC_COMM_WORLD,"*************************************************************************\n");
            return false;
        }
    }
    else{
        if(_IsGmshMesh){
            return ReadMeshFromGmsh();
        }
        else if(_IsAbaqusMesh){
            return ReadMeshFromAbaqus();
        }
        else{
            // PetscPrintf(PETSC_COMM_WORLD,"*************************************************************************\n");
            PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsuppoted mesh import(only gmsh is supported)        !!!  ***\n");
            // PetscPrintf(PETSC_COMM_WORLD,"*************************************************************************\n");
            return false;
        }
    }
    
    return false;
}