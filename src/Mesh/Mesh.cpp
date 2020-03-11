//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "Mesh/Mesh.h"

Mesh::Mesh(){
    _IsMeshCreated=false;
    _NodeCoords.clear();
    _ElmtVTKCellType.clear();
    _ElmtConn.clear();
    _BulkElmtConn.clear();
    _LineElmtConn.clear();_SurfaceElmtConn.clear();
    _BulkMeshType=MeshType::EDGE2;
    _SurfaceMeshType=MeshType::NULLTYPE;
    _LineMeshType=MeshType::NULLTYPE;
    _IsBuiltInMesh=true;_IsGmshMesh=false;

    _MeshFileName.clear();
    _GmshFileName.clear();

    _Nx=2;_Ny=2;_Nz=2;
    _Xmin=0.0;_Xmax=1.0;
    _Ymin=0.0;_Ymax=1.0;
    _Zmin=0.0;_Zmax=1.0;

    _nMaxDim=1;_nMinDim=0;
    _nElmts=0;_nBulkElmts=0;_nNodes=0;
    _nOrder=0;

    _TotalVolume=0.0;
    _ElmtVolume.clear();

    _PhysicGroupNameList.clear();
    _PhysicNameToIDList.clear();
    _PhysicIDToNameList.clear();
    _PhysicNameToElmtIndexSet.clear();
    _nPhysicGroups=0;
    _BulkMeshTypeName.clear();
}
//*******************************************
void Mesh::Reset(){
    _IsMeshCreated=false;
    _ElmtVTKCellType.clear();
    _NodeCoords.clear();
    _ElmtConn.clear();
    _BulkElmtConn.clear();
    _LineElmtConn.clear();_SurfaceElmtConn.clear();
    _BulkMeshType=MeshType::EDGE2;
    _SurfaceMeshType=MeshType::NULLTYPE;
    _LineMeshType=MeshType::NULLTYPE;
    _IsBuiltInMesh=true;_IsGmshMesh=false;

    _MeshFileName.clear();
    _GmshFileName.clear();

    _Nx=2;_Ny=2;_Nz=2;
    _Xmin=0.0;_Xmax=1.0;
    _Ymin=0.0;_Ymax=1.0;
    _Zmin=0.0;_Zmax=1.0;

    _nMaxDim=1;_nMinDim=0;

    _nElmts=0;_nBulkElmts=0;_nNodes=0;
    _nOrder=0;

    _TotalVolume=0.0;
    _ElmtVolume.clear();

    _PhysicGroupNameList.clear();
    _PhysicNameToIDList.clear();
    _PhysicIDToNameList.clear();
    _PhysicNameToElmtIndexSet.clear();
    _nPhysicGroups=0;
    _BulkMeshTypeName.clear();
}
//*********************************************
void Mesh::SetMeshType(string meshtype){
    if(meshtype=="edge2"){
        _BulkMeshTypeName=meshtype;
        _BulkMeshType=MeshType::EDGE2;
    }
    else if(meshtype=="edge3"){
        _BulkMeshTypeName=meshtype;
        _BulkMeshType=MeshType::EDGE3;
    }
    else if(meshtype=="edge4"){
        _BulkMeshTypeName=meshtype;
        _BulkMeshType=MeshType::EDGE4;
    }
    else if(meshtype=="quad4"){
        _BulkMeshTypeName=meshtype;
        _BulkMeshType=MeshType::QUAD4;
    }
    else if(meshtype=="quad8"){
        _BulkMeshTypeName=meshtype;
        _BulkMeshType=MeshType::QUAD8;
    }
    else if(meshtype=="quad9"){
        _BulkMeshTypeName=meshtype;
        _BulkMeshType=MeshType::QUAD9;
    }
    else if(meshtype=="hex8"){
        _BulkMeshTypeName=meshtype;
        _BulkMeshType=MeshType::HEX8;
    }
    else if(meshtype=="hex20"){
        _BulkMeshTypeName=meshtype;
        _BulkMeshType=MeshType::HEX20;
    }
    else if(meshtype=="hex27"){
        _BulkMeshTypeName=meshtype;
        _BulkMeshType=MeshType::HEX27;
    }
    else{
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported built-in meshtype !!!                          ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***        currently, AsFem support:                                  ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***                              1D: edge2, edge3, edge4              ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***                              2D: quad4, quad8, quad9              ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***                              3D: hex8, hex20, hex27               ***\n");
        Msg_AsFem_Exit();
        abort();
    }
}