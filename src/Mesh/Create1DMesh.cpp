//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "Mesh/Mesh.h"

bool Mesh::Create1DMesh(){
    _IsMeshCreated=false;

    _nMaxDim=1;_nMinDim=0;
    if(_BulkMeshType==MeshType::EDGE2){
        _LineMeshType=MeshType::EDGE2;
        _SurfaceMeshType=MeshType::NULLTYPE;
        _BulkElmtVTKCellType=3;
        _nNodesPerBulkElmt=2;
        _nNodesPerLineElmt=0;
        _nNodesPerSurfaceElmt=0;
        _nOrder=1;
        _BulkMeshTypeName="edge2";
    }
    else if(_BulkMeshType==MeshType::EDGE3){
        _LineMeshType=MeshType::EDGE3;
        _SurfaceMeshType=MeshType::NULLTYPE;
        _BulkElmtVTKCellType=4;
        _nNodesPerBulkElmt=3;
        _nNodesPerLineElmt=0;
        _nNodesPerSurfaceElmt=0;
        _nOrder=2;
        _BulkMeshTypeName="edge3";
    }
    else if(_BulkMeshType==MeshType::EDGE4){
        _LineMeshType=MeshType::EDGE4;
        _SurfaceMeshType=MeshType::NULLTYPE;
        _BulkElmtVTKCellType=4;
        _nNodesPerBulkElmt=4;
        _nNodesPerLineElmt=0;
        _nNodesPerSurfaceElmt=0;
        _nOrder=3;
        _BulkMeshTypeName="edge4";
    }
    else if(_BulkMeshType==MeshType::EDGE5){
        _LineMeshType=MeshType::EDGE5;
        _SurfaceMeshType=MeshType::NULLTYPE;
        _BulkElmtVTKCellType=4;
        _nNodesPerBulkElmt=5;
        _nNodesPerLineElmt=0;
        _nNodesPerSurfaceElmt=0;
        _nOrder=4;
        _BulkMeshTypeName="edge5";
    }
    else{
        PetscPrintf(PETSC_COMM_WORLD,"*************************************************************************\n");
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported mesh type for 1D case !!!                      ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***        currently, only edge2,edge3,edge4 and edge5 are supported!!***\n");
        PetscPrintf(PETSC_COMM_WORLD,"*************************************************************************\n");
        _IsMeshCreated=false;
        return _IsMeshCreated;
    }

    _nBulkElmts=_Nx;
    _nElmts=_nBulkElmts+2;
    _nNodes=_nBulkElmts*_nOrder+1;

    PetscReal dx=(_Xmax-_Xmin)/(_nNodes-1);

    _TotalVolume=_Xmax-_Xmin;
    _ElmtVolume.resize(_nElmts,_TotalVolume/_nBulkElmts);
    _ElmtVolume[0]=0.0;// the 'volume' of point is zero
    _ElmtVolume[1]=0.0;// the 'volume' of point is zero

    _ElmtConn.resize(_nElmts,vector<PetscInt>(0));
    _ElmtVTKCellType.resize(_nElmts,0);
    _NodeCoords.resize(_nNodes*4,0.0);
    for(PetscInt i=0;i<_nNodes;++i){
        _NodeCoords[i*4+0]=1.0;
        _NodeCoords[i*4+1]=_Xmin+i*dx;
        _NodeCoords[i*4+2]=0.0;
        _NodeCoords[i*4+3]=0.0;
    }

    // generate the element connectivity information
    _ElmtConn[0].push_back(1);// first component is length
    _ElmtConn[0].push_back(1);// second one is the node index(start from 1, not zero!!!)

    _ElmtConn[1].push_back(1);
    _ElmtConn[1].push_back(_nNodes);

    vector<PetscInt> left,right;
    left.clear();left.push_back(1);// element id index
    right.clear();right.push_back(2);//element id index, "right" is the second element in global array

    _PhysicNameToElmtIndexSet.clear();
    _PhysicNameToElmtIndexSet.push_back(make_pair("left",left));
    _PhysicNameToElmtIndexSet.push_back(make_pair("right",right));

    vector<PetscInt> tempconn;
    tempconn.clear();
    for(PetscInt e=0;e<_nBulkElmts;++e){
        _ElmtConn[2+e].clear();
        _ElmtConn[2+e].push_back(_nNodesPerBulkElmt);
        tempconn.push_back(2+e+1);
        for(PetscInt j=1;j<=_nNodesPerBulkElmt;++j){
            _ElmtConn[2+e].push_back(e*_nOrder+j);
        }
        _ElmtVTKCellType[2+e]=_BulkElmtVTKCellType;
    }
    _PhysicNameToElmtIndexSet.push_back(make_pair("alldomain",tempconn));


    _PhysicGroupNameList.clear();
    _PhysicGroupNameList.push_back("left");
    _PhysicGroupNameList.push_back("right");
    _PhysicGroupNameList.push_back("alldomain");

    _PhysicNameToIDList.clear();
    _PhysicNameToIDList.push_back(make_pair("left",1));
    _PhysicNameToIDList.push_back(make_pair("right",2));
    _PhysicNameToIDList.push_back(make_pair("alldomain",3));

    _PhysicIDToNameList.clear();
    _PhysicIDToNameList.push_back(make_pair(1,"left"));
    _PhysicIDToNameList.push_back(make_pair(2,"right"));
    _PhysicIDToNameList.push_back(make_pair(3,"alldomain"));

    _nPhysicGroups=3;

    _PhysicNameToDimList.clear();
    _PhysicNameToDimList.push_back(make_pair("left",0));
    _PhysicNameToDimList.push_back(make_pair("right",0));
    _PhysicNameToDimList.push_back(make_pair("alldomain",1));

    _PhysicNameToElmtNodesNumList.clear();
    _PhysicNameToElmtNodesNumList.push_back(make_pair("left",1));
    _PhysicNameToElmtNodesNumList.push_back(make_pair("right",1));
    _PhysicNameToElmtNodesNumList.push_back(make_pair("alldomain",_nNodesPerBulkElmt));

    _IsMeshCreated=true;
    return _IsMeshCreated;
}