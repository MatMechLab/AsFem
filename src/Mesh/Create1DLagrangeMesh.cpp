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
//+++ Purpose: Implement the mesh generation for 1D lagrange mesh
//+++          the related mesh type is: edge2,edge3,edge4...
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Mesh/LagrangeMesh.h"

bool LagrangeMesh::Create1DLagrangeMesh(){
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
        MessagePrinter::PrintErrorTxt("unsupported mesh type for 1D case. Currently, only edge2,edge3,edge4 and edge5 are supported.");
        _IsMeshCreated=false;
        return _IsMeshCreated;
    }
    
    _nBulkElmts=_Nx;
    _nElmts=_nBulkElmts+2;
    _nNodes=_nBulkElmts*_nOrder+1;

    double dx=(_Xmax-_Xmin)/(_nNodes-1);
    double dy=(_Ymax-_Ymin)/(_nNodes-1);
    double dz=(_Zmax-_Zmin)/(_nNodes-1);

    _TotalVolume=_Xmax-_Xmin;
    _ElmtVolume.resize(_nElmts,_TotalVolume/_nBulkElmts);
    _ElmtVolume[0]=0.0;// the 'volume' of point is zero
    _ElmtVolume[1]=0.0;// the 'volume' of point is zero

    _ElmtConn.resize(_nElmts,vector<PetscInt>(0));
    _ElmtVTKCellTypeList.resize(_nElmts,0);
    _NodeCoords.resize(_nNodes*3,0.0);
    for(PetscInt i=0;i<_nNodes;++i){
        _NodeCoords[i*3+1-1]=_Xmin+i*dx;
        _NodeCoords[i*3+2-1]=_Ymin+i*dy;
        _NodeCoords[i*3+3-1]=_Zmin+i*dz;
    }

    // generate the element connectivity information
    _ElmtConn[0].push_back(1);// first component is length
    _ElmtConn[0].push_back(1);// second one is the node index(start from 1, not zero!!!)

    _ElmtConn[1].push_back(1);
    _ElmtConn[1].push_back(_nNodes);

    vector<int> left,right;
    left.clear();left.push_back(1);// element id index
    right.clear();right.push_back(2);//element id index, "right" is the second element in global array

    _PhysicalName2ElmtIDsList.clear();
    _PhysicalName2ElmtIDsList.push_back(make_pair("left",left));
    _PhysicalName2ElmtIDsList.push_back(make_pair("right",right));

    vector<int> tempconn;
    tempconn.clear();
    for(int e=0;e<_nBulkElmts;++e){
        _ElmtConn[2+e].clear();
        _ElmtConn[2+e].push_back(_nNodesPerBulkElmt);
        tempconn.push_back(2+e+1);
        for(int j=1;j<=_nNodesPerBulkElmt;++j){
            _ElmtConn[2+e].push_back(e*_nOrder+j);
        }
        _ElmtVTKCellTypeList[2+e]=_BulkElmtVTKCellType;
    }
    _PhysicalName2ElmtIDsList.push_back(make_pair("alldomain",tempconn));

    //*** generate the nodeset physical information
    vector<int> leftnodeids,rightnodeids,allnodeids;
    leftnodeids.clear();leftnodeids.push_back(1);
    rightnodeids.clear();rightnodeids.push_back(_nNodes);
    allnodeids.resize(_nNodes,0);
    iota(allnodeids.begin(),allnodeids.end(),1);
    
    _PhysicalName2NodeIDsList.clear();
    _PhysicalName2NodeIDsList.push_back(make_pair("left",leftnodeids));
    _PhysicalName2NodeIDsList.push_back(make_pair("right",rightnodeids));
    _PhysicalName2NodeIDsList.push_back(make_pair("alldomain",allnodeids));

    _nPhysicalGroups=3;
    _PhysicalGroupNameList.clear();
    _PhysicalGroupNameList.push_back("left");
    _PhysicalGroupNameList.push_back("right");
    _PhysicalGroupNameList.push_back("alldomain");

    _PhysicalGroupDimList.clear();
    _PhysicalGroupDimList.push_back(0);
    _PhysicalGroupDimList.push_back(0);
    _PhysicalGroupDimList.push_back(1);

    _PhysicalGroupIDList.clear();
    _PhysicalGroupIDList.push_back(1);
    _PhysicalGroupIDList.push_back(2);
    _PhysicalGroupIDList.push_back(3);

    _PhysicalGroupName2IDList.clear();
    _PhysicalGroupName2IDList.push_back(make_pair("left",     1));
    _PhysicalGroupName2IDList.push_back(make_pair("right",    2));
    _PhysicalGroupName2IDList.push_back(make_pair("alldomain",3));

    _PhysicalGroupID2NameList.clear();
    _PhysicalGroupID2NameList.push_back(make_pair(1,"left"));
    _PhysicalGroupID2NameList.push_back(make_pair(2,"right"));
    _PhysicalGroupID2NameList.push_back(make_pair(3,"alldomain"));


    _PhysicalGroupName2DimList.push_back(make_pair("left",     0));
    _PhysicalGroupName2DimList.push_back(make_pair("right",    0));
    _PhysicalGroupName2DimList.push_back(make_pair("alldomain",1));

    _PhysicalGroupName2NodesNumPerElmtList.clear();
    _PhysicalGroupName2NodesNumPerElmtList.push_back(make_pair("left", 1));
    _PhysicalGroupName2NodesNumPerElmtList.push_back(make_pair("right",1));
    _PhysicalGroupName2NodesNumPerElmtList.push_back(make_pair("alldomain",_nNodesPerBulkElmt));

    _IsMeshCreated=true;
    return _IsMeshCreated;
}