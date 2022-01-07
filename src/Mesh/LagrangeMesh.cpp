//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.06.28
//+++ Purpose: Define bulk mesh class for the bulk domain(except 
//+++          the interface domain)
//+++          This class offers the functions to get as well as
//+++          create the mesh and the related geometric info
//+++          In order to import mesh from external file, one should
//+++          call the meshio class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Mesh/LagrangeMesh.h"

LagrangeMesh::LagrangeMesh(){
    _IsMeshCreated=false;
    _nNodes=0;_nElmts=0;
    _nNodesPerBulkElmt=0;_nNodesPerSurfaceElmt=0;_nNodesPerLineElmt=0;
    _nBulkElmts=0;_nSurfaceElmts=0;_nLineElmts=0;
    _nMaxDim=0;_nMinDim=0;
    _Nx=0;_Ny=0;_Nz=0;
    _Xmax=0.0;_Xmin=0.0;_Ymax=0.0;_Ymin=0.0;_Zmin=0.0;_Zmax=0.0;
    _nOrder=1;
    _BulkMeshType=MeshType::QUAD4;
    _SurfaceMeshType=MeshType::QUAD4;
    _LineMeshType=MeshType::EDGE2;
    _NodeCoords.clear();// store all the nodes/controlpts' coordinate
    _ElmtConn.clear();  // store all the element(bulk+surface+line+node elements)
    _ElmtVolume.clear();// it could be: volume(3D), area(2D) or length(1D)

    _ElmtVTKCellTypeList.clear();
    _ElmtPhyIDList.clear();
    _ElmtDimList.clear();
    _ElmtMeshTypeList.clear();
    _BulkElmtVTKCellType=0;
    _BulkMeshTypeName="quad4";
    _TotalVolume=0.0;

    //************************************************************
    //*** for the basic physical group information
    //************************************************************
    _nPhysicalGroups=0;
    _PhysicalGroupNameList.clear();
    _PhysicalGroupIDList.clear();
    _PhysicalGroupDimList.clear();
    _PhysicalGroupName2DimList.clear();
    _PhysicalGroupID2NameList.clear();
    _PhysicalGroupName2IDList.clear();
    _PhysicalGroupName2NodesNumPerElmtList.clear();
    _PhysicalName2ElmtIDsList.clear();

    //*** for node set physical groups
    _nNodeSetPhysicalGroups=0;
    _NodeSetPhysicalGroupNameList.clear();
    _NodeSetPhysicalGroupIDList.clear();
    _NodeSetPhysicalGroupID2NameList.clear();
    _NodeSetPhysicalGroupName2IDList.clear();
    _NodeSetPhysicalName2NodeIDsList.clear();
    
}

//**********************************
void LagrangeMesh::SetBulkMeshMeshTypeName(string meshname){
    if(meshname.find("edge2")!=string::npos){
        _BulkMeshTypeName="edge2";
        _BulkMeshType=MeshType::EDGE2;
    }
    else if(meshname.find("edge3")!=string::npos){
        _BulkMeshTypeName="edge3";
        _BulkMeshType=MeshType::EDGE3;
    }
    else if(meshname.find("edge4")!=string::npos){
        _BulkMeshTypeName="edge4";
        _BulkMeshType=MeshType::EDGE4;
    }
    else if(meshname.find("edge5")!=string::npos){
        _BulkMeshTypeName="edge5";
        _BulkMeshType=MeshType::EDGE5;
    }
    else if(meshname.find("quad4")!=string::npos){
        _BulkMeshTypeName="quad4";
        _BulkMeshType=MeshType::QUAD4;
    }
    else if(meshname.find("quad8")!=string::npos){
        _BulkMeshTypeName="quad8";
        _BulkMeshType=MeshType::QUAD8;
    }
    else if(meshname.find("tri3")!=string::npos){
        _BulkMeshTypeName="tri3";
        _BulkMeshType=MeshType::TRI3;
    }
    else if(meshname.find("tri6")!=string::npos){
        _BulkMeshTypeName="tri6";
        _BulkMeshType=MeshType::TRI6;
    }
    else if(meshname.find("quad9")!=string::npos){
        _BulkMeshTypeName="quad9";
        _BulkMeshType=MeshType::QUAD9;
    }
    else if(meshname.find("hex8")!=string::npos){
        _BulkMeshTypeName="hex8";
        _BulkMeshType=MeshType::HEX8;
    }
    else if(meshname.find("hex20")!=string::npos){
        _BulkMeshTypeName="hex20";
        _BulkMeshType=MeshType::HEX20;
    }
    else if(meshname.find("hex27")!=string::npos){
        _BulkMeshTypeName="hex27";
        _BulkMeshType=MeshType::HEX27;
    }
    else if(meshname.find("tet4")!=string::npos){
        _BulkMeshTypeName="tet4";
        _BulkMeshType=MeshType::TET4;
    }
    else if(meshname.find("tet10")!=string::npos){
        _BulkMeshTypeName="tet10";
        _BulkMeshType=MeshType::TET10;
    }
    else{
        MessagePrinter::PrintErrorTxt("unsupported mesh type setting");
        MessagePrinter::AsFem_Exit();
    }
}