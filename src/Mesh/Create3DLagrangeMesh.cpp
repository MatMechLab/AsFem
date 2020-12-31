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
//+++ Purpose: Implement the mesh generation for 3D lagrange mesh
//+++          the related mesh type is: hex8,hex20,hex27
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Mesh/LagrangeMesh.h"

bool LagrangeMesh::Create3DLagrangeMesh(){
    _nNodes=0;_nElmts=0;_nBulkElmts=0;
    _nMinDim=2;_nMaxDim=3;
    _nPhysicalGroups=0;
    _NodeCoords.clear();
    _ElmtConn.clear();
    _PhysicalGroupNameList.clear();
    _PhysicalGroupDimList.clear();
    _PhysicalGroupIDList.clear();

    _PhysicalGroupID2NameList.clear();
    _PhysicalGroupName2IDList.clear();
    _PhysicalGroupName2DimList.clear();

    _PhysicalGroupName2NodesNumPerElmtList.clear();
    _PhysicalName2ElmtIDsList.clear();
    _PhysicalName2NodeIDsList.clear();


    _nNodesPerBulkElmt=0;_nNodesPerSurfaceElmt=0;_nNodesPerLineElmt=0;

    vector<int> tempconn;

    _IsMeshCreated=false;

    //******************************************
    //*** Start to generate mesh
    //******************************************
    double dx,dy,dz;
    int e,i,j,k,kk;
    int i1,i2,i3,i4,i5,i6,i7,i8,i9;
    int i10,i11,i12,i13,i14,i15,i16,i17,i18,i19;
    int i20,i21,i22,i23,i24,i25,i26,i27;
    int nNodesPerBCElmt=4;
    int VTKCellType;

    if(_BulkMeshType==MeshType::HEX8){
        _nOrder=1;
        dx=(_Xmax-_Xmin)/_Nx;
        dy=(_Ymax-_Ymin)/_Ny;
        dz=(_Zmax-_Zmin)/_Nz;

        _nBulkElmts=_Nx*_Ny*_Nz;
        _nElmts=_nBulkElmts+2*(_Nx*_Ny+_Nx*_Nz+_Ny*_Nz);
        _nSurfaceElmts=2*(_Nx*_Ny+_Nx*_Nz+_Ny*_Nz);
        _nLineElmts=0;
        _nNodes=(_Nx+1)*(_Ny+1)*(_Nz+1);
        _nNodesPerBulkElmt=8;
        nNodesPerBCElmt=4;

        _BulkMeshType=MeshType::HEX8;
        _SurfaceMeshType=MeshType::QUAD4;
        _LineMeshType=MeshType::EDGE2;

        _nNodesPerBulkElmt=8;
        _nNodesPerSurfaceElmt=4;
        _nNodesPerLineElmt=2;

        VTKCellType=12;

        _NodeCoords.resize(_nNodes*3,0.0);
        for(k=1;k<=_Nz+1;k++){
            for(j=1;j<=_Ny+1;++j){
                for(i=1;i<=_Nx+1;++i){
                    kk=(j-1)*(_Nx+1)+i+(k-1)*(_Nx+1)*(_Ny+1);
                    _NodeCoords[(kk-1)*3+1-1]=_Xmin+(i-1)*dx;
                    _NodeCoords[(kk-1)*3+2-1]=_Ymin+(j-1)*dy;
                    _NodeCoords[(kk-1)*3+3-1]=_Zmin+(k-1)*dz;
                }
            }
        }
        // Create Connectivity matrix
        _ElmtConn.resize(_nElmts,vector<int>(0));
        _ElmtVolume.resize(_nElmts,0.0);
        tempconn.clear();
        _ElmtVTKCellTypeList.resize(_nElmts,0);
        kk=0;
        for(k=1;k<=_Nz;k++){
            for(j=1;j<=_Ny;j++){
                for(i=1;i<=_Nx;i++){
                    e=(j-1)*_Nx+i+(k-1)*_Nx*_Ny;
                    i1=(j-1)*(_Nx+1)+i+(k-1)*(_Nx+1)*(_Ny+1);
                    i2=i1+1;
                    i3=i2+_Nx+1;
                    i4=i3-1;
                    i5=i1+(_Nx+1)*(_Ny+1);
                    i6=i2+(_Nx+1)*(_Ny+1);
                    i7=i3+(_Nx+1)*(_Ny+1);
                    i8=i4+(_Nx+1)*(_Ny+1);

                    _ElmtConn[e-1+_nElmts-_nBulkElmts].clear();
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(_nNodesPerBulkElmt);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i1);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i2);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i3);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i4);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i5);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i6);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i7);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i8);
                    _ElmtVTKCellTypeList[e-1+_nElmts-_nBulkElmts]=VTKCellType;

                    _ElmtVolume[e-1+_nElmts-_nBulkElmts]=dx*dy*dz;

                    tempconn.push_back(e+_nElmts-_nBulkElmts);
                }
            }
        }
    }
    else if(_BulkMeshType==MeshType::HEX20){
        _nOrder=2;
        // for 2D-8 nodes mesh
        dx=(_Xmax-_Xmin)/(2.0*_Nx);
        dy=(_Ymax-_Ymin)/(2.0*_Ny);
        dz=(_Zmax-_Zmin)/(2.0*_Nz);

        _nBulkElmts=_Nx*_Ny*_Nz;
        _nElmts=_nBulkElmts+2*(_Nx*_Ny+_Nx*_Nz+_Ny*_Nz);
        _nSurfaceElmts=2*(_Nx*_Ny+_Nx*_Nz+_Ny*_Nz);
        _nLineElmts=0;
        int nLayer1Nodes=(2*_Nx+1)*(2*_Ny+1)-_Nx*_Ny;// for norm layer
        int nLayer2Nodes=(_Nx+1)*(_Ny+1);          // for middle layer

        _nNodes=nLayer1Nodes*(_Nz+1)+nLayer2Nodes*_Nz;
        _nNodesPerBulkElmt=20;
        nNodesPerBCElmt=8;

        _BulkMeshType=MeshType::HEX20;
        _SurfaceMeshType=MeshType::QUAD8;
        _LineMeshType=MeshType::EDGE3;

        _nNodesPerBulkElmt=20;
        _nNodesPerSurfaceElmt=8;
        _nNodesPerLineElmt=3;

        VTKCellType=25;

        _NodeCoords.resize(_nNodes*3,0.0);
        for(k=1;k<=_Nz;++k){
            // First for normal layer
            for(j=1;j<=_Ny;++j){
                // for bottom line of each element
                for(i=1;i<=2*_Nx+1;i++){
                    kk=(j-1)*(2*_Nx+1+_Nx+1)+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                    _NodeCoords[(kk-1)*3+1-1]=_Xmin+(i-1)*dx;
                    _NodeCoords[(kk-1)*3+2-1]=_Ymin+(j-1)*2*dy;
                    _NodeCoords[(kk-1)*3+3-1]=_Zmin+(k-1)*2*dz;
                }
                // for middle line of each element
                for(i=1;i<=_Nx+1;++i){
                    kk=(j-1)*(2*_Nx+1+_Nx+1)+2*_Nx+1+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                    _NodeCoords[(kk-1)*3+1-1]=_Xmin+(i-1)*2*dx;
                    _NodeCoords[(kk-1)*3+2-1]=_Ymin+(j-1)*2*dy+dy;
                    _NodeCoords[(kk-1)*3+3-1]=_Zmin+(k-1)*2*dz;
                }
            }
            // for top line
            j=_Ny+1;
            for(i=1;i<=2*_Nx+1;i++){
                kk=(j-1)*(2*_Nx+1+_Nx+1)+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                _NodeCoords[(kk-1)*3+1-1]=_Xmin+(i-1)*dx;
                _NodeCoords[(kk-1)*3+2-1]=_Ymin+(j-1)*2*dy;
                _NodeCoords[(kk-1)*3+3-1]=_Zmin+(k-1)*2*dz;
            }
            // Then for middle type layer
            for(j=1;j<=_Ny+1;++j){
                for(i=1;i<=_Nx+1;++i){
                    kk=(j-1)*(_Nx+1)+i+(k-1)*(nLayer1Nodes+nLayer2Nodes)+nLayer1Nodes;
                    _NodeCoords[(kk-1)*3+1-1]=_Xmin+(i-1)*2*dx;
                    _NodeCoords[(kk-1)*3+2-1]=_Ymin+(j-1)*2*dy;
                    _NodeCoords[(kk-1)*3+3-1]=_Zmin+(k-1)*2*dz+dz;
                }
            }
        }
        // for the last top layer
        k=_Nz+1;
        for(j=1;j<=_Ny;++j){
            // for bottom line of each element
            for(i=1;i<=2*_Nx+1;i++){
                kk=(j-1)*(2*_Nx+1+_Nx+1)+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                _NodeCoords[(kk-1)*3+1-1]=_Xmin+(i-1)*dx;
                _NodeCoords[(kk-1)*3+2-1]=_Ymin+(j-1)*2*dy;
                _NodeCoords[(kk-1)*3+3-1]=_Zmin+(k-1)*2*dz;
            }
            // for middle line of each element
            for(i=1;i<=_Nx+1;++i){
                kk=(j-1)*(2*_Nx+1+_Nx+1)+2*_Nx+1+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                _NodeCoords[(kk-1)*3+1-1]=_Xmin+(i-1)*2*dx;
                _NodeCoords[(kk-1)*3+2-1]=_Ymin+(j-1)*2*dy+dy;
                _NodeCoords[(kk-1)*3+3-1]=_Zmin+(k-1)*2*dz;
            }
        }
        // for top line
        j=_Ny+1;
        for(i=1;i<=2*_Nx+1;i++){
            kk=(j-1)*(2*_Nx+1+_Nx+1)+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
            _NodeCoords[(kk-1)*3+1-1]=_Xmin+(i-1)*dx;
            _NodeCoords[(kk-1)*3+2-1]=_Ymin+(j-1)*2*dy;
            _NodeCoords[(kk-1)*3+3-1]=_Zmin+(k-1)*2*dz;
        }
        //***************************************
        // Create Connectivity matrix
        //***************************************
        tempconn.clear();
        _ElmtConn.resize(_nElmts,vector<int>(0));
        _ElmtVolume.resize(_nElmts,0.0);
        _ElmtVTKCellTypeList.resize(_nElmts,0);
        for(k=1;k<=_Nz;++k){
            for(j=1;j<=_Ny;++j){
                for(i=1;i<=_Nx;++i){
                    e=(j-1)*_Nx+i+(k-1)*_Nx*_Ny;
                    i1=(j-1)*(2*_Nx+1+_Nx+1)+2*i-1+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                    i2=i1+2;
                    i3=i2+(2*_Nx+1+_Nx+1);
                    i4=i3-2;

                    i5=i1+nLayer1Nodes+nLayer2Nodes;
                    i6=i2+nLayer1Nodes+nLayer2Nodes;
                    i7=i3+nLayer1Nodes+nLayer2Nodes;
                    i8=i4+nLayer1Nodes+nLayer2Nodes;

                    i9 =i1+1;
                    i10=i2+(2*_Nx+1-i);
                    i11=i3-1;
                    i12=i10-1;

                    i13= i9+nLayer1Nodes+nLayer2Nodes;
                    i14=i10+nLayer1Nodes+nLayer2Nodes;
                    i15=i11+nLayer1Nodes+nLayer2Nodes;
                    i16=i12+nLayer1Nodes+nLayer2Nodes;

                    i17=i1+nLayer1Nodes-(i-1+(j-1)*(_Nx+_Nx+1));
                    i18=i17+1;
                    i19=i18+_Nx+1;
                    i20=i19-1;

                    _ElmtConn[e-1+_nElmts-_nBulkElmts].clear();
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(_nNodesPerBulkElmt);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i1);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i2);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i3);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i4);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i5);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i6);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i7);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i8);

                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i9);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i10);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i11);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i12);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i13);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i14);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i15);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i16);

                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i17);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i18);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i19);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i20);

                    _ElmtVTKCellTypeList[e-1+_nElmts-_nBulkElmts]=VTKCellType;

                    _ElmtVolume[e-1+_nElmts-_nBulkElmts]=2*dx*2*dy*2*dz;

                    tempconn.push_back(e+_nElmts-_nBulkElmts);

                }
            }
        }
    }
    else if(_BulkMeshType==MeshType::HEX27){
        _nOrder=2;
        // for 3D-27 nodes mesh
        dx=(_Xmax-_Xmin)/(2.0*_Nx);
        dy=(_Ymax-_Ymin)/(2.0*_Ny);
        dz=(_Zmax-_Zmin)/(2.0*_Nz);


        _nBulkElmts=_Nx*_Ny*_Nz;
        _nElmts=_nBulkElmts+2*(_Nx*_Ny+_Nx*_Nz+_Ny*_Nz);
        _nSurfaceElmts=2*(_Nx*_Ny+_Nx*_Nz+_Ny*_Nz);
        _nLineElmts=0;
        int nLayerNodes=(2*_Nx+1)*(2*_Ny+1);// for norm layer

        _nNodes=(2*_Nz+1)*nLayerNodes;
        _nNodesPerBulkElmt=27;
        nNodesPerBCElmt=9;


        _BulkMeshType=MeshType::HEX27;
        _SurfaceMeshType=MeshType::QUAD9;
        _LineMeshType=MeshType::EDGE3;

        _nNodesPerBulkElmt=27;
        _nNodesPerSurfaceElmt=9;
        _nNodesPerLineElmt=3;

        VTKCellType=29;

        _NodeCoords.resize(3*_nNodes,0.0);
        for(k=1;k<=_Nz;++k){
            // For first layer
            for(j=1;j<=2*_Ny+1;++j){
                // for bottom line of each element
                for(i=1;i<=2*_Nx+1;++i){
                    kk=(j-1)*(2*_Nx+1)+i+(k-1)*2*nLayerNodes;
                    _NodeCoords[(kk-1)*3+1-1]=_Xmin+(i-1)*dx;
                    _NodeCoords[(kk-1)*3+2-1]=_Ymin+(j-1)*dy;
                    _NodeCoords[(kk-1)*3+3-1]=_Zmin+(k-1)*2*dz;
                }
            }
            // Then for second layer
            for(j=1;j<=2*_Ny+1;++j){
                for(i=1;i<=2*_Nx+1;++i){
                    kk=(j-1)*(2*_Nx+1)+i+(k-1)*2*nLayerNodes+nLayerNodes;
                    _NodeCoords[(kk-1)*3+1-1]=_Xmin+(i-1)*dx;
                    _NodeCoords[(kk-1)*3+2-1]=_Ymin+(j-1)*dy;
                    _NodeCoords[(kk-1)*3+3-1]=_Zmin+(k-1)*2*dz+dz;
                }
            }
        }
        // for the last top layer
        k=_Nz+1;
        for(j=1;j<=2*_Ny+1;++j){
            // for bottom line of each element
            for(i=1;i<=2*_Nx+1;++i){
                kk=(j-1)*(2*_Nx+1)+i+(k-1)*2*nLayerNodes;
                _NodeCoords[(kk-1)*3+1-1]=_Xmin+(i-1)*dx;
                _NodeCoords[(kk-1)*3+2-1]=_Ymin+(j-1)*dy;
                _NodeCoords[(kk-1)*3+3-1]=_Zmin+(k-1)*2*dz;
            }
        }
        // Create Connectivity matrix
        tempconn.clear();
        _ElmtConn.resize(_nElmts,vector<int>(0));
        _ElmtVolume.resize(_nElmts,0.0);
        _ElmtVTKCellTypeList.resize(_nElmts,0);
        for(k=1;k<=_Nz;++k){
            for(j=1;j<=_Ny;++j){
                for(i=1;i<=_Nx;++i){
                    e=(j-1)*_Nx+i+(k-1)*_Nx*_Ny;
                    i1=(j-1)*2*(2*_Nx+1)+2*i-1+(k-1)*2*nLayerNodes;
                    i2=i1+2;
                    i3=i2+(2*_Nx+1)*2;
                    i4=i3-2;

                    i5=i1+2*nLayerNodes;
                    i6=i2+2*nLayerNodes;
                    i7=i3+2*nLayerNodes;
                    i8=i4+2*nLayerNodes;

                    i9 =i1+1;
                    i10=i2+(2*_Nx+1);
                    i11=i3-1;
                    i12=i1+(2*_Nx+1);

                    i13=i5+1;
                    i14=i6+(2*_Nx+1);
                    i15=i7-1;
                    i16=i5+(2*_Nx+1);

                    i17=i1+nLayerNodes;
                    i18=i2+nLayerNodes;
                    i19=i3+nLayerNodes;
                    i20=i4+nLayerNodes;

                    i21=i17+(2*_Nx+1);
                    i22=i21+2;

                    //i23=i20+1;
                    //i24=i17+1;

                    i23=i17+1;
                    i24=i20+1;

                    i25=i12+1;
                    i26=i16+1;

                    i27=i21+1;

                    _ElmtConn[e-1+_nElmts-_nBulkElmts].clear();
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(_nNodesPerBulkElmt);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i1);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i2);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i3);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i4);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i5);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i6);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i7);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i8);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i9);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i10);

                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i11);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i12);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i13);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i14);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i15);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i16);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i17);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i18);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i19);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i20);

                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i21);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i22);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i23);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i24);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i25);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i26);
                    _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i27);

                    _ElmtVTKCellTypeList[e-1+_nElmts-_nBulkElmts]=VTKCellType;

                    _ElmtVolume[e-1+_nElmts-_nBulkElmts]=2*dx*2*dy*2*dz;

                    tempconn.push_back(e+_nElmts-_nBulkElmts);
                }
            }
        }
    }

    //*********************************************
    //*** Now we split the boundary mesh
    //*********************************************
    _nPhysicalGroups=6+1;
    vector<int> left,right,bottom,top,back,front;
    vector<int> leftnodeids,rightnodeids;
    vector<int> bottomnodeids,topnodeids;
    vector<int> backnodeids,frontnodeids,allnodeids;

    _PhysicalGroupNameList.clear();
    _PhysicalGroupNameList.push_back("left");
    _PhysicalGroupNameList.push_back("right");
    _PhysicalGroupNameList.push_back("bottom");
    _PhysicalGroupNameList.push_back("top");
    _PhysicalGroupNameList.push_back("back");
    _PhysicalGroupNameList.push_back("front");
    _PhysicalGroupNameList.push_back("alldomain");

    _PhysicalGroupDimList.clear();
    _PhysicalGroupDimList.push_back(2);
    _PhysicalGroupDimList.push_back(2);
    _PhysicalGroupDimList.push_back(2);
    _PhysicalGroupDimList.push_back(2);
    _PhysicalGroupDimList.push_back(2);
    _PhysicalGroupDimList.push_back(2);
    _PhysicalGroupDimList.push_back(3);

    _PhysicalGroupIDList.clear();
    _PhysicalGroupIDList.push_back(1);
    _PhysicalGroupIDList.push_back(2);
    _PhysicalGroupIDList.push_back(3);
    _PhysicalGroupIDList.push_back(4);
    _PhysicalGroupIDList.push_back(5);
    _PhysicalGroupIDList.push_back(6);
    _PhysicalGroupIDList.push_back(7);

    //******************************
    _PhysicalGroupName2IDList.clear();
    _PhysicalGroupName2IDList.push_back(make_pair("left",     1));
    _PhysicalGroupName2IDList.push_back(make_pair("right",    2));
    _PhysicalGroupName2IDList.push_back(make_pair("bottom",   3));
    _PhysicalGroupName2IDList.push_back(make_pair("top",      4));
    _PhysicalGroupName2IDList.push_back(make_pair("back",     5));
    _PhysicalGroupName2IDList.push_back(make_pair("front",    6));
    _PhysicalGroupName2IDList.push_back(make_pair("alldomain",7));

    _PhysicalGroupID2NameList.clear();
    _PhysicalGroupID2NameList.push_back(make_pair(1,"left"));
    _PhysicalGroupID2NameList.push_back(make_pair(2,"right"));
    _PhysicalGroupID2NameList.push_back(make_pair(3,"bottom"));
    _PhysicalGroupID2NameList.push_back(make_pair(4,"top"));
    _PhysicalGroupID2NameList.push_back(make_pair(5,"back"));
    _PhysicalGroupID2NameList.push_back(make_pair(6,"front"));
    _PhysicalGroupID2NameList.push_back(make_pair(7,"alldomain"));

    _PhysicalGroupName2DimList.clear();
    _PhysicalGroupName2DimList.push_back(make_pair("left",     2));
    _PhysicalGroupName2DimList.push_back(make_pair("right",    2));
    _PhysicalGroupName2DimList.push_back(make_pair("bottom",   2));
    _PhysicalGroupName2DimList.push_back(make_pair("top",      2));
    _PhysicalGroupName2DimList.push_back(make_pair("back",     2));
    _PhysicalGroupName2DimList.push_back(make_pair("front",    2));
    _PhysicalGroupName2DimList.push_back(make_pair("alldomain",3));

    
    //***********************************************
    int nBCElmts=0;
    // for leftside
    left.clear();
    leftnodeids.clear();
    i=1;
    for(k=1;k<=_Nz;k++){
        for(j=1;j<=_Ny;j++){
            e=(j-1)*_Nx+i+(k-1)*_Nx*_Ny;
            _ElmtConn[nBCElmts].clear();
            if(nNodesPerBCElmt==4){
                // hex8 case
                // must out of plane
                _ElmtConn[nBCElmts].push_back(4);
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,1));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,5));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,8));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,4));
                _ElmtVTKCellTypeList[nBCElmts]=9;

                leftnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,1));
                leftnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,5));
                leftnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,8));
                leftnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,4));
            }
            else if(nNodesPerBCElmt==8){
                _ElmtConn[nBCElmts].push_back(8);
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,1));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,5));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,8));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,4));

                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,17));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,16));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,20));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,12));

                _ElmtVTKCellTypeList[nBCElmts]=23;

                leftnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,1));
                leftnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,5));
                leftnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,8));
                leftnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,4));

                leftnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,17));
                leftnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,16));
                leftnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,20));
                leftnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,12));
            }
            else if(nNodesPerBCElmt==9){
                _ElmtConn[nBCElmts].push_back(9);
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,1));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,5));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,8));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,4));

                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,17));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,16));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,20));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,12));

                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,21));

                _ElmtVTKCellTypeList[nBCElmts]=28;

                leftnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,1));
                leftnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,5));
                leftnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,8));
                leftnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,4));

                leftnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,17));
                leftnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,16));
                leftnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,20));
                leftnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,12));

                leftnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,21));
            }
            nBCElmts+=1;
            left.push_back(nBCElmts);
        }
    }
    sort(leftnodeids.begin(),leftnodeids.end());
    leftnodeids.erase( unique(leftnodeids.begin(),leftnodeids.end()), leftnodeids.end());

    // For right side
    right.clear();
    rightnodeids.clear();
    i=_Nx;
    for(k=1;k<=_Nz;k++){
        for(j=1;j<=_Ny;j++){
            e=(j-1)*_Nx+i+(k-1)*_Nx*_Ny;
            _ElmtConn[nBCElmts].clear();
            if(nNodesPerBCElmt==4){
                // hex8 case
                // must out of plane
                _ElmtConn[nBCElmts].push_back(4);
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,2));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,3));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,7));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,6));
                _ElmtVTKCellTypeList[nBCElmts]=9;
                
                rightnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,2));
                rightnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,3));
                rightnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,7));
                rightnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,6));
            }
            else if(nNodesPerBCElmt==8){
                _ElmtConn[nBCElmts].push_back(8);
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,2));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,3));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,7));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,6));

                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,10));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,19));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,14));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,18));
                _ElmtVTKCellTypeList[nBCElmts]=23;

                rightnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,2));
                rightnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,3));
                rightnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,7));
                rightnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,6));

                rightnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,10));
                rightnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,19));
                rightnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,14));
                rightnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,18));
            }
            else if(nNodesPerBCElmt==9){
                _ElmtConn[nBCElmts].push_back(9);
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,2));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,3));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,7));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,6));

                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,10));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,19));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,14));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,18));

                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,22));
                _ElmtVTKCellTypeList[nBCElmts]=28;

                rightnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,2));
                rightnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,3));
                rightnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,7));
                rightnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,6));

                rightnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,10));
                rightnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,19));
                rightnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,14));
                rightnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,18));

                rightnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,22));
            }
            nBCElmts+=1;
            right.push_back(nBCElmts);
        }
    }
    sort(rightnodeids.begin(),rightnodeids.end());
    rightnodeids.erase( unique(rightnodeids.begin(),rightnodeids.end()), rightnodeids.end());

    // For bottom edge
    bottom.clear();
    bottomnodeids.clear();
    j=1;
    for(k=1;k<=_Nz;k++){
        for(i=1;i<=_Nx;i++){
            e=(j-1)*_Nx+i+(k-1)*_Nx*_Ny;
            _ElmtConn[nBCElmts].clear();
            if(nNodesPerBCElmt==4){
                // hex8 case
                // must out of plane
                _ElmtConn[nBCElmts].push_back(4);
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,1));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,2));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,6));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,5));
                _ElmtVTKCellTypeList[nBCElmts]=9;

                bottomnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,1));
                bottomnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,2));
                bottomnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,6));
                bottomnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,5));
            }
            else if(nNodesPerBCElmt==8){
                _ElmtConn[nBCElmts].push_back(8);
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,1));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,2));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,6));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,5));

                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,9));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,18));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,13));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,17));
                _ElmtVTKCellTypeList[nBCElmts]=23;

                bottomnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,1));
                bottomnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,2));
                bottomnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,6));
                bottomnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,5));

                bottomnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,9));
                bottomnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,18));
                bottomnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,13));
                bottomnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,17));
            }
            else if(nNodesPerBCElmt==9){
                _ElmtConn[nBCElmts].push_back(9);
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,1));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,2));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,6));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,5));

                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,9));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,18));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,13));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,17));

                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,23));
                _ElmtVTKCellTypeList[nBCElmts]=28;

                bottomnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,1));
                bottomnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,2));
                bottomnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,6));
                bottomnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,5));

                bottomnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,9));
                bottomnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,18));
                bottomnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,13));
                bottomnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,17));

                bottomnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,23));
            }
            nBCElmts+=1;
            bottom.push_back(nBCElmts);
        }
    }
    sort(bottomnodeids.begin(),bottomnodeids.end());
    bottomnodeids.erase( unique(bottomnodeids.begin(),bottomnodeids.end()), bottomnodeids.end());

    // For top edge
    top.clear();
    topnodeids.clear();
    j=_Ny;
    for(k=1;k<=_Nz;k++){
        for(i=1;i<=_Nx;i++){
            e=(j-1)*_Nx+i+(k-1)*_Nx*_Ny;
            _ElmtConn[nBCElmts].clear();
            if(nNodesPerBCElmt==4){
                // hex8 case
                // must out of plane
                _ElmtConn[nBCElmts].push_back(4);
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,4));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,8));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,7));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,3));

                _ElmtVTKCellTypeList[nBCElmts]=9;

                topnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,4));
                topnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,8));
                topnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,7));
                topnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,3));
            }
            else if(nNodesPerBCElmt==8){
                _ElmtConn[nBCElmts].push_back(8);
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,4));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,8));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,7));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,3));

                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,20));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,15));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,19));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,11));

                _ElmtVTKCellTypeList[nBCElmts]=23;

                topnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,4));
                topnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,8));
                topnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,7));
                topnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,3));

                topnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,20));
                topnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,15));
                topnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,19));
                topnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,11));
            }
            else if(nNodesPerBCElmt==9){
                _ElmtConn[nBCElmts].push_back(9);
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,4));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,8));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,7));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,3));

                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,20));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,15));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,19));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,11));

                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,24));

                _ElmtVTKCellTypeList[nBCElmts]=28;

                topnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,4));
                topnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,8));
                topnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,7));
                topnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,3));

                topnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,20));
                topnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,15));
                topnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,19));
                topnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,11));

                topnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,24));
            }
            nBCElmts+=1;
            top.push_back(nBCElmts);
        }
    }
    sort(topnodeids.begin(),topnodeids.end());
    topnodeids.erase( unique(topnodeids.begin(),topnodeids.end()), topnodeids.end());

    // For back edge
    back.clear();
    backnodeids.clear();
    k=1;
    for(j=1;j<=_Ny;j++){
        for(i=1;i<=_Nx;i++){
            e=(j-1)*_Nx+i+(k-1)*_Nx*_Ny;
            _ElmtConn[nBCElmts].clear();
            if(nNodesPerBCElmt==4){
                // hex8 case
                // must out of plane
                _ElmtConn[nBCElmts].push_back(4);
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,1));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,4));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,3));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,2));

                _ElmtVTKCellTypeList[nBCElmts]=9;

                backnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,1));
                backnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,4));
                backnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,3));
                backnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,2));
            }
            else if(nNodesPerBCElmt==8){
                _ElmtConn[nBCElmts].push_back(8);
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,1));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,4));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,3));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,2));

                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,12));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,11));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,10));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,9));

                _ElmtVTKCellTypeList[nBCElmts]=23;

                backnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,1));
                backnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,4));
                backnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,3));
                backnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,2));

                backnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,12));
                backnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,11));
                backnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,10));
                backnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,9));
            }
            else if(nNodesPerBCElmt==9){
                _ElmtConn[nBCElmts].push_back(9);
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,1));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,4));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,3));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,2));

                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,12));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,11));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,10));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,9));

                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,25));

                _ElmtVTKCellTypeList[nBCElmts]=28;

                backnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,1));
                backnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,4));
                backnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,3));
                backnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,2));

                backnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,12));
                backnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,11));
                backnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,10));
                backnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,9));

                backnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,25));
            }
            nBCElmts+=1;
            back.push_back(nBCElmts);
        }
    }
    sort(backnodeids.begin(),backnodeids.end());
    backnodeids.erase( unique(backnodeids.begin(),backnodeids.end()), backnodeids.end());

    // For front edge
    front.clear();
    frontnodeids.clear();
    k=_Nz;
    for(j=1;j<=_Ny;j++){
        for(i=1;i<=_Nx;i++){
            e=(j-1)*_Nx+i+(k-1)*_Nx*_Ny;
            _ElmtConn[nBCElmts].clear();
            if(nNodesPerBCElmt==4){
                // hex8 case
                // must out of plane
                _ElmtConn[nBCElmts].push_back(4);
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,5));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,6));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,7));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,8));

                _ElmtVTKCellTypeList[nBCElmts]=9;

                frontnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,5));
                frontnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,6));
                frontnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,7));
                frontnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,8));
            }
            else if(nNodesPerBCElmt==8){
                _ElmtConn[nBCElmts].push_back(8);
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,5));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,6));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,7));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,8));

                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,13));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,14));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,15));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,16));

                _ElmtVTKCellTypeList[nBCElmts]=23;

                frontnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,5));
                frontnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,6));
                frontnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,7));
                frontnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,8));

                frontnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,13));
                frontnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,14));
                frontnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,15));
                frontnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,16));
            }
            else if(nNodesPerBCElmt==9){
                _ElmtConn[nBCElmts].push_back(9);
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,5));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,6));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,7));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,8));

                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,13));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,14));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,15));
                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,16));

                _ElmtConn[nBCElmts].push_back(GetBulkMeshIthBulkElmtJthNodeID(e,26));

                _ElmtVTKCellTypeList[nBCElmts]=28;

                frontnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,5));
                frontnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,6));
                frontnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,7));
                frontnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,8));

                frontnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,13));
                frontnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,14));
                frontnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,15));
                frontnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,16));

                frontnodeids.push_back(GetBulkMeshIthBulkElmtJthNodeID(e,26));
            }
            nBCElmts+=1;
            front.push_back(nBCElmts);
        }
    }
    sort(frontnodeids.begin(),frontnodeids.end());
    frontnodeids.erase( unique(frontnodeids.begin(),frontnodeids.end()), frontnodeids.end());
    // add all node ids
    allnodeids.resize(_nNodes,0);
    iota(allnodeids.begin(),allnodeids.end(),1);

    _PhysicalName2NodeIDsList.clear();
    _PhysicalName2NodeIDsList.push_back(make_pair("left",leftnodeids));
    _PhysicalName2NodeIDsList.push_back(make_pair("right",rightnodeids));
    _PhysicalName2NodeIDsList.push_back(make_pair("bottom",bottomnodeids));
    _PhysicalName2NodeIDsList.push_back(make_pair("top",topnodeids));
    _PhysicalName2NodeIDsList.push_back(make_pair("back",backnodeids));
    _PhysicalName2NodeIDsList.push_back(make_pair("front",frontnodeids));
    _PhysicalName2NodeIDsList.push_back(make_pair("alldomain",allnodeids));

    _PhysicalName2ElmtIDsList.clear();
    _PhysicalName2ElmtIDsList.push_back(make_pair("left",left));
    _PhysicalName2ElmtIDsList.push_back(make_pair("right",right));

    _PhysicalName2ElmtIDsList.push_back(make_pair("bottom",bottom));
    _PhysicalName2ElmtIDsList.push_back(make_pair("top",top));

    _PhysicalName2ElmtIDsList.push_back(make_pair("back",back));
    _PhysicalName2ElmtIDsList.push_back(make_pair("front",front));

    _PhysicalName2ElmtIDsList.push_back(make_pair("alldomain",tempconn));

    //****************************************
    

    _PhysicalGroupName2NodesNumPerElmtList.clear();
    _PhysicalGroupName2NodesNumPerElmtList.push_back(make_pair("left",_nNodesPerSurfaceElmt));
    _PhysicalGroupName2NodesNumPerElmtList.push_back(make_pair("right",_nNodesPerSurfaceElmt));
    _PhysicalGroupName2NodesNumPerElmtList.push_back(make_pair("bottom",_nNodesPerSurfaceElmt));
    _PhysicalGroupName2NodesNumPerElmtList.push_back(make_pair("top",_nNodesPerSurfaceElmt));
    _PhysicalGroupName2NodesNumPerElmtList.push_back(make_pair("back",_nNodesPerSurfaceElmt));
    _PhysicalGroupName2NodesNumPerElmtList.push_back(make_pair("front",_nNodesPerSurfaceElmt));
    _PhysicalGroupName2NodesNumPerElmtList.push_back(make_pair("alldomain",_nNodesPerBulkElmt));

    _IsMeshCreated=true;

    return _IsMeshCreated;
}