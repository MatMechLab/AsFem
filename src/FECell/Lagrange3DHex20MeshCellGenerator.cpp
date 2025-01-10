//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2024.01.31
//+++ Purpose: the 3D hex20 lagrange FE cell generator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FECell/Lagrange3DHex20MeshCellGenerator.h"
#include "FECell/SingleMeshCell.h"

Lagrange3DHex20MeshCellGenerator::Lagrange3DHex20MeshCellGenerator(){
    m_mesh_generated=false;
}
Lagrange3DHex20MeshCellGenerator::~Lagrange3DHex20MeshCellGenerator(){
    m_mesh_generated=false;
}
//*********************************************************
bool Lagrange3DHex20MeshCellGenerator::generateFECell(FECellData &t_celldata){
    int rank;
    t_celldata.MaxDim=3;
    t_celldata.MinDim=2;

    t_celldata.ActiveDofsNum=0;
    t_celldata.TotalDofsNum=0;
    t_celldata.MaxDofsPerNode=0;

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if(rank==0){
        vector<SingleMeshCell> leftconn,rightconn;
        vector<SingleMeshCell> bottomconn,topconn;
        vector<SingleMeshCell> backconn,frontconn;
        vector<int> leftnodes,rightnodes;
        vector<int> bottomnodes,topnodes;
        vector<int> backnodes,frontnodes;
        // only create mesh on the master rank, then distributed them into different ranks !
        double dx,dy,dz;
        int i,j,k,kk,e;
        int i1,i2,i3,i4,i5,i6,i7,i8,i9,i10;
        int i11,i12,i13,i14,i15,i16,i17,i18,i19,i20;
        
        t_celldata.MeshOrder=2;
        t_celldata.BulkMeshTypeName="hex20";

        dx=(t_celldata.Xmax-t_celldata.Xmin)/(2.0*t_celldata.Nx);
        dy=(t_celldata.Ymax-t_celldata.Ymin)/(2.0*t_celldata.Ny);
        dz=(t_celldata.Zmax-t_celldata.Zmin)/(2.0*t_celldata.Nz);
            
        t_celldata.BulkElmtsNum=t_celldata.Nx*t_celldata.Ny*t_celldata.Nz;
        t_celldata.LineElmtsNum=0;
        t_celldata.SurfElmtsNum=2*(t_celldata.Nx*t_celldata.Ny
                                  +t_celldata.Nx*t_celldata.Nz
                                  +t_celldata.Ny*t_celldata.Nz);
        t_celldata.ElmtsNum=t_celldata.BulkElmtsNum
                           +t_celldata.SurfElmtsNum
                           +t_celldata.LineElmtsNum;
            
        int nLayer1Nodes=(2*t_celldata.Nx+1)*(2*t_celldata.Ny+1)-t_celldata.Nx*t_celldata.Ny;// for norm layer
        int nLayer2Nodes=(t_celldata.Nx+1)*(t_celldata.Ny+1);          // for middle layer

        t_celldata.NodesNum=nLayer1Nodes*(t_celldata.Nz+1)+nLayer2Nodes*t_celldata.Nz;
        t_celldata.NodesNumPerBulkElmt=20;
        t_celldata.NodesNumPerSurfElmt=8;
        t_celldata.NodesNumPerLineElmt=3;

        t_celldata.BulkElmtVTKCellType=25;

        t_celldata.LineElmtVTKCellType=4;
        t_celldata.LineElmtMeshType=MeshType::EDGE3;

        t_celldata.SurfElmtVTKCellType=23;
        t_celldata.SurfElmtMeshType=MeshType::QUAD8;
            
        // for the coordinates of each node
        t_celldata.NodeCoords_Global.resize(t_celldata.NodesNum*3,0.0);
        leftnodes.clear();rightnodes.clear();
        bottomnodes.clear();rightnodes.clear();
        backnodes.clear();frontnodes.clear();
        int _Nz,_Ny,_Nx;
        _Nz=t_celldata.Nz;_Ny=t_celldata.Ny;_Nx=t_celldata.Nx;
        for(k=1;k<=_Nz;++k){
            // First for normal layer
            for(j=1;j<=_Ny;++j){
                // for bottom line of each element
                for(i=1;i<=2*_Nx+1;i++){
                    kk=(j-1)*(2*_Nx+1+_Nx+1)+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                    t_celldata.NodeCoords_Global[(kk-1)*3+1-1]=t_celldata.Xmin+(i-1)*dx;
                    t_celldata.NodeCoords_Global[(kk-1)*3+2-1]=t_celldata.Ymin+(j-1)*2*dy;
                    t_celldata.NodeCoords_Global[(kk-1)*3+3-1]=t_celldata.Zmin+(k-1)*2*dz;
                }
                // for middle line of each element
                for(i=1;i<=_Nx+1;++i){
                    kk=(j-1)*(2*_Nx+1+_Nx+1)+2*_Nx+1+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                    t_celldata.NodeCoords_Global[(kk-1)*3+1-1]=t_celldata.Xmin+(i-1)*2*dx;
                    t_celldata.NodeCoords_Global[(kk-1)*3+2-1]=t_celldata.Ymin+(j-1)*2*dy+dy;
                    t_celldata.NodeCoords_Global[(kk-1)*3+3-1]=t_celldata.Zmin+(k-1)*2*dz;
                }
            }
            // for top line
            j=_Ny+1;
            for(i=1;i<=2*_Nx+1;i++){
                kk=(j-1)*(2*_Nx+1+_Nx+1)+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                t_celldata.NodeCoords_Global[(kk-1)*3+1-1]=t_celldata.Xmin+(i-1)*dx;
                t_celldata.NodeCoords_Global[(kk-1)*3+2-1]=t_celldata.Ymin+(j-1)*2*dy;
                t_celldata.NodeCoords_Global[(kk-1)*3+3-1]=t_celldata.Zmin+(k-1)*2*dz;
            }
            // Then for middle type layer
            for(j=1;j<=_Ny+1;++j){
                for(i=1;i<=_Nx+1;++i){
                    kk=(j-1)*(_Nx+1)+i+(k-1)*(nLayer1Nodes+nLayer2Nodes)+nLayer1Nodes;
                    t_celldata.NodeCoords_Global[(kk-1)*3+1-1]=t_celldata.Xmin+(i-1)*2*dx;
                    t_celldata.NodeCoords_Global[(kk-1)*3+2-1]=t_celldata.Ymin+(j-1)*2*dy;
                    t_celldata.NodeCoords_Global[(kk-1)*3+3-1]=t_celldata.Zmin+(k-1)*2*dz+dz;
                }
            }
        }
        // for the last top layer
        k=_Nz+1;
        for(j=1;j<=_Ny;++j){
            // for bottom line of each element
            for(i=1;i<=2*_Nx+1;i++){
                kk=(j-1)*(2*_Nx+1+_Nx+1)+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                t_celldata.NodeCoords_Global[(kk-1)*3+1-1]=t_celldata.Xmin+(i-1)*dx;
                t_celldata.NodeCoords_Global[(kk-1)*3+2-1]=t_celldata.Ymin+(j-1)*2*dy;
                t_celldata.NodeCoords_Global[(kk-1)*3+3-1]=t_celldata.Zmin+(k-1)*2*dz;
            }
            // for middle line of each element
            for(i=1;i<=_Nx+1;++i){
                kk=(j-1)*(2*_Nx+1+_Nx+1)+2*_Nx+1+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                t_celldata.NodeCoords_Global[(kk-1)*3+1-1]=t_celldata.Xmin+(i-1)*2*dx;
                t_celldata.NodeCoords_Global[(kk-1)*3+2-1]=t_celldata.Ymin+(j-1)*2*dy+dy;
                t_celldata.NodeCoords_Global[(kk-1)*3+3-1]=t_celldata.Zmin+(k-1)*2*dz;
            }
        }
        // for top line
        j=_Ny+1;
        for(i=1;i<=2*_Nx+1;i++){
            kk=(j-1)*(2*_Nx+1+_Nx+1)+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
            t_celldata.NodeCoords_Global[(kk-1)*3+1-1]=t_celldata.Xmin+(i-1)*dx;
            t_celldata.NodeCoords_Global[(kk-1)*3+2-1]=t_celldata.Ymin+(j-1)*2*dy;
            t_celldata.NodeCoords_Global[(kk-1)*3+3-1]=t_celldata.Zmin+(k-1)*2*dz;
        }

        // end-of-node-generation
        t_celldata.MeshCell_Total.resize(t_celldata.BulkElmtsNum);
        leftconn.resize(t_celldata.Ny*t_celldata.Nz);
        rightconn.resize(t_celldata.Ny*t_celldata.Nz);
        //
        bottomconn.resize(t_celldata.Nx*t_celldata.Nz);
        topconn.resize(t_celldata.Nx*t_celldata.Nz);
        //
        backconn.resize(t_celldata.Nx*t_celldata.Ny);
        frontconn.resize(t_celldata.Nx*t_celldata.Ny);

        t_celldata.PhyID2BulkFECellIDMap_Global.clear();
        t_celldata.PhyName2BulkFECellIDMap_Global.clear();
            
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

                    t_celldata.PhyName2BulkFECellIDMap_Global["alldomain"].push_back(e);
                    t_celldata.PhyID2BulkFECellIDMap_Global[0].push_back(e);

                    t_celldata.MeshCell_Total[e-1].Dim=3;
                    t_celldata.MeshCell_Total[e-1].NodesNumPerElmt=20;
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords.resize(20);
                    t_celldata.MeshCell_Total[e-1].VTKCellType=t_celldata.BulkElmtVTKCellType;

                    t_celldata.MeshCell_Total[e-1].ElmtConn.clear();
                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i1);
                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i2);
                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i3);
                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i4);
                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i5);
                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i6);
                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i7);
                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i8);

                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i9);
                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i10);
                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i11);
                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i12);
                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i13);
                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i14);
                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i15);
                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i16);

                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i17);
                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i18);
                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i19);
                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i20);

                    // for element's nodecoords
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(1,1)=t_celldata.NodeCoords_Global[(i1-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(1,2)=t_celldata.NodeCoords_Global[(i1-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(1,3)=t_celldata.NodeCoords_Global[(i1-1)*3+3-1];
                    //
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(2,1)=t_celldata.NodeCoords_Global[(i2-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(2,2)=t_celldata.NodeCoords_Global[(i2-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(2,3)=t_celldata.NodeCoords_Global[(i2-1)*3+3-1];
                    //
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(3,1)=t_celldata.NodeCoords_Global[(i3-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(3,2)=t_celldata.NodeCoords_Global[(i3-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(3,3)=t_celldata.NodeCoords_Global[(i3-1)*3+3-1];
                    //
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(4,1)=t_celldata.NodeCoords_Global[(i4-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(4,2)=t_celldata.NodeCoords_Global[(i4-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(4,3)=t_celldata.NodeCoords_Global[(i4-1)*3+3-1];
                    //
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(5,1)=t_celldata.NodeCoords_Global[(i5-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(5,2)=t_celldata.NodeCoords_Global[(i5-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(5,3)=t_celldata.NodeCoords_Global[(i5-1)*3+3-1];
                    //
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(6,1)=t_celldata.NodeCoords_Global[(i6-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(6,2)=t_celldata.NodeCoords_Global[(i6-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(6,3)=t_celldata.NodeCoords_Global[(i6-1)*3+3-1];
                    //
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(7,1)=t_celldata.NodeCoords_Global[(i7-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(7,2)=t_celldata.NodeCoords_Global[(i7-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(7,3)=t_celldata.NodeCoords_Global[(i7-1)*3+3-1];
                    //
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(8,1)=t_celldata.NodeCoords_Global[(i8-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(8,2)=t_celldata.NodeCoords_Global[(i8-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(8,3)=t_celldata.NodeCoords_Global[(i8-1)*3+3-1];
                    //
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(9,1)=t_celldata.NodeCoords_Global[(i9-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(9,2)=t_celldata.NodeCoords_Global[(i9-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(9,3)=t_celldata.NodeCoords_Global[(i9-1)*3+3-1];
                    //
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(10,1)=t_celldata.NodeCoords_Global[(i10-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(10,2)=t_celldata.NodeCoords_Global[(i10-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(10,3)=t_celldata.NodeCoords_Global[(i10-1)*3+3-1];
                    //
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(11,1)=t_celldata.NodeCoords_Global[(i11-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(11,2)=t_celldata.NodeCoords_Global[(i11-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(11,3)=t_celldata.NodeCoords_Global[(i11-1)*3+3-1];
                    //
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(12,1)=t_celldata.NodeCoords_Global[(i12-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(12,2)=t_celldata.NodeCoords_Global[(i12-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(12,3)=t_celldata.NodeCoords_Global[(i12-1)*3+3-1];
                    //
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(13,1)=t_celldata.NodeCoords_Global[(i13-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(13,2)=t_celldata.NodeCoords_Global[(i13-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(13,3)=t_celldata.NodeCoords_Global[(i13-1)*3+3-1];
                    //
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(14,1)=t_celldata.NodeCoords_Global[(i14-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(14,2)=t_celldata.NodeCoords_Global[(i14-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(14,3)=t_celldata.NodeCoords_Global[(i14-1)*3+3-1];
                    //
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(15,1)=t_celldata.NodeCoords_Global[(i15-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(15,2)=t_celldata.NodeCoords_Global[(i15-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(15,3)=t_celldata.NodeCoords_Global[(i15-1)*3+3-1];
                    //
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(16,1)=t_celldata.NodeCoords_Global[(i16-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(16,2)=t_celldata.NodeCoords_Global[(i16-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(16,3)=t_celldata.NodeCoords_Global[(i16-1)*3+3-1];
                    //
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(17,1)=t_celldata.NodeCoords_Global[(i17-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(17,2)=t_celldata.NodeCoords_Global[(i17-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(17,3)=t_celldata.NodeCoords_Global[(i17-1)*3+3-1];
                    //
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(18,1)=t_celldata.NodeCoords_Global[(i18-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(18,2)=t_celldata.NodeCoords_Global[(i18-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(18,3)=t_celldata.NodeCoords_Global[(i18-1)*3+3-1];
                    //
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(19,1)=t_celldata.NodeCoords_Global[(i19-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(19,2)=t_celldata.NodeCoords_Global[(i19-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(19,3)=t_celldata.NodeCoords_Global[(i19-1)*3+3-1];
                    //
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(20,1)=t_celldata.NodeCoords_Global[(i20-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(20,2)=t_celldata.NodeCoords_Global[(i20-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(20,3)=t_celldata.NodeCoords_Global[(i20-1)*3+3-1];

                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords0=t_celldata.MeshCell_Total[e-1].ElmtNodeCoords;
                        
                    if(i==1){
                        // for left bc elements
                        leftconn[(k-1)*t_celldata.Ny+j-1].Dim=2;
                        leftconn[(k-1)*t_celldata.Ny+j-1].NodesNumPerElmt=8;
                        leftconn[(k-1)*t_celldata.Ny+j-1].VTKCellType=23;
                            
                        leftconn[(k-1)*t_celldata.Ny+j-1].PhysicalGroupNums=1;
                        leftconn[(k-1)*t_celldata.Ny+j-1].PhysicalIDList.clear();
                        leftconn[(k-1)*t_celldata.Ny+j-1].PhysicalIDList.push_back(1);
                        leftconn[(k-1)*t_celldata.Ny+j-1].PhysicalNameList.clear();
                        leftconn[(k-1)*t_celldata.Ny+j-1].PhysicalNameList.push_back("left");

                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtConn.clear();
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtConn.push_back(i1);
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtConn.push_back(i5);
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtConn.push_back(i8);
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtConn.push_back(i4);
                        //
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtConn.push_back(i17);
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtConn.push_back(i16);
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtConn.push_back(i20);
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtConn.push_back(i12);

                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords.resize(8);
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(1,1)=t_celldata.NodeCoords_Global[(i1-1)*3+1-1];
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(1,2)=t_celldata.NodeCoords_Global[(i1-1)*3+2-1];
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(1,3)=t_celldata.NodeCoords_Global[(i1-1)*3+3-1];
                        //
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(2,1)=t_celldata.NodeCoords_Global[(i5-1)*3+1-1];
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(2,2)=t_celldata.NodeCoords_Global[(i5-1)*3+2-1];
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(2,3)=t_celldata.NodeCoords_Global[(i5-1)*3+3-1];
                        //
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(3,1)=t_celldata.NodeCoords_Global[(i8-1)*3+1-1];
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(3,2)=t_celldata.NodeCoords_Global[(i8-1)*3+2-1];
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(3,3)=t_celldata.NodeCoords_Global[(i8-1)*3+3-1];
                        //
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(4,1)=t_celldata.NodeCoords_Global[(i4-1)*3+1-1];
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(4,2)=t_celldata.NodeCoords_Global[(i4-1)*3+2-1];
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(4,3)=t_celldata.NodeCoords_Global[(i4-1)*3+3-1];
                        ////
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(5,1)=t_celldata.NodeCoords_Global[(i17-1)*3+1-1];
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(5,2)=t_celldata.NodeCoords_Global[(i17-1)*3+2-1];
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(5,3)=t_celldata.NodeCoords_Global[(i17-1)*3+3-1];
                        //
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(6,1)=t_celldata.NodeCoords_Global[(i16-1)*3+1-1];
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(6,2)=t_celldata.NodeCoords_Global[(i16-1)*3+2-1];
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(6,3)=t_celldata.NodeCoords_Global[(i16-1)*3+3-1];
                        //
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(7,1)=t_celldata.NodeCoords_Global[(i20-1)*3+1-1];
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(7,2)=t_celldata.NodeCoords_Global[(i20-1)*3+2-1];
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(7,3)=t_celldata.NodeCoords_Global[(i20-1)*3+3-1];
                        //
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(8,1)=t_celldata.NodeCoords_Global[(i12-1)*3+1-1];
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(8,2)=t_celldata.NodeCoords_Global[(i12-1)*3+2-1];
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(8,3)=t_celldata.NodeCoords_Global[(i12-1)*3+3-1];
                        //
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords0=leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords;

                        leftnodes.push_back(i1);
                        leftnodes.push_back(i5);
                        leftnodes.push_back(i8);
                        leftnodes.push_back(i4);
                        leftnodes.push_back(i17);
                        leftnodes.push_back(i16);
                        leftnodes.push_back(i20);
                        leftnodes.push_back(i12);
                    }
                    if(i==t_celldata.Nx){
                        // for right bc elements
                        rightconn[(k-1)*t_celldata.Ny+j-1].Dim=2;
                        rightconn[(k-1)*t_celldata.Ny+j-1].NodesNumPerElmt=8;
                        rightconn[(k-1)*t_celldata.Ny+j-1].VTKCellType=23;
                            
                        rightconn[(k-1)*t_celldata.Ny+j-1].PhysicalGroupNums=1;
                        rightconn[(k-1)*t_celldata.Ny+j-1].PhysicalIDList.clear();
                        rightconn[(k-1)*t_celldata.Ny+j-1].PhysicalIDList.push_back(2);
                        rightconn[(k-1)*t_celldata.Ny+j-1].PhysicalNameList.clear();
                        rightconn[(k-1)*t_celldata.Ny+j-1].PhysicalNameList.push_back("right");

                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtConn.clear();
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtConn.push_back(i2);
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtConn.push_back(i3);
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtConn.push_back(i7);
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtConn.push_back(i6);
                        //
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtConn.push_back(i10);
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtConn.push_back(i19);
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtConn.push_back(i14);
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtConn.push_back(i18);

                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords.resize(8);
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(1,1)=t_celldata.NodeCoords_Global[(i2-1)*3+1-1];
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(1,2)=t_celldata.NodeCoords_Global[(i2-1)*3+2-1];
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(1,3)=t_celldata.NodeCoords_Global[(i2-1)*3+3-1];
                        //
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(2,1)=t_celldata.NodeCoords_Global[(i3-1)*3+1-1];
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(2,2)=t_celldata.NodeCoords_Global[(i3-1)*3+2-1];
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(2,3)=t_celldata.NodeCoords_Global[(i3-1)*3+3-1];
                        //
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(3,1)=t_celldata.NodeCoords_Global[(i7-1)*3+1-1];
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(3,2)=t_celldata.NodeCoords_Global[(i7-1)*3+2-1];
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(3,3)=t_celldata.NodeCoords_Global[(i7-1)*3+3-1];
                        //
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(4,1)=t_celldata.NodeCoords_Global[(i6-1)*3+1-1];
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(4,2)=t_celldata.NodeCoords_Global[(i6-1)*3+2-1];
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(4,3)=t_celldata.NodeCoords_Global[(i6-1)*3+3-1];
                        ////
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(5,1)=t_celldata.NodeCoords_Global[(i10-1)*3+1-1];
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(5,2)=t_celldata.NodeCoords_Global[(i10-1)*3+2-1];
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(5,3)=t_celldata.NodeCoords_Global[(i10-1)*3+3-1];
                        //
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(6,1)=t_celldata.NodeCoords_Global[(i19-1)*3+1-1];
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(6,2)=t_celldata.NodeCoords_Global[(i19-1)*3+2-1];
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(6,3)=t_celldata.NodeCoords_Global[(i19-1)*3+3-1];
                        //
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(7,1)=t_celldata.NodeCoords_Global[(i14-1)*3+1-1];
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(7,2)=t_celldata.NodeCoords_Global[(i14-1)*3+2-1];
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(7,3)=t_celldata.NodeCoords_Global[(i14-1)*3+3-1];
                        //
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(8,1)=t_celldata.NodeCoords_Global[(i18-1)*3+1-1];
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(8,2)=t_celldata.NodeCoords_Global[(i18-1)*3+2-1];
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(8,3)=t_celldata.NodeCoords_Global[(i18-1)*3+3-1];
                        //
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords0=rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords;

                        rightnodes.push_back(i2);
                        rightnodes.push_back(i3);
                        rightnodes.push_back(i7);
                        rightnodes.push_back(i6);
                        rightnodes.push_back(i10);
                        rightnodes.push_back(i19);
                        rightnodes.push_back(i14);
                        rightnodes.push_back(i18);
                    }
                    if(j==1){
                        // for bottom bc elements
                        bottomconn[(k-1)*t_celldata.Nx+i-1].Dim=2;
                        bottomconn[(k-1)*t_celldata.Nx+i-1].NodesNumPerElmt=8;
                        bottomconn[(k-1)*t_celldata.Nx+i-1].VTKCellType=23;
                            
                        bottomconn[(k-1)*t_celldata.Nx+i-1].PhysicalGroupNums=1;
                        bottomconn[(k-1)*t_celldata.Nx+i-1].PhysicalIDList.clear();
                        bottomconn[(k-1)*t_celldata.Nx+i-1].PhysicalIDList.push_back(3);
                        bottomconn[(k-1)*t_celldata.Nx+i-1].PhysicalNameList.clear();
                        bottomconn[(k-1)*t_celldata.Nx+i-1].PhysicalNameList.push_back("bottom");

                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtConn.clear();
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i1);
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i2);
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i6);
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i5);
                        //
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i9);
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i18);
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i13);
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i17);

                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords.resize(8);
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(1,1)=t_celldata.NodeCoords_Global[(i1-1)*3+1-1];
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(1,2)=t_celldata.NodeCoords_Global[(i1-1)*3+2-1];
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(1,3)=t_celldata.NodeCoords_Global[(i1-1)*3+3-1];
                        //
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(2,1)=t_celldata.NodeCoords_Global[(i2-1)*3+1-1];
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(2,2)=t_celldata.NodeCoords_Global[(i2-1)*3+2-1];
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(2,3)=t_celldata.NodeCoords_Global[(i2-1)*3+3-1];
                        //
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(3,1)=t_celldata.NodeCoords_Global[(i6-1)*3+1-1];
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(3,2)=t_celldata.NodeCoords_Global[(i6-1)*3+2-1];
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(3,3)=t_celldata.NodeCoords_Global[(i6-1)*3+3-1];
                        //
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(4,1)=t_celldata.NodeCoords_Global[(i5-1)*3+1-1];
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(4,2)=t_celldata.NodeCoords_Global[(i5-1)*3+2-1];
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(4,3)=t_celldata.NodeCoords_Global[(i5-1)*3+3-1];
                        ////
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(5,1)=t_celldata.NodeCoords_Global[(i9-1)*3+1-1];
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(5,2)=t_celldata.NodeCoords_Global[(i9-1)*3+2-1];
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(5,3)=t_celldata.NodeCoords_Global[(i9-1)*3+3-1];
                        //
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(6,1)=t_celldata.NodeCoords_Global[(i18-1)*3+1-1];
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(6,2)=t_celldata.NodeCoords_Global[(i18-1)*3+2-1];
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(6,3)=t_celldata.NodeCoords_Global[(i18-1)*3+3-1];
                        //
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(7,1)=t_celldata.NodeCoords_Global[(i13-1)*3+1-1];
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(7,2)=t_celldata.NodeCoords_Global[(i13-1)*3+2-1];
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(7,3)=t_celldata.NodeCoords_Global[(i13-1)*3+3-1];
                        //
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(8,1)=t_celldata.NodeCoords_Global[(i17-1)*3+1-1];
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(8,2)=t_celldata.NodeCoords_Global[(i17-1)*3+2-1];
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(8,3)=t_celldata.NodeCoords_Global[(i17-1)*3+3-1];
                        //
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords0=bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords;

                        bottomnodes.push_back(i1);
                        bottomnodes.push_back(i2);
                        bottomnodes.push_back(i6);
                        bottomnodes.push_back(i5);
                        bottomnodes.push_back(i9 );
                        bottomnodes.push_back(i18);
                        bottomnodes.push_back(i13);
                        bottomnodes.push_back(i17);
                    }
                    if(j==t_celldata.Ny){
                        // for top bc elements
                        topconn[(k-1)*t_celldata.Nx+i-1].Dim=2;
                        topconn[(k-1)*t_celldata.Nx+i-1].NodesNumPerElmt=8;
                        topconn[(k-1)*t_celldata.Nx+i-1].VTKCellType=23;
                            
                        topconn[(k-1)*t_celldata.Nx+i-1].PhysicalGroupNums=1;
                        topconn[(k-1)*t_celldata.Nx+i-1].PhysicalIDList.clear();
                        topconn[(k-1)*t_celldata.Nx+i-1].PhysicalIDList.push_back(4);
                        topconn[(k-1)*t_celldata.Nx+i-1].PhysicalNameList.clear();
                        topconn[(k-1)*t_celldata.Nx+i-1].PhysicalNameList.push_back("top");

                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtConn.clear();
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i4);
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i8);
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i7);
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i3);
                        //
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i20);
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i15);
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i19);
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i11);

                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords.resize(8);
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(1,1)=t_celldata.NodeCoords_Global[(i4-1)*3+1-1];
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(1,2)=t_celldata.NodeCoords_Global[(i4-1)*3+2-1];
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(1,3)=t_celldata.NodeCoords_Global[(i4-1)*3+3-1];
                        //
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(2,1)=t_celldata.NodeCoords_Global[(i8-1)*3+1-1];
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(2,2)=t_celldata.NodeCoords_Global[(i8-1)*3+2-1];
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(2,3)=t_celldata.NodeCoords_Global[(i8-1)*3+3-1];
                        //
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(3,1)=t_celldata.NodeCoords_Global[(i7-1)*3+1-1];
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(3,2)=t_celldata.NodeCoords_Global[(i7-1)*3+2-1];
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(3,3)=t_celldata.NodeCoords_Global[(i7-1)*3+3-1];
                        //
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(4,1)=t_celldata.NodeCoords_Global[(i3-1)*3+1-1];
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(4,2)=t_celldata.NodeCoords_Global[(i3-1)*3+2-1];
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(4,3)=t_celldata.NodeCoords_Global[(i3-1)*3+3-1];
                        ////
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(5,1)=t_celldata.NodeCoords_Global[(i20-1)*3+1-1];
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(5,2)=t_celldata.NodeCoords_Global[(i20-1)*3+2-1];
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(5,3)=t_celldata.NodeCoords_Global[(i20-1)*3+3-1];
                        //
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(6,1)=t_celldata.NodeCoords_Global[(i15-1)*3+1-1];
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(6,2)=t_celldata.NodeCoords_Global[(i15-1)*3+2-1];
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(6,3)=t_celldata.NodeCoords_Global[(i15-1)*3+3-1];
                        //
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(7,1)=t_celldata.NodeCoords_Global[(i19-1)*3+1-1];
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(7,2)=t_celldata.NodeCoords_Global[(i19-1)*3+2-1];
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(7,3)=t_celldata.NodeCoords_Global[(i19-1)*3+3-1];
                        //
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(8,1)=t_celldata.NodeCoords_Global[(i11-1)*3+1-1];
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(8,2)=t_celldata.NodeCoords_Global[(i11-1)*3+2-1];
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(8,3)=t_celldata.NodeCoords_Global[(i11-1)*3+3-1];
                        //
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords0=topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords;

                        topnodes.push_back(i4);
                        topnodes.push_back(i8);
                        topnodes.push_back(i7);
                        topnodes.push_back(i3);
                        topnodes.push_back(i20);
                        topnodes.push_back(i15);
                        topnodes.push_back(i19);
                        topnodes.push_back(i11);
                    }
                    if(k==1){
                        // for back bc elements
                        backconn[(j-1)*t_celldata.Nx+i-1].Dim=2;
                        backconn[(j-1)*t_celldata.Nx+i-1].NodesNumPerElmt=8;
                        backconn[(j-1)*t_celldata.Nx+i-1].VTKCellType=23;
                            
                        backconn[(j-1)*t_celldata.Nx+i-1].PhysicalGroupNums=1;
                        backconn[(j-1)*t_celldata.Nx+i-1].PhysicalIDList.clear();
                        backconn[(j-1)*t_celldata.Nx+i-1].PhysicalIDList.push_back(5);
                        backconn[(j-1)*t_celldata.Nx+i-1].PhysicalNameList.clear();
                        backconn[(j-1)*t_celldata.Nx+i-1].PhysicalNameList.push_back("back");

                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtConn.clear();
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i1);
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i4);
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i3);
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i2);
                        //
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i12);
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i11);
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i10);
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i9);

                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords.resize(8);
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(1,1)=t_celldata.NodeCoords_Global[(i1-1)*3+1-1];
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(1,2)=t_celldata.NodeCoords_Global[(i1-1)*3+2-1];
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(1,3)=t_celldata.NodeCoords_Global[(i1-1)*3+3-1];
                        //
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(2,1)=t_celldata.NodeCoords_Global[(i4-1)*3+1-1];
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(2,2)=t_celldata.NodeCoords_Global[(i4-1)*3+2-1];
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(2,3)=t_celldata.NodeCoords_Global[(i4-1)*3+3-1];
                        //
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(3,1)=t_celldata.NodeCoords_Global[(i3-1)*3+1-1];
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(3,2)=t_celldata.NodeCoords_Global[(i3-1)*3+2-1];
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(3,3)=t_celldata.NodeCoords_Global[(i3-1)*3+3-1];
                        //
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(4,1)=t_celldata.NodeCoords_Global[(i2-1)*3+1-1];
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(4,2)=t_celldata.NodeCoords_Global[(i2-1)*3+2-1];
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(4,3)=t_celldata.NodeCoords_Global[(i2-1)*3+3-1];
                        ////
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(5,1)=t_celldata.NodeCoords_Global[(i12-1)*3+1-1];
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(5,2)=t_celldata.NodeCoords_Global[(i12-1)*3+2-1];
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(5,3)=t_celldata.NodeCoords_Global[(i12-1)*3+3-1];
                        //
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(6,1)=t_celldata.NodeCoords_Global[(i11-1)*3+1-1];
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(6,2)=t_celldata.NodeCoords_Global[(i11-1)*3+2-1];
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(6,3)=t_celldata.NodeCoords_Global[(i11-1)*3+3-1];
                        //
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(7,1)=t_celldata.NodeCoords_Global[(i10-1)*3+1-1];
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(7,2)=t_celldata.NodeCoords_Global[(i10-1)*3+2-1];
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(7,3)=t_celldata.NodeCoords_Global[(i10-1)*3+3-1];
                        //
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(8,1)=t_celldata.NodeCoords_Global[(i9-1)*3+1-1];
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(8,2)=t_celldata.NodeCoords_Global[(i9-1)*3+2-1];
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(8,3)=t_celldata.NodeCoords_Global[(i9-1)*3+3-1];
                        //
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords0=backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords;

                        backnodes.push_back(i1);
                        backnodes.push_back(i4);
                        backnodes.push_back(i3);
                        backnodes.push_back(i2);
                        backnodes.push_back(i12);
                        backnodes.push_back(i11);
                        backnodes.push_back(i10);
                        backnodes.push_back(i9 );
                    }
                    if(k==t_celldata.Nz){
                        // for front bc elements
                        frontconn[(j-1)*t_celldata.Nx+i-1].Dim=2;
                        frontconn[(j-1)*t_celldata.Nx+i-1].NodesNumPerElmt=8;
                        frontconn[(j-1)*t_celldata.Nx+i-1].VTKCellType=23;
                            
                        frontconn[(j-1)*t_celldata.Nx+i-1].PhysicalGroupNums=1;
                        frontconn[(j-1)*t_celldata.Nx+i-1].PhysicalIDList.clear();
                        frontconn[(j-1)*t_celldata.Nx+i-1].PhysicalIDList.push_back(6);
                        frontconn[(j-1)*t_celldata.Nx+i-1].PhysicalNameList.clear();
                        frontconn[(j-1)*t_celldata.Nx+i-1].PhysicalNameList.push_back("front");

                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtConn.clear();
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i5);
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i6);
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i7);
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i8);
                        //
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i13);
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i14);
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i15);
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtConn.push_back(i16);

                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords.resize(8);
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(1,1)=t_celldata.NodeCoords_Global[(i5-1)*3+1-1];
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(1,2)=t_celldata.NodeCoords_Global[(i5-1)*3+2-1];
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(1,3)=t_celldata.NodeCoords_Global[(i5-1)*3+3-1];
                        //
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(2,1)=t_celldata.NodeCoords_Global[(i6-1)*3+1-1];
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(2,2)=t_celldata.NodeCoords_Global[(i6-1)*3+2-1];
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(2,3)=t_celldata.NodeCoords_Global[(i6-1)*3+3-1];
                        //
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(3,1)=t_celldata.NodeCoords_Global[(i7-1)*3+1-1];
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(3,2)=t_celldata.NodeCoords_Global[(i7-1)*3+2-1];
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(3,3)=t_celldata.NodeCoords_Global[(i7-1)*3+3-1];
                        //
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(4,1)=t_celldata.NodeCoords_Global[(i8-1)*3+1-1];
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(4,2)=t_celldata.NodeCoords_Global[(i8-1)*3+2-1];
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(4,3)=t_celldata.NodeCoords_Global[(i8-1)*3+3-1];
                        ////
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(5,1)=t_celldata.NodeCoords_Global[(i13-1)*3+1-1];
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(5,2)=t_celldata.NodeCoords_Global[(i13-1)*3+2-1];
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(5,3)=t_celldata.NodeCoords_Global[(i13-1)*3+3-1];
                        //
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(6,1)=t_celldata.NodeCoords_Global[(i14-1)*3+1-1];
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(6,2)=t_celldata.NodeCoords_Global[(i14-1)*3+2-1];
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(6,3)=t_celldata.NodeCoords_Global[(i14-1)*3+3-1];
                        //
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(7,1)=t_celldata.NodeCoords_Global[(i15-1)*3+1-1];
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(7,2)=t_celldata.NodeCoords_Global[(i15-1)*3+2-1];
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(7,3)=t_celldata.NodeCoords_Global[(i15-1)*3+3-1];
                        //
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(8,1)=t_celldata.NodeCoords_Global[(i16-1)*3+1-1];
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(8,2)=t_celldata.NodeCoords_Global[(i16-1)*3+2-1];
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(8,3)=t_celldata.NodeCoords_Global[(i16-1)*3+3-1];
                        //
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords0=frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords;

                        frontnodes.push_back(i5);
                        frontnodes.push_back(i6);
                        frontnodes.push_back(i7);
                        frontnodes.push_back(i8);
                        frontnodes.push_back(i13);
                        frontnodes.push_back(i14);
                        frontnodes.push_back(i15);
                        frontnodes.push_back(i16);
                    }
                }
            }
        } // end-of-element-generation-loop
    
        // remove the duplicate node id, for nodal type set, you don't need these duplicate node ids!
        sort(leftnodes.begin(),leftnodes.end());
        leftnodes.erase(unique(leftnodes.begin(),leftnodes.end()),leftnodes.end());
        // for rightnodes
        sort(rightnodes.begin(),rightnodes.end());
        rightnodes.erase(unique(rightnodes.begin(),rightnodes.end()),rightnodes.end());
        // for bottomnodes
        sort(bottomnodes.begin(),bottomnodes.end());
        bottomnodes.erase(unique(bottomnodes.begin(),bottomnodes.end()),bottomnodes.end());
        // for topnodes
        sort(topnodes.begin(),topnodes.end());
        topnodes.erase(unique(topnodes.begin(),topnodes.end()),topnodes.end());
        // for backnodes
        sort(backnodes.begin(),backnodes.end());
        backnodes.erase(unique(backnodes.begin(),backnodes.end()),backnodes.end());
        // for frontnodes
        sort(frontnodes.begin(),frontnodes.end());
        frontnodes.erase(unique(frontnodes.begin(),frontnodes.end()),frontnodes.end());

        // setup the physical group information
        t_celldata.PhyGroupNum_Global=1+6;
        t_celldata.PhyDimVector_Global.resize(1+6,0);
        t_celldata.PhyIDVector_Global.resize(1+6,0);
        t_celldata.PhyNameVector_Global.resize(1+6);
        t_celldata.PhyGroupElmtsNumVector_Global.resize(1+6,0);

        // for physical dim vector
        t_celldata.PhyDimVector_Global[0]=3;
        t_celldata.PhyDimVector_Global[1]=2;
        t_celldata.PhyDimVector_Global[2]=2;
        t_celldata.PhyDimVector_Global[3]=2;
        t_celldata.PhyDimVector_Global[4]=2;
        t_celldata.PhyDimVector_Global[5]=2;
        t_celldata.PhyDimVector_Global[6]=2;

        // for physical id vector
        t_celldata.PhyIDVector_Global[0]=0;
        t_celldata.PhyIDVector_Global[1]=1;
        t_celldata.PhyIDVector_Global[2]=2;
        t_celldata.PhyIDVector_Global[3]=3;
        t_celldata.PhyIDVector_Global[4]=4;
        t_celldata.PhyIDVector_Global[5]=5;
        t_celldata.PhyIDVector_Global[6]=6;

        // for physical name vector
        t_celldata.PhyNameVector_Global[0]="alldomain";
        t_celldata.PhyNameVector_Global[1]="left";
        t_celldata.PhyNameVector_Global[2]="right";
        t_celldata.PhyNameVector_Global[3]="bottom";
        t_celldata.PhyNameVector_Global[4]="top";
        t_celldata.PhyNameVector_Global[5]="back";
        t_celldata.PhyNameVector_Global[6]="front";

        // for phy group elmts num vector
        t_celldata.PhyGroupElmtsNumVector_Global[0]=static_cast<int>(t_celldata.BulkElmtsNum);
        t_celldata.PhyGroupElmtsNumVector_Global[1]=static_cast<int>(leftconn.size());
        t_celldata.PhyGroupElmtsNumVector_Global[2]=static_cast<int>(rightconn.size());
        t_celldata.PhyGroupElmtsNumVector_Global[3]=static_cast<int>(bottomconn.size());
        t_celldata.PhyGroupElmtsNumVector_Global[4]=static_cast<int>(topconn.size());
        t_celldata.PhyGroupElmtsNumVector_Global[5]=static_cast<int>(backconn.size());
        t_celldata.PhyGroupElmtsNumVector_Global[6]=static_cast<int>(frontconn.size());

        /**
         * setup id<---->name map
        */
        // id--->name map
        t_celldata.PhyID2NameMap_Global[0]="alldomain";
        t_celldata.PhyID2NameMap_Global[1]="left";
        t_celldata.PhyID2NameMap_Global[2]="right";
        t_celldata.PhyID2NameMap_Global[3]="bottom";
        t_celldata.PhyID2NameMap_Global[4]="top";
        t_celldata.PhyID2NameMap_Global[5]="back";
        t_celldata.PhyID2NameMap_Global[6]="front";
        // name--->id map
        t_celldata.PhyName2IDMap_Global["alldomain"]=0;
        t_celldata.PhyName2IDMap_Global["left"]=1;
        t_celldata.PhyName2IDMap_Global["right"]=2;
        t_celldata.PhyName2IDMap_Global["bottom"]=3;
        t_celldata.PhyName2IDMap_Global["top"]=4;
        t_celldata.PhyName2IDMap_Global["back"]=5;
        t_celldata.PhyName2IDMap_Global["front"]=6;

        /**
         * Setup the global mapping, this should only be nonzero on master rank !!!
        */
        t_celldata.PhyName2MeshCellVectorMap_Global["alldomain"]=t_celldata.MeshCell_Total;
        t_celldata.PhyName2MeshCellVectorMap_Global["left"]=leftconn;
        t_celldata.PhyName2MeshCellVectorMap_Global["right"]=rightconn;
        t_celldata.PhyName2MeshCellVectorMap_Global["bottom"]=bottomconn;
        t_celldata.PhyName2MeshCellVectorMap_Global["top"]=topconn;
        t_celldata.PhyName2MeshCellVectorMap_Global["back"]=backconn;
        t_celldata.PhyName2MeshCellVectorMap_Global["front"]=frontconn;

        t_celldata.PhyID2MeshCellVectorMap_Global[0]=t_celldata.MeshCell_Total;
        t_celldata.PhyID2MeshCellVectorMap_Global[1]=leftconn;
        t_celldata.PhyID2MeshCellVectorMap_Global[2]=rightconn;
        t_celldata.PhyID2MeshCellVectorMap_Global[3]=bottomconn;
        t_celldata.PhyID2MeshCellVectorMap_Global[4]=topconn;
        t_celldata.PhyID2MeshCellVectorMap_Global[5]=backconn;
        t_celldata.PhyID2MeshCellVectorMap_Global[6]=frontconn;

        /**
         * Setup the nodal physical info group
        */
        t_celldata.NodalPhyGroupNum_Global=6;
        t_celldata.NodalPhyIDVector_Global.resize(6);
        t_celldata.NodalPhyNameVector_Global.resize(6);
        t_celldata.NodalPhyGroupNodesNumVector_Global.resize(6);

        t_celldata.NodalPhyIDVector_Global[0]=10001;
        t_celldata.NodalPhyIDVector_Global[1]=10002;
        t_celldata.NodalPhyIDVector_Global[2]=10003;
        t_celldata.NodalPhyIDVector_Global[3]=10004;
        t_celldata.NodalPhyIDVector_Global[4]=10005;
        t_celldata.NodalPhyIDVector_Global[5]=10006;

        t_celldata.NodalPhyNameVector_Global[0]="leftnodes";
        t_celldata.NodalPhyNameVector_Global[1]="rightnodes";
        t_celldata.NodalPhyNameVector_Global[2]="bottomnodes";
        t_celldata.NodalPhyNameVector_Global[3]="topnodes";
        t_celldata.NodalPhyNameVector_Global[4]="backnodes";
        t_celldata.NodalPhyNameVector_Global[5]="frontnodes";

        t_celldata.NodalPhyGroupNodesNumVector_Global[0]=static_cast<int>(leftnodes.size());
        t_celldata.NodalPhyGroupNodesNumVector_Global[1]=static_cast<int>(rightnodes.size());
        t_celldata.NodalPhyGroupNodesNumVector_Global[2]=static_cast<int>(bottomnodes.size());
        t_celldata.NodalPhyGroupNodesNumVector_Global[3]=static_cast<int>(topnodes.size());
        t_celldata.NodalPhyGroupNodesNumVector_Global[4]=static_cast<int>(backnodes.size());
        t_celldata.NodalPhyGroupNodesNumVector_Global[5]=static_cast<int>(frontnodes.size());
        
        t_celldata.NodalPhyID2NameMap_Global[10001]="leftnodes";
        t_celldata.NodalPhyID2NameMap_Global[10002]="rightnodes";
        t_celldata.NodalPhyID2NameMap_Global[10003]="bottomnodes";
        t_celldata.NodalPhyID2NameMap_Global[10004]="topnodes";
        t_celldata.NodalPhyID2NameMap_Global[10005]="backnodes";
        t_celldata.NodalPhyID2NameMap_Global[10006]="frontnodes";
        
        t_celldata.NodalPhyName2IDMap_Global["leftnodes"]=10001;
        t_celldata.NodalPhyName2IDMap_Global["rightnodes"]=10002;
        t_celldata.NodalPhyName2IDMap_Global["bottomnodes"]=10003;
        t_celldata.NodalPhyName2IDMap_Global["topnodes"]=10004;
        t_celldata.NodalPhyName2IDMap_Global["backnodes"]=10005;
        t_celldata.NodalPhyName2IDMap_Global["frontnodes"]=10006;

        t_celldata.NodalPhyName2NodeIDVecMap_Global["leftnodes"]=leftnodes;
        t_celldata.NodalPhyName2NodeIDVecMap_Global["rightnodes"]=rightnodes;
        t_celldata.NodalPhyName2NodeIDVecMap_Global["bottomnodes"]=bottomnodes;
        t_celldata.NodalPhyName2NodeIDVecMap_Global["topnodes"]=topnodes;
        t_celldata.NodalPhyName2NodeIDVecMap_Global["backnodes"]=backnodes;
        t_celldata.NodalPhyName2NodeIDVecMap_Global["frontnodes"]=frontnodes;

        // for mesh cell partition info
        t_celldata.BulkCellPartionInfo_Global.resize(t_celldata.BulkElmtsNum,0);

    }// end-of-if(rank==0)
    else{
        t_celldata.MeshOrder=2;
        t_celldata.BulkMeshTypeName="hex20";

        t_celldata.BulkElmtsNum=t_celldata.Nx*t_celldata.Ny*t_celldata.Nz;
        t_celldata.LineElmtsNum=0;
        t_celldata.SurfElmtsNum=2*(t_celldata.Nx*t_celldata.Ny
                                  +t_celldata.Nx*t_celldata.Nz
                                  +t_celldata.Ny*t_celldata.Nz);
        t_celldata.ElmtsNum=t_celldata.BulkElmtsNum
                           +t_celldata.SurfElmtsNum
                           +t_celldata.LineElmtsNum;

        int nLayer1Nodes=(2*t_celldata.Nx+1)*(2*t_celldata.Ny+1)-t_celldata.Nx*t_celldata.Ny;// for norm layer
        int nLayer2Nodes=(t_celldata.Nx+1)*(t_celldata.Ny+1);          // for middle layer

        t_celldata.NodesNum=nLayer1Nodes*(t_celldata.Nz+1)+nLayer2Nodes*t_celldata.Nz;
        t_celldata.NodesNumPerBulkElmt=20;
        t_celldata.NodesNumPerSurfElmt=8;
        t_celldata.NodesNumPerLineElmt=3;

        t_celldata.BulkElmtVTKCellType=25;

        t_celldata.LineElmtVTKCellType=4;
        t_celldata.LineElmtMeshType=MeshType::EDGE3;

        t_celldata.SurfElmtVTKCellType=23;
        t_celldata.SurfElmtMeshType=MeshType::QUAD8;
        // here one must setup the basic phy group info, which is required by the 'readElmtBlock' function in readInputFile class !!!
        // you don't need to partition the mesh, instead, the basic phy group info should be setup
        // setup the physical group information
        t_celldata.PhyGroupNum_Global=1+6;
        t_celldata.PhyDimVector_Global.resize(1+6,0);
        t_celldata.PhyIDVector_Global.resize(1+6,0);
        t_celldata.PhyNameVector_Global.resize(1+6);
        t_celldata.PhyGroupElmtsNumVector_Global.resize(1+6,0);

        // for physical dim vector
        t_celldata.PhyDimVector_Global[0]=3;
        t_celldata.PhyDimVector_Global[1]=2;
        t_celldata.PhyDimVector_Global[2]=2;
        t_celldata.PhyDimVector_Global[3]=2;
        t_celldata.PhyDimVector_Global[4]=2;
        t_celldata.PhyDimVector_Global[5]=2;
        t_celldata.PhyDimVector_Global[6]=2;

        // for physical id vector
        t_celldata.PhyIDVector_Global[0]=0;
        t_celldata.PhyIDVector_Global[1]=1;
        t_celldata.PhyIDVector_Global[2]=2;
        t_celldata.PhyIDVector_Global[3]=3;
        t_celldata.PhyIDVector_Global[4]=4;
        t_celldata.PhyIDVector_Global[5]=5;
        t_celldata.PhyIDVector_Global[6]=6;

        // for physical name vector
        t_celldata.PhyNameVector_Global[0]="alldomain";
        t_celldata.PhyNameVector_Global[1]="left";
        t_celldata.PhyNameVector_Global[2]="right";
        t_celldata.PhyNameVector_Global[3]="bottom";
        t_celldata.PhyNameVector_Global[4]="top";
        t_celldata.PhyNameVector_Global[5]="back";
        t_celldata.PhyNameVector_Global[6]="front";
    }

    return true;
}