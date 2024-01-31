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
//+++ Purpose: the 3D hex8 lagrange FE cell generator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FECell/Lagrange3DHex8MeshCellGenerator.h"
#include "FECell/SingleMeshCell.h"
#include "MPIUtils/MPIDataBus.h"

Lagrange3DHex8MeshCellGenerator::Lagrange3DHex8MeshCellGenerator(){
    m_mesh_generated=false;
}
Lagrange3DHex8MeshCellGenerator::~Lagrange3DHex8MeshCellGenerator(){
    m_mesh_generated=false;
}
//*********************************************************
bool Lagrange3DHex8MeshCellGenerator::generateFECell(FECellData &t_celldata){
    int rank,size;
    t_celldata.MaxDim=3;
    t_celldata.MinDim=2;

    t_celldata.ActiveDofsNum=0;
    t_celldata.TotalDofsNum=0;
    t_celldata.MaxDofsPerNode=0;

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
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
        int i1,i2,i3,i4,i5,i6,i7,i8;
        
        // generate the mesh cell for hex8 mesh
        t_celldata.MeshOrder=1;
        t_celldata.BulkMeshTypeName="hex8";
        dx=(t_celldata.Xmax-t_celldata.Xmin)/t_celldata.Nx;
        dy=(t_celldata.Ymax-t_celldata.Ymin)/t_celldata.Ny;
        dz=(t_celldata.Zmax-t_celldata.Zmin)/t_celldata.Nz;
            
        t_celldata.BulkElmtsNum=t_celldata.Nx*t_celldata.Ny*t_celldata.Nz;
        t_celldata.LineElmtsNum=0;
        t_celldata.SurfElmtsNum=2*(t_celldata.Nx*t_celldata.Ny
                                  +t_celldata.Nx*t_celldata.Nz
                                  +t_celldata.Ny*t_celldata.Nz);
        t_celldata.ElmtsNum=t_celldata.BulkElmtsNum
                           +t_celldata.SurfElmtsNum
                           +t_celldata.LineElmtsNum;
                               
        t_celldata.NodesNum=(t_celldata.Nx+1)*(t_celldata.Ny+1)*(t_celldata.Nz+1);
        t_celldata.NodesNumPerBulkElmt=8;
        t_celldata.NodesNumPerSurfElmt=4;
        t_celldata.NodesNumPerLineElmt=2;
            
        t_celldata.BulkElmtVTKCellType=12;
        t_celldata.BulkElmtMeshType=MeshType::HEX8;

        t_celldata.LineElmtVTKCellType=3;
        t_celldata.LineElmtMeshType=MeshType::EDGE2;
            
        t_celldata.SurfElmtVTKCellType=9;
        t_celldata.SurfElmtMeshType=MeshType::QUAD4;

        vector<double> nodecoords;
        nodecoords.resize(t_celldata.NodesNum*3,0.0);
        leftnodes.clear();rightnodes.clear();
        bottomnodes.clear();rightnodes.clear();
        backnodes.clear();frontnodes.clear();
        for(k=1;k<=t_celldata.Nz+1;k++){
            for(j=1;j<=t_celldata.Ny+1;j++){
                for(i=1;i<=t_celldata.Nx+1;i++){
                    kk=(j-1)*(t_celldata.Nx+1)+i+(k-1)*(t_celldata.Nx+1)*(t_celldata.Ny+1);
                        
                    nodecoords[(kk-1)*3+1-1]=t_celldata.Xmin+(i-1)*dx;
                    nodecoords[(kk-1)*3+2-1]=t_celldata.Ymin+(j-1)*dy;
                    nodecoords[(kk-1)*3+3-1]=t_celldata.Zmin+(k-1)*dz;
                        
                    if(i==1){
                        // for left side nodes
                        leftnodes.push_back(kk);// global id, start from 1
                    }
                    if(i==t_celldata.Nx+1){
                        // for right side nodes
                        rightnodes.push_back(kk);
                    }
                    if(j==1){
                        // for bottom side nodes
                        bottomnodes.push_back(kk);
                    }
                    if(j==t_celldata.Ny+1){
                        // for top side nodes
                        topnodes.push_back(kk);
                    }
                    if(k==1){
                        // for back side nodes
                        backnodes.push_back(kk);
                    }
                    if(k==t_celldata.Nz+1){
                        // for front side nodes
                        frontnodes.push_back(kk);
                    }
                }
            }
        }// end-of-node-generation

        t_celldata.NodeCoords_Global.clear();
        for(const auto &it:nodecoords) t_celldata.NodeCoords_Global.push_back(it);

        // for the connectivity information of bulk elements
        t_celldata.MeshCell_Total.resize(t_celldata.BulkElmtsNum);
        leftconn.resize(t_celldata.Ny*t_celldata.Nz);
        rightconn.resize(t_celldata.Ny*t_celldata.Nz);
        //
        bottomconn.resize(t_celldata.Nx*t_celldata.Nz);
        topconn.resize(t_celldata.Nx*t_celldata.Nz);
        //
        backconn.resize(t_celldata.Nx*t_celldata.Ny);
        frontconn.resize(t_celldata.Nx*t_celldata.Ny);
            
        leftnodes.clear();rightnodes.clear();
        bottomnodes.clear();rightnodes.clear();
        backnodes.clear();frontnodes.clear();
            
        for(k=1;k<=t_celldata.Nz;k++){
            for(j=1;j<=t_celldata.Ny;j++){
                for(i=1;i<=t_celldata.Nx;i++){
                    e=(j-1)*t_celldata.Nx+i+(k-1)*t_celldata.Nx*t_celldata.Ny;
                    i1=(j-1)*(t_celldata.Nx+1)+i+(k-1)*(t_celldata.Nx+1)*(t_celldata.Ny+1);
                    i2=i1+1;
                    i3=i2+t_celldata.Nx+1;
                    i4=i3-1;
                    i5=i1+(t_celldata.Nx+1)*(t_celldata.Ny+1);
                    i6=i2+(t_celldata.Nx+1)*(t_celldata.Ny+1);
                    i7=i3+(t_celldata.Nx+1)*(t_celldata.Ny+1);
                    i8=i4+(t_celldata.Nx+1)*(t_celldata.Ny+1);

                    t_celldata.MeshCell_Total[e-1].Dim=3;
                    t_celldata.MeshCell_Total[e-1].NodesNumPerElmt=8;
                    t_celldata.MeshCell_Total[e-1].VTKCellType=t_celldata.BulkElmtVTKCellType;
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords.resize(8);
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords0.resize(8);

                    t_celldata.MeshCell_Total[e-1].ElmtConn.clear();
                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i1);
                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i2);
                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i3);
                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i4);
                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i5);
                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i6);
                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i7);
                    t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i8);
                    // assign node coords
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(1,1)=nodecoords[(i1-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(1,2)=nodecoords[(i1-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(1,3)=nodecoords[(i1-1)*3+3-1];
                    //
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(2,1)=nodecoords[(i2-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(2,2)=nodecoords[(i2-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(2,3)=nodecoords[(i2-1)*3+3-1];
                    //
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(3,1)=nodecoords[(i3-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(3,2)=nodecoords[(i3-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(3,3)=nodecoords[(i3-1)*3+3-1];
                    //
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(4,1)=nodecoords[(i4-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(4,2)=nodecoords[(i4-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(4,3)=nodecoords[(i4-1)*3+3-1];
                    //
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(5,1)=nodecoords[(i5-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(5,2)=nodecoords[(i5-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(5,3)=nodecoords[(i5-1)*3+3-1];
                    //
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(6,1)=nodecoords[(i6-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(6,2)=nodecoords[(i6-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(6,3)=nodecoords[(i6-1)*3+3-1];
                    //
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(7,1)=nodecoords[(i7-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(7,2)=nodecoords[(i7-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(7,3)=nodecoords[(i7-1)*3+3-1];
                    //
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(8,1)=nodecoords[(i8-1)*3+1-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(8,2)=nodecoords[(i8-1)*3+2-1];
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords(8,3)=nodecoords[(i8-1)*3+3-1];
                        
                    t_celldata.MeshCell_Total[e-1].ElmtNodeCoords0=t_celldata.MeshCell_Total[e-1].ElmtNodeCoords;
                        
                    if(i==1){
                        // for left bc elements
                        leftconn[(k-1)*t_celldata.Ny+j-1].Dim=2;
                        leftconn[(k-1)*t_celldata.Ny+j-1].NodesNumPerElmt=4;
                        leftconn[(k-1)*t_celldata.Ny+j-1].VTKCellType=9;

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

                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords.resize(4);
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(1,1)=nodecoords[(i1-1)*3+1-1];
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(1,2)=nodecoords[(i1-1)*3+2-1];
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(1,3)=nodecoords[(i1-1)*3+3-1];
                        //
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(2,1)=nodecoords[(i5-1)*3+1-1];
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(2,2)=nodecoords[(i5-1)*3+2-1];
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(2,3)=nodecoords[(i5-1)*3+3-1];
                        //
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(3,1)=nodecoords[(i8-1)*3+1-1];
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(3,2)=nodecoords[(i8-1)*3+2-1];
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(3,3)=nodecoords[(i8-1)*3+3-1];
                        //
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(4,1)=nodecoords[(i4-1)*3+1-1];
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(4,2)=nodecoords[(i4-1)*3+2-1];
                        leftconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(4,3)=nodecoords[(i4-1)*3+3-1];

                        leftnodes.push_back(i1);
                        leftnodes.push_back(i5);
                        leftnodes.push_back(i8);
                        leftnodes.push_back(i4);
                    }
                    if(i==t_celldata.Nx){
                        // for right bc elements
                        rightconn[(k-1)*t_celldata.Ny+j-1].Dim=2;
                        rightconn[(k-1)*t_celldata.Ny+j-1].NodesNumPerElmt=4;
                        rightconn[(k-1)*t_celldata.Ny+j-1].VTKCellType=9;
                            
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

                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords.resize(4);
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(1,1)=nodecoords[(i2-1)*3+1-1];
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(1,2)=nodecoords[(i2-1)*3+2-1];
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(1,3)=nodecoords[(i2-1)*3+3-1];
                        //
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(2,1)=nodecoords[(i3-1)*3+1-1];
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(2,2)=nodecoords[(i3-1)*3+2-1];
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(2,3)=nodecoords[(i3-1)*3+3-1];
                        //
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(3,1)=nodecoords[(i7-1)*3+1-1];
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(3,2)=nodecoords[(i7-1)*3+2-1];
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(3,3)=nodecoords[(i7-1)*3+3-1];
                        //
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(4,1)=nodecoords[(i6-1)*3+1-1];
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(4,2)=nodecoords[(i6-1)*3+2-1];
                        rightconn[(k-1)*t_celldata.Ny+j-1].ElmtNodeCoords(4,3)=nodecoords[(i6-1)*3+3-1];

                        rightnodes.push_back(i2);
                        rightnodes.push_back(i3);
                        rightnodes.push_back(i7);
                        rightnodes.push_back(i6);
                    }
                    if(j==1){
                        // for bottom bc elements
                        bottomconn[(k-1)*t_celldata.Nx+i-1].Dim=2;
                        bottomconn[(k-1)*t_celldata.Nx+i-1].NodesNumPerElmt=4;
                        bottomconn[(k-1)*t_celldata.Nx+i-1].VTKCellType=9;
                            
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

                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords.resize(4);
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(1,1)=nodecoords[(i1-1)*3+1-1];
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(1,2)=nodecoords[(i1-1)*3+2-1];
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(1,3)=nodecoords[(i1-1)*3+3-1];
                        //
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(2,1)=nodecoords[(i2-1)*3+1-1];
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(2,2)=nodecoords[(i2-1)*3+2-1];
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(2,3)=nodecoords[(i2-1)*3+3-1];
                        //
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(3,1)=nodecoords[(i6-1)*3+1-1];
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(3,2)=nodecoords[(i6-1)*3+2-1];
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(3,3)=nodecoords[(i6-1)*3+3-1];
                        //
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(4,1)=nodecoords[(i5-1)*3+1-1];
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(4,2)=nodecoords[(i5-1)*3+2-1];
                        bottomconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(4,3)=nodecoords[(i5-1)*3+3-1];

                        bottomnodes.push_back(i1);
                        bottomnodes.push_back(i2);
                        bottomnodes.push_back(i6);
                        bottomnodes.push_back(i5);
                    }
                    if(j==t_celldata.Ny){
                        // for top bc elements
                        topconn[(k-1)*t_celldata.Nx+i-1].Dim=2;
                        topconn[(k-1)*t_celldata.Nx+i-1].NodesNumPerElmt=4;
                        topconn[(k-1)*t_celldata.Nx+i-1].VTKCellType=9;
                            
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

                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords.resize(4);
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(1,1)=nodecoords[(i4-1)*3+1-1];
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(1,2)=nodecoords[(i4-1)*3+2-1];
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(1,3)=nodecoords[(i4-1)*3+3-1];
                        //
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(2,1)=nodecoords[(i8-1)*3+1-1];
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(2,2)=nodecoords[(i8-1)*3+2-1];
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(2,3)=nodecoords[(i8-1)*3+3-1];
                        //
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(3,1)=nodecoords[(i7-1)*3+1-1];
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(3,2)=nodecoords[(i7-1)*3+2-1];
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(3,3)=nodecoords[(i7-1)*3+3-1];
                        //
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(4,1)=nodecoords[(i3-1)*3+1-1];
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(4,2)=nodecoords[(i3-1)*3+2-1];
                        topconn[(k-1)*t_celldata.Nx+i-1].ElmtNodeCoords(4,3)=nodecoords[(i3-1)*3+3-1];

                        topnodes.push_back(i4);
                        topnodes.push_back(i8);
                        topnodes.push_back(i7);
                        topnodes.push_back(i3);
                    }
                    if(k==1){
                        // for back bc elements
                        backconn[(j-1)*t_celldata.Nx+i-1].Dim=2;
                        backconn[(j-1)*t_celldata.Nx+i-1].NodesNumPerElmt=4;
                        backconn[(j-1)*t_celldata.Nx+i-1].VTKCellType=9;
                            
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

                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords.resize(4);
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(1,1)=nodecoords[(i1-1)*3+1-1];
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(1,2)=nodecoords[(i1-1)*3+2-1];
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(1,3)=nodecoords[(i1-1)*3+3-1];
                        //
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(2,1)=nodecoords[(i4-1)*3+1-1];
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(2,2)=nodecoords[(i4-1)*3+2-1];
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(2,3)=nodecoords[(i4-1)*3+3-1];
                        //
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(3,1)=nodecoords[(i3-1)*3+1-1];
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(3,2)=nodecoords[(i3-1)*3+2-1];
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(3,3)=nodecoords[(i3-1)*3+3-1];
                        //
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(4,1)=nodecoords[(i2-1)*3+1-1];
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(4,2)=nodecoords[(i2-1)*3+2-1];
                        backconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(4,3)=nodecoords[(i2-1)*3+3-1];

                        backnodes.push_back(i1);
                        backnodes.push_back(i4);
                        backnodes.push_back(i3);
                        backnodes.push_back(i2);
                    }
                    if(k==t_celldata.Nz){
                        // for front bc elements
                        frontconn[(j-1)*t_celldata.Nx+i-1].Dim=2;
                        frontconn[(j-1)*t_celldata.Nx+i-1].NodesNumPerElmt=4;
                        frontconn[(j-1)*t_celldata.Nx+i-1].VTKCellType=9;
                            
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

                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords.resize(4);
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(1,1)=nodecoords[(i5-1)*3+1-1];
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(1,2)=nodecoords[(i5-1)*3+2-1];
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(1,3)=nodecoords[(i5-1)*3+3-1];
                        //
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(2,1)=nodecoords[(i6-1)*3+1-1];
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(2,2)=nodecoords[(i6-1)*3+2-1];
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(2,3)=nodecoords[(i6-1)*3+3-1];
                        //
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(3,1)=nodecoords[(i7-1)*3+1-1];
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(3,2)=nodecoords[(i7-1)*3+2-1];
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(3,3)=nodecoords[(i7-1)*3+3-1];
                        //
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(4,1)=nodecoords[(i8-1)*3+1-1];
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(4,2)=nodecoords[(i8-1)*3+2-1];
                        frontconn[(j-1)*t_celldata.Nx+i-1].ElmtNodeCoords(4,3)=nodecoords[(i8-1)*3+3-1];

                        frontnodes.push_back(i5);
                        frontnodes.push_back(i6);
                        frontnodes.push_back(i7);
                        frontnodes.push_back(i8);
                    }
                }
            }
        }// end-of-element-generation-loop
    
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

        /**
         * Now we start to distribute the gloabl mesh into different ranks
        */
        // send out the physical group info
        int cpuid;

        // send out the total mesh cell
        int iStart,iEnd,ranksize;
        vector<SingleMeshCell> LocalCellVector;
        vector<int> nodeids;

        t_celldata.PhyID2MeshCellVectorMap_Local.clear();
        t_celldata.PhyName2MeshCellVectorMap_Local.clear();

        for(cpuid=0;cpuid<size;cpuid++){
            ranksize=t_celldata.BulkElmtsNum/size;
            iStart=cpuid*ranksize;
            iEnd=(cpuid+1)*ranksize;
            if(cpuid==size-1) iEnd=t_celldata.BulkElmtsNum;

            LocalCellVector.clear();
            for(int e=iStart;e<iEnd;e++){
                LocalCellVector.push_back(t_celldata.MeshCell_Total[e]);
            }
            if(cpuid==0){
                t_celldata.MeshCell_Local=LocalCellVector;
                t_celldata.PhyID2MeshCellVectorMap_Local[0]=LocalCellVector;
                t_celldata.PhyName2MeshCellVectorMap_Local["alldomain"]=t_celldata.MeshCell_Local;// each local rank share the same phy name as the master rank!!!
            }
            else{
                MPIDataBus::sendMeshCellToOthers(LocalCellVector,1000*cpuid,cpuid);
                MPIDataBus::sendPhyID2MeshCellMapToOthers(0,LocalCellVector,1000*cpuid+20,cpuid);
                MPIDataBus::sendPhyName2MeshCellMapToOthers("alldomain",LocalCellVector,1000*cpuid+40,cpuid);
            }

            // for left conn
            ranksize=static_cast<int>(leftconn.size())/size;
            iStart=cpuid*ranksize;
            iEnd=(cpuid+1)*ranksize;
            if(cpuid==size-1) iEnd=static_cast<int>(leftconn.size());
            LocalCellVector.clear();
            for(int e=iStart;e<iEnd;e++){
                LocalCellVector.push_back(leftconn[e]);
            }
            if(cpuid==0){
                t_celldata.PhyID2MeshCellVectorMap_Local[1]=LocalCellVector;
                t_celldata.PhyName2MeshCellVectorMap_Local["left"]=LocalCellVector;
            }
            else{
                MPIDataBus::sendPhyID2MeshCellMapToOthers(1,LocalCellVector,1000*cpuid+60,cpuid);
                MPIDataBus::sendPhyName2MeshCellMapToOthers("left",LocalCellVector,1000*cpuid+80,cpuid);
            }
            // for right conn
            ranksize=static_cast<int>(rightconn.size())/size;
            iStart=cpuid*ranksize;
            iEnd=(cpuid+1)*ranksize;
            if(cpuid==size-1) iEnd=static_cast<int>(rightconn.size());
            LocalCellVector.clear();
            for(int e=iStart;e<iEnd;e++){
                LocalCellVector.push_back(rightconn[e]);
            }
            if(cpuid==0){
                t_celldata.PhyID2MeshCellVectorMap_Local[2]=LocalCellVector;
                t_celldata.PhyName2MeshCellVectorMap_Local["right"]=LocalCellVector;
            }
            else{
                MPIDataBus::sendPhyID2MeshCellMapToOthers(2,LocalCellVector,1000*cpuid+100,cpuid);
                MPIDataBus::sendPhyName2MeshCellMapToOthers("right",LocalCellVector,1000*cpuid+120,cpuid);
            }

            // for bottom conn
            ranksize=static_cast<int>(bottomconn.size())/size;
            iStart=cpuid*ranksize;
            iEnd=(cpuid+1)*ranksize;
            if(cpuid==size-1) iEnd=static_cast<int>(bottomconn.size());
            LocalCellVector.clear();
            for(int e=iStart;e<iEnd;e++){
                LocalCellVector.push_back(bottomconn[e]);
            }
            if(cpuid==0){
                t_celldata.PhyID2MeshCellVectorMap_Local[3]=LocalCellVector;
                t_celldata.PhyName2MeshCellVectorMap_Local["bottom"]=LocalCellVector;
            }
            else{
                MPIDataBus::sendPhyID2MeshCellMapToOthers(3,LocalCellVector,1000*cpuid+140,cpuid);
                MPIDataBus::sendPhyName2MeshCellMapToOthers("bottom",LocalCellVector,1000*cpuid+160,cpuid);
            }
            // for top conn
            ranksize=static_cast<int>(topconn.size())/size;
            iStart=cpuid*ranksize;
            iEnd=(cpuid+1)*ranksize;
            if(cpuid==size-1) iEnd=static_cast<int>(topconn.size());
            LocalCellVector.clear();
            for(int e=iStart;e<iEnd;e++){
                LocalCellVector.push_back(topconn[e]);
            }
            if(cpuid==0){
                t_celldata.PhyID2MeshCellVectorMap_Local[4]=LocalCellVector;
                t_celldata.PhyName2MeshCellVectorMap_Local["top"]=LocalCellVector;
            }
            else{
                MPIDataBus::sendPhyID2MeshCellMapToOthers(4,LocalCellVector,1000*cpuid+180,cpuid);
                MPIDataBus::sendPhyName2MeshCellMapToOthers("top",LocalCellVector,1000*cpuid+200,cpuid);
            }

            // for back conn
            ranksize=static_cast<int>(backconn.size())/size;
            iStart=cpuid*ranksize;
            iEnd=(cpuid+1)*ranksize;
            if(cpuid==size-1) iEnd=static_cast<int>(backconn.size());
            LocalCellVector.clear();
            for(int e=iStart;e<iEnd;e++){
                LocalCellVector.push_back(backconn[e]);
            }
            if(cpuid==0){
                t_celldata.PhyID2MeshCellVectorMap_Local[5]=LocalCellVector;
                t_celldata.PhyName2MeshCellVectorMap_Local["back"]=LocalCellVector;
            }
            else{
                MPIDataBus::sendPhyID2MeshCellMapToOthers(5,LocalCellVector,1000*cpuid+220,cpuid);
                MPIDataBus::sendPhyName2MeshCellMapToOthers("back",LocalCellVector,1000*cpuid+240,cpuid);
            }
            // for front conn
            ranksize=static_cast<int>(frontconn.size())/size;
            iStart=cpuid*ranksize;
            iEnd=(cpuid+1)*ranksize;
            if(cpuid==size-1) iEnd=static_cast<int>(frontconn.size());
            LocalCellVector.clear();
            for(int e=iStart;e<iEnd;e++){
                LocalCellVector.push_back(frontconn[e]);
            }
            if(cpuid==0){
                t_celldata.PhyID2MeshCellVectorMap_Local[6]=LocalCellVector;
                t_celldata.PhyName2MeshCellVectorMap_Local["front"]=LocalCellVector;
            }
            else{
                MPIDataBus::sendPhyID2MeshCellMapToOthers(6,LocalCellVector,1000*cpuid+260,cpuid);
                MPIDataBus::sendPhyName2MeshCellMapToOthers("front",LocalCellVector,1000*cpuid+280,cpuid);
            }

            //***************************************************
            // for nodal physical groups
            //***************************************************
            // for leftnodes
            ranksize=static_cast<int>(leftnodes.size())/size;
            iStart=cpuid*ranksize;
            iEnd=(cpuid+1)*ranksize;
            if(cpuid==size-1) iEnd=static_cast<int>(leftnodes.size());
            nodeids.clear();
            for(int e=iStart;e<iEnd;e++){
                nodeids.push_back(leftnodes[e]);
            }
            if(cpuid==0){
                t_celldata.NodalPhyName2NodeIDVecMap_Local["leftnodes"]=nodeids;
            }
            else{
                MPIDataBus::sendPhyName2NodeIDVecMapToOthers("leftnodes",nodeids,1000*cpuid+300,cpuid);
            }
            // for rightnodes
            ranksize=static_cast<int>(rightnodes.size())/size;
            iStart=cpuid*ranksize;
            iEnd=(cpuid+1)*ranksize;
            if(cpuid==size-1) iEnd=static_cast<int>(rightnodes.size());
            nodeids.clear();
            for(int e=iStart;e<iEnd;e++){
                nodeids.push_back(rightnodes[e]);
            }
            if(cpuid==0){
                t_celldata.NodalPhyName2NodeIDVecMap_Local["rightnodes"]=nodeids;
            }
            else{
                MPIDataBus::sendPhyName2NodeIDVecMapToOthers("rightnodes",nodeids,1000*cpuid+320,cpuid);
            }
            // for bottomnodes
            ranksize=static_cast<int>(bottomnodes.size())/size;
            iStart=cpuid*ranksize;
            iEnd=(cpuid+1)*ranksize;
            if(cpuid==size-1) iEnd=static_cast<int>(bottomnodes.size());
            nodeids.clear();
            for(int e=iStart;e<iEnd;e++){
                nodeids.push_back(bottomnodes[e]);
            }
            if(cpuid==0){
                t_celldata.NodalPhyName2NodeIDVecMap_Local["bottomnodes"]=nodeids;
            }
            else{
                MPIDataBus::sendPhyName2NodeIDVecMapToOthers("bottomnodes",nodeids,1000*cpuid+340,cpuid);
            }
            // for topnodes
            ranksize=static_cast<int>(topnodes.size())/size;
            iStart=cpuid*ranksize;
            iEnd=(cpuid+1)*ranksize;
            if(cpuid==size-1) iEnd=static_cast<int>(topnodes.size());
            nodeids.clear();
            for(int e=iStart;e<iEnd;e++){
                nodeids.push_back(topnodes[e]);
            }
            if(cpuid==0){
                t_celldata.NodalPhyName2NodeIDVecMap_Local["topnodes"]=nodeids;
            }
            else{
                MPIDataBus::sendPhyName2NodeIDVecMapToOthers("topnodes",nodeids,1000*cpuid+360,cpuid);
            }
            // for backnodes
            ranksize=static_cast<int>(backnodes.size())/size;
            iStart=cpuid*ranksize;
            iEnd=(cpuid+1)*ranksize;
            if(cpuid==size-1) iEnd=static_cast<int>(backnodes.size());
            nodeids.clear();
            for(int e=iStart;e<iEnd;e++){
                nodeids.push_back(backnodes[e]);
            }
            if(cpuid==0){
                t_celldata.NodalPhyName2NodeIDVecMap_Local["backnodes"]=nodeids;
            }
            else{
                MPIDataBus::sendPhyName2NodeIDVecMapToOthers("backnodes",nodeids,1000*cpuid+380,cpuid);
            }
            // for frontnodes
            ranksize=static_cast<int>(frontnodes.size())/size;
            iStart=cpuid*ranksize;
            iEnd=(cpuid+1)*ranksize;
            if(cpuid==size-1) iEnd=static_cast<int>(frontnodes.size());
            nodeids.clear();
            for(int e=iStart;e<iEnd;e++){
                nodeids.push_back(frontnodes[e]);
            }
            if(cpuid==0){
                t_celldata.NodalPhyName2NodeIDVecMap_Local["frontnodes"]=nodeids;
            }
            else{
                MPIDataBus::sendPhyName2NodeIDVecMapToOthers("frontnodes",nodeids,1000*cpuid+400,cpuid);
            }
        }// end-of-cpuid-loop

    }// end-of-if(rank==0)
    else{
        // now we distribute the global mesh into different ranks
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
        
        /**
         * Receive message from master rank
        */
        MPIDataBus::receiveMeshCellFromMaster(t_celldata.MeshCell_Local,1000*rank);
        MPIDataBus::receivePhyID2MeshCellMapFromMaster(t_celldata.PhyID2MeshCellVectorMap_Local,1000*rank+20);
        MPIDataBus::receivePhyName2MeshCellMapFromMaster(t_celldata.PhyName2MeshCellVectorMap_Local,1000*rank+40);

        // for leftconn
        MPIDataBus::receivePhyID2MeshCellMapFromMaster(t_celldata.PhyID2MeshCellVectorMap_Local,1000*rank+60);
        MPIDataBus::receivePhyName2MeshCellMapFromMaster(t_celldata.PhyName2MeshCellVectorMap_Local,1000*rank+80);
        // for rightconn
        MPIDataBus::receivePhyID2MeshCellMapFromMaster(t_celldata.PhyID2MeshCellVectorMap_Local,1000*rank+100);
        MPIDataBus::receivePhyName2MeshCellMapFromMaster(t_celldata.PhyName2MeshCellVectorMap_Local,1000*rank+120);
        // for bottomconn
        MPIDataBus::receivePhyID2MeshCellMapFromMaster(t_celldata.PhyID2MeshCellVectorMap_Local,1000*rank+140);
        MPIDataBus::receivePhyName2MeshCellMapFromMaster(t_celldata.PhyName2MeshCellVectorMap_Local,1000*rank+160);
        // for topconn
        MPIDataBus::receivePhyID2MeshCellMapFromMaster(t_celldata.PhyID2MeshCellVectorMap_Local,1000*rank+180);
        MPIDataBus::receivePhyName2MeshCellMapFromMaster(t_celldata.PhyName2MeshCellVectorMap_Local,1000*rank+200);
        // for backconn
        MPIDataBus::receivePhyID2MeshCellMapFromMaster(t_celldata.PhyID2MeshCellVectorMap_Local,1000*rank+220);
        MPIDataBus::receivePhyName2MeshCellMapFromMaster(t_celldata.PhyName2MeshCellVectorMap_Local,1000*rank+240);
        // for frontconn
        MPIDataBus::receivePhyID2MeshCellMapFromMaster(t_celldata.PhyID2MeshCellVectorMap_Local,1000*rank+260);
        MPIDataBus::receivePhyName2MeshCellMapFromMaster(t_celldata.PhyName2MeshCellVectorMap_Local,1000*rank+280);

        //*** for nodal physical group
        // for leftnodes
        MPIDataBus::receivePhyName2NodeIDVecMapFromMaster(t_celldata.NodalPhyName2NodeIDVecMap_Local,1000*rank+300);
        // for rightnodes
        MPIDataBus::receivePhyName2NodeIDVecMapFromMaster(t_celldata.NodalPhyName2NodeIDVecMap_Local,1000*rank+320);
        // for bottomnodes
        MPIDataBus::receivePhyName2NodeIDVecMapFromMaster(t_celldata.NodalPhyName2NodeIDVecMap_Local,1000*rank+340);
        // for topnodes
        MPIDataBus::receivePhyName2NodeIDVecMapFromMaster(t_celldata.NodalPhyName2NodeIDVecMap_Local,1000*rank+360);
        // for backnodes
        MPIDataBus::receivePhyName2NodeIDVecMapFromMaster(t_celldata.NodalPhyName2NodeIDVecMap_Local,1000*rank+380);
        // for frontnodes
        MPIDataBus::receivePhyName2NodeIDVecMapFromMaster(t_celldata.NodalPhyName2NodeIDVecMap_Local,1000*rank+400);
    }

    return true;
}