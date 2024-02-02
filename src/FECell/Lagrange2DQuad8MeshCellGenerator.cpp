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
//+++ Date   : 2024.02.02
//+++ Purpose: the 2D quad8 lagrange FE cell generator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FECell/Lagrange2DQuad8MeshCellGenerator.h"
#include "FECell/SingleMeshCell.h"
#include "MPIUtils/MPIDataBus.h"

Lagrange2DQuad8MeshCellGenerator::Lagrange2DQuad8MeshCellGenerator(){
    m_mesh_generated=false;
}
Lagrange2DQuad8MeshCellGenerator::~Lagrange2DQuad8MeshCellGenerator(){
    m_mesh_generated=false;
}
//*********************************************************
bool Lagrange2DQuad8MeshCellGenerator::generateFECell(FECellData &t_celldata){
    int rank,size;
    t_celldata.MaxDim=2;
    t_celldata.MinDim=1;

    t_celldata.ActiveDofsNum=0;
    t_celldata.TotalDofsNum=0;
    t_celldata.MaxDofsPerNode=0;

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    if(rank==0){
        vector<SingleMeshCell> leftconn,rightconn;
        vector<SingleMeshCell> bottomconn,topconn;
        vector<int> leftnodes,rightnodes;
        vector<int> bottomnodes,topnodes;

        // only create mesh on the master rank, then distributed them into different ranks !
        double dx,dy;
        int i,j,k,e;
        int i1,i2,i3,i4;
        int i5,i6,i7,i8;
        
        // generate the mesh cell for quad8 mesh
        t_celldata.MeshOrder=2;
        t_celldata.BulkMeshTypeName="quad8";
        dx=(t_celldata.Xmax-t_celldata.Xmin)/(2.0*t_celldata.Nx);
        dy=(t_celldata.Ymax-t_celldata.Ymin)/(2.0*t_celldata.Ny);
            
        t_celldata.BulkElmtsNum=t_celldata.Nx*t_celldata.Ny;
        t_celldata.LineElmtsNum=2*(t_celldata.Nx+t_celldata.Ny);
        t_celldata.SurfElmtsNum=0;
        t_celldata.ElmtsNum=t_celldata.BulkElmtsNum
                           +t_celldata.SurfElmtsNum
                           +t_celldata.LineElmtsNum;
                               
        t_celldata.NodesNum=(2*t_celldata.Nx+1)*(2*t_celldata.Ny+1)-t_celldata.BulkElmtsNum;
        t_celldata.NodesNumPerBulkElmt=8;
        t_celldata.NodesNumPerSurfElmt=0;
        t_celldata.NodesNumPerLineElmt=3;
            
        t_celldata.BulkElmtVTKCellType=23;
        t_celldata.BulkElmtMeshType=MeshType::QUAD8;

        t_celldata.LineElmtVTKCellType=4;
        t_celldata.LineElmtMeshType=MeshType::EDGE3;
            
        t_celldata.SurfElmtVTKCellType=0;
        t_celldata.SurfElmtMeshType=MeshType::NULLTYPE;

        vector<double> nodecoords;
        nodecoords.resize(t_celldata.NodesNum*3,0.0);
        leftnodes.clear();rightnodes.clear();
        bottomnodes.clear();rightnodes.clear();
        for(j=1;j<=t_celldata.Ny;j++){
            // for bottom line of each element
            for(i=1;i<=2*t_celldata.Nx+1;i++){
                k=(j-1)*(2*t_celldata.Nx+1+t_celldata.Nx+1)+i;
                nodecoords[(k-1)*3+1-1]=t_celldata.Xmin+(i-1)*dx;
                nodecoords[(k-1)*3+2-1]=t_celldata.Ymin+(j-1)*2*dy;
                nodecoords[(k-1)*3+3-1]=0.0;
            }
            // for middle line of each element
            for(i=1;i<=t_celldata.Nx+1;i++){
                k=(j-1)*(2*t_celldata.Nx+1+t_celldata.Nx+1)+2*t_celldata.Nx+1+i;
                nodecoords[(k-1)*3+1-1]=t_celldata.Xmin+(i-1)*2*dx;
                nodecoords[(k-1)*3+2-1]=t_celldata.Ymin+(j-1)*2*dy+dy;
                nodecoords[(k-1)*3+3-1]=0.0;
            }
        }
        // for the last top line
        j=t_celldata.Ny+1;
        for(i=1;i<=2*t_celldata.Nx+1;i++){
            k=(j-1)*(2*t_celldata.Nx+1+t_celldata.Nx+1)+i;
            nodecoords[(k-1)*3+1-1]=t_celldata.Xmin+(i-1)*dx;
            nodecoords[(k-1)*3+2-1]=t_celldata.Ymin+(j-1)*2*dy;
            nodecoords[(k-1)*3+3-1]=0.0;
        }// end-of-node-generation

        t_celldata.NodeCoords_Global.clear();
        for(const auto &it:nodecoords) t_celldata.NodeCoords_Global.push_back(it);

        // for the connectivity information of bulk elements
        t_celldata.MeshCell_Total.resize(t_celldata.BulkElmtsNum);
        leftconn.resize(t_celldata.Ny);
        rightconn.resize(t_celldata.Ny);
        //
        bottomconn.resize(t_celldata.Nx);
        topconn.resize(t_celldata.Nx);
            
        leftnodes.clear();rightnodes.clear();
        bottomnodes.clear();rightnodes.clear();

        for(j=1;j<=t_celldata.Ny;j++){
            for(i=1;i<=t_celldata.Nx;i++){
                e=(j-1)*t_celldata.Nx+i;
                i1=(j-1)*(2*t_celldata.Nx+1+t_celldata.Nx+1)+2*i-1;
                i2=i1+2;
                i3=i2+(2*t_celldata.Nx+1+t_celldata.Nx+1);
                i4=i3-2;

                i5=i1+1;
                i6=i2+(2*t_celldata.Nx+1)-i;
                i7=i3-1;
                i8=i1+(2*t_celldata.Nx+1)-(i-1);

                t_celldata.MeshCell_Total[e-1].Dim=2;
                t_celldata.MeshCell_Total[e-1].VTKCellType=23;
                t_celldata.MeshCell_Local[e-1].NodesNumPerElmt=8;
                t_celldata.MeshCell_Total[e-1].ElmtConn.clear();
                t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i1);
                t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i2);
                t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i3);
                t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i4);

                t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i5);
                t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i6);
                t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i7);
                t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i8);

                t_celldata.MeshCell_Total[e-1].ElmtNodeCoords.resize(8);
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

                // for the boundary element and nodes
                // the layout of your quad8 should be:
                // 4--7--3
                // |     |
                // 8     6
                // |     |
                // 1--5--2
                if(j==1){
                    // for bottom bc elements
                    bottomconn[i-1].Dim=1;
                    bottomconn[i-1].NodesNumPerElmt=3;
                    bottomconn[i-1].VTKCellType=4;

                    bottomconn[i-1].ElmtConn.clear();
                    bottomconn[i-1].ElmtConn.push_back(i1);
                    bottomconn[i-1].ElmtConn.push_back(i5);
                    bottomconn[i-1].ElmtConn.push_back(i2);

                    bottomnodes.push_back(i1);
                    bottomnodes.push_back(i5);
                    bottomnodes.push_back(i2);
                }
                if(j==t_celldata.Ny){
                    // for top bc elements
                    topconn[i-1].Dim=1;
                    topconn[i-1].NodesNumPerElmt=3;
                    topconn[i-1].VTKCellType=4;

                    topconn[i-1].ElmtConn.clear();
                    topconn[i-1].ElmtConn.push_back(i3);
                    topconn[i-1].ElmtConn.push_back(i7);
                    topconn[i-1].ElmtConn.push_back(i4);

                    topnodes.push_back(i3);
                    topnodes.push_back(i7);
                    topnodes.push_back(i4);
                }
                if(i==1){
                    // for left bc elements
                    leftconn[j-1].Dim=1;
                    leftconn[j-1].NodesNumPerElmt=3;
                    leftconn[j-1].VTKCellType=4;

                    leftconn[j-1].ElmtConn.clear();
                    leftconn[j-1].ElmtConn.push_back(i4);
                    leftconn[j-1].ElmtConn.push_back(i8);
                    leftconn[j-1].ElmtConn.push_back(i1);

                    leftnodes.push_back(i4);
                    leftnodes.push_back(i8);
                    leftnodes.push_back(i1);
                }
                if(i==t_celldata.Nx){
                    // for right bc elements
                    rightconn[j-1].Dim=1;
                    rightconn[j-1].NodesNumPerElmt=3;
                    rightconn[j-1].VTKCellType=4;

                    rightconn[j-1].ElmtConn.clear();
                    rightconn[j-1].ElmtConn.push_back(i2);
                    rightconn[j-1].ElmtConn.push_back(i6);
                    rightconn[j-1].ElmtConn.push_back(i3);

                    rightnodes.push_back(i2);
                    rightnodes.push_back(i6);
                    rightnodes.push_back(i3);
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

        // setup the physical group information
        t_celldata.PhyGroupNum_Global=1+4;
        t_celldata.PhyDimVector_Global.resize(1+4,0);
        t_celldata.PhyIDVector_Global.resize(1+4,0);
        t_celldata.PhyNameVector_Global.resize(1+4);
        t_celldata.PhyGroupElmtsNumVector_Global.resize(1+4,0);

        // for physical dim vector
        t_celldata.PhyDimVector_Global[0]=2;
        t_celldata.PhyDimVector_Global[1]=1;
        t_celldata.PhyDimVector_Global[2]=1;
        t_celldata.PhyDimVector_Global[3]=1;
        t_celldata.PhyDimVector_Global[4]=1;

        // for physical id vector
        t_celldata.PhyIDVector_Global[0]=0;
        t_celldata.PhyIDVector_Global[1]=1;
        t_celldata.PhyIDVector_Global[2]=2;
        t_celldata.PhyIDVector_Global[3]=3;
        t_celldata.PhyIDVector_Global[4]=4;

        // for physical name vector
        t_celldata.PhyNameVector_Global[0]="alldomain";
        t_celldata.PhyNameVector_Global[1]="left";
        t_celldata.PhyNameVector_Global[2]="right";
        t_celldata.PhyNameVector_Global[3]="bottom";
        t_celldata.PhyNameVector_Global[4]="top";
       
        // for phy group elmts num vector
        t_celldata.PhyGroupElmtsNumVector_Global[0]=static_cast<int>(t_celldata.BulkElmtsNum);
        t_celldata.PhyGroupElmtsNumVector_Global[1]=static_cast<int>(leftconn.size());
        t_celldata.PhyGroupElmtsNumVector_Global[2]=static_cast<int>(rightconn.size());
        t_celldata.PhyGroupElmtsNumVector_Global[3]=static_cast<int>(bottomconn.size());
        t_celldata.PhyGroupElmtsNumVector_Global[4]=static_cast<int>(topconn.size());

        /**
         * setup id<---->name map
        */
        // id--->name map
        t_celldata.PhyID2NameMap_Global[0]="alldomain";
        t_celldata.PhyID2NameMap_Global[1]="left";
        t_celldata.PhyID2NameMap_Global[2]="right";
        t_celldata.PhyID2NameMap_Global[3]="bottom";
        t_celldata.PhyID2NameMap_Global[4]="top";
        // name--->id map
        t_celldata.PhyName2IDMap_Global["alldomain"]=0;
        t_celldata.PhyName2IDMap_Global["left"]=1;
        t_celldata.PhyName2IDMap_Global["right"]=2;
        t_celldata.PhyName2IDMap_Global["bottom"]=3;
        t_celldata.PhyName2IDMap_Global["top"]=4;

        /**
         * Setup the global mapping, this should only be nonzero on master rank !!!
        */
        t_celldata.PhyName2MeshCellVectorMap_Global["alldomain"]=t_celldata.MeshCell_Total;
        t_celldata.PhyName2MeshCellVectorMap_Global["left"]=leftconn;
        t_celldata.PhyName2MeshCellVectorMap_Global["right"]=rightconn;
        t_celldata.PhyName2MeshCellVectorMap_Global["bottom"]=bottomconn;
        t_celldata.PhyName2MeshCellVectorMap_Global["top"]=topconn;

        t_celldata.PhyID2MeshCellVectorMap_Global[0]=t_celldata.MeshCell_Total;
        t_celldata.PhyID2MeshCellVectorMap_Global[1]=leftconn;
        t_celldata.PhyID2MeshCellVectorMap_Global[2]=rightconn;
        t_celldata.PhyID2MeshCellVectorMap_Global[3]=bottomconn;
        t_celldata.PhyID2MeshCellVectorMap_Global[4]=topconn;

        /**
         * Setup the nodal physical info group
        */
        t_celldata.NodalPhyGroupNum_Global=4;
        t_celldata.NodalPhyIDVector_Global.resize(4);
        t_celldata.NodalPhyNameVector_Global.resize(4);
        t_celldata.NodalPhyGroupNodesNumVector_Global.resize(4);

        t_celldata.NodalPhyIDVector_Global[0]=10001;
        t_celldata.NodalPhyIDVector_Global[1]=10002;
        t_celldata.NodalPhyIDVector_Global[2]=10003;
        t_celldata.NodalPhyIDVector_Global[3]=10004;

        t_celldata.NodalPhyNameVector_Global[0]="leftnodes";
        t_celldata.NodalPhyNameVector_Global[1]="rightnodes";
        t_celldata.NodalPhyNameVector_Global[2]="bottomnodes";
        t_celldata.NodalPhyNameVector_Global[3]="topnodes";

        t_celldata.NodalPhyGroupNodesNumVector_Global[0]=static_cast<int>(leftnodes.size());
        t_celldata.NodalPhyGroupNodesNumVector_Global[1]=static_cast<int>(rightnodes.size());
        t_celldata.NodalPhyGroupNodesNumVector_Global[2]=static_cast<int>(bottomnodes.size());
        t_celldata.NodalPhyGroupNodesNumVector_Global[3]=static_cast<int>(topnodes.size());
        
        t_celldata.NodalPhyID2NameMap_Global[10001]="leftnodes";
        t_celldata.NodalPhyID2NameMap_Global[10002]="rightnodes";
        t_celldata.NodalPhyID2NameMap_Global[10003]="bottomnodes";
        t_celldata.NodalPhyID2NameMap_Global[10004]="topnodes";
        
        t_celldata.NodalPhyName2IDMap_Global["leftnodes"]=10001;
        t_celldata.NodalPhyName2IDMap_Global["rightnodes"]=10002;
        t_celldata.NodalPhyName2IDMap_Global["bottomnodes"]=10003;
        t_celldata.NodalPhyName2IDMap_Global["topnodes"]=10004;

        t_celldata.NodalPhyName2NodeIDVecMap_Global["leftnodes"]=leftnodes;
        t_celldata.NodalPhyName2NodeIDVecMap_Global["rightnodes"]=rightnodes;
        t_celldata.NodalPhyName2NodeIDVecMap_Global["bottomnodes"]=bottomnodes;
        t_celldata.NodalPhyName2NodeIDVecMap_Global["topnodes"]=topnodes;

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
        }// end-of-cpuid-loop

    }// end-of-if(rank==0)
    else{
        // now we distribute the global mesh into different ranks
        // setup the physical group information
        t_celldata.PhyGroupNum_Global=1+4;
        t_celldata.PhyDimVector_Global.resize(1+4,0);
        t_celldata.PhyIDVector_Global.resize(1+4,0);
        t_celldata.PhyNameVector_Global.resize(1+4);
        t_celldata.PhyGroupElmtsNumVector_Global.resize(1+4,0);

        // for physical dim vector
        t_celldata.PhyDimVector_Global[0]=2;
        t_celldata.PhyDimVector_Global[1]=1;
        t_celldata.PhyDimVector_Global[2]=1;
        t_celldata.PhyDimVector_Global[3]=1;
        t_celldata.PhyDimVector_Global[4]=1;

        // for physical id vector
        t_celldata.PhyIDVector_Global[0]=0;
        t_celldata.PhyIDVector_Global[1]=1;
        t_celldata.PhyIDVector_Global[2]=2;
        t_celldata.PhyIDVector_Global[3]=3;
        t_celldata.PhyIDVector_Global[4]=4;

        // for physical name vector
        t_celldata.PhyNameVector_Global[0]="alldomain";
        t_celldata.PhyNameVector_Global[1]="left";
        t_celldata.PhyNameVector_Global[2]="right";
        t_celldata.PhyNameVector_Global[3]="bottom";
        t_celldata.PhyNameVector_Global[4]="top";

        /**
         * setup id<---->name map
        */
        // id--->name map
        t_celldata.PhyID2NameMap_Global[0]="alldomain";
        t_celldata.PhyID2NameMap_Global[1]="left";
        t_celldata.PhyID2NameMap_Global[2]="right";
        t_celldata.PhyID2NameMap_Global[3]="bottom";
        t_celldata.PhyID2NameMap_Global[4]="top";
        // name--->id map
        t_celldata.PhyName2IDMap_Global["alldomain"]=0;
        t_celldata.PhyName2IDMap_Global["left"]=1;
        t_celldata.PhyName2IDMap_Global["right"]=2;
        t_celldata.PhyName2IDMap_Global["bottom"]=3;
        t_celldata.PhyName2IDMap_Global["top"]=4;

        /**
         * Setup the nodal physical info group
        */
        t_celldata.NodalPhyGroupNum_Global=4;
        t_celldata.NodalPhyIDVector_Global.resize(4);
        t_celldata.NodalPhyNameVector_Global.resize(4);
        t_celldata.NodalPhyGroupNodesNumVector_Global.resize(4);

        t_celldata.NodalPhyIDVector_Global[0]=10001;
        t_celldata.NodalPhyIDVector_Global[1]=10002;
        t_celldata.NodalPhyIDVector_Global[2]=10003;
        t_celldata.NodalPhyIDVector_Global[3]=10004;

        t_celldata.NodalPhyNameVector_Global[0]="leftnodes";
        t_celldata.NodalPhyNameVector_Global[1]="rightnodes";
        t_celldata.NodalPhyNameVector_Global[2]="bottomnodes";
        t_celldata.NodalPhyNameVector_Global[3]="topnodes";
        
        t_celldata.NodalPhyID2NameMap_Global[10001]="leftnodes";
        t_celldata.NodalPhyID2NameMap_Global[10002]="rightnodes";
        t_celldata.NodalPhyID2NameMap_Global[10003]="bottomnodes";
        t_celldata.NodalPhyID2NameMap_Global[10004]="topnodes";
        
        t_celldata.NodalPhyName2IDMap_Global["leftnodes"]=10001;
        t_celldata.NodalPhyName2IDMap_Global["rightnodes"]=10002;
        t_celldata.NodalPhyName2IDMap_Global["bottomnodes"]=10003;
        t_celldata.NodalPhyName2IDMap_Global["topnodes"]=10004;
        
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
        
        //*** for nodal physical group
        // for leftnodes
        MPIDataBus::receivePhyName2NodeIDVecMapFromMaster(t_celldata.NodalPhyName2NodeIDVecMap_Local,1000*rank+300);
        // for rightnodes
        MPIDataBus::receivePhyName2NodeIDVecMapFromMaster(t_celldata.NodalPhyName2NodeIDVecMap_Local,1000*rank+320);
        // for bottomnodes
        MPIDataBus::receivePhyName2NodeIDVecMapFromMaster(t_celldata.NodalPhyName2NodeIDVecMap_Local,1000*rank+340);
        // for topnodes
        MPIDataBus::receivePhyName2NodeIDVecMapFromMaster(t_celldata.NodalPhyName2NodeIDVecMap_Local,1000*rank+360);
    }

    return true;
}