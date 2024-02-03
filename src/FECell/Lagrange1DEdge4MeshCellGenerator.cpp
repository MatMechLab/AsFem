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
//+++ Date   : 2024.02.03
//+++ Purpose: the 1D edge3 lagrange FE cell generator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FECell/Lagrange1DEdge4MeshCellGenerator.h"
#include "FECell/SingleMeshCell.h"
#include "MPIUtils/MPIDataBus.h"

Lagrange1DEdge4MeshCellGenerator::Lagrange1DEdge4MeshCellGenerator(){
    m_mesh_generated=false;
}
Lagrange1DEdge4MeshCellGenerator::~Lagrange1DEdge4MeshCellGenerator(){
    m_mesh_generated=false;
}
//*********************************************************
bool Lagrange1DEdge4MeshCellGenerator::generateFECell(FECellData &t_celldata){
    int rank,size;
    t_celldata.MaxDim=1;
    t_celldata.MinDim=0;

    t_celldata.ActiveDofsNum=0;
    t_celldata.TotalDofsNum=0;
    t_celldata.MaxDofsPerNode=0;

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    if(rank==0){
        vector<SingleMeshCell> leftconn,rightconn;
        vector<int> leftnodes,rightnodes;

        // only create mesh on the master rank, then distributed them into different ranks !
        double dx;
        int i,j,e;
        int i1,i2,i3,i4;
        
        // generate the mesh cell for edge mesh
        t_celldata.MeshOrder=3;
        t_celldata.BulkMeshTypeName="edge4";

        t_celldata.BulkElmtsNum=t_celldata.Nx;
        t_celldata.LineElmtsNum=0;
        t_celldata.SurfElmtsNum=0;
        t_celldata.ElmtsNum=t_celldata.BulkElmtsNum
                           +t_celldata.SurfElmtsNum
                           +t_celldata.LineElmtsNum;
                               
        t_celldata.NodesNum=t_celldata.Nx*t_celldata.MeshOrder+1;
        t_celldata.NodesNumPerBulkElmt=4;
        t_celldata.NodesNumPerSurfElmt=0;
        t_celldata.NodesNumPerLineElmt=0;

        dx=(t_celldata.Xmax-t_celldata.Xmin)/(t_celldata.NodesNum-1);
            
        t_celldata.BulkElmtVTKCellType=4;
        t_celldata.BulkElmtMeshType=MeshType::EDGE4;

        t_celldata.LineElmtVTKCellType=0;
        t_celldata.LineElmtMeshType=MeshType::NULLTYPE;
            
        t_celldata.SurfElmtVTKCellType=0;
        t_celldata.SurfElmtMeshType=MeshType::NULLTYPE;

        vector<double> nodecoords;
        nodecoords.resize(t_celldata.NodesNum*3,0.0);
        leftnodes.clear();rightnodes.clear();
        for(i=1;i<=t_celldata.NodesNum;i++){
            nodecoords[(i-1)*3+1-1]=t_celldata.Xmin+(i-1)*dx;
            nodecoords[(i-1)*3+2-1]=0.0;
            nodecoords[(i-1)*3+3-1]=0.0;
        }// end-of-node-generation

        t_celldata.NodeCoords_Global.clear();
        for(const auto &it:nodecoords) t_celldata.NodeCoords_Global.push_back(it);

        // for the connectivity information of bulk elements
        t_celldata.MeshCell_Total.resize(t_celldata.BulkElmtsNum);
        leftconn.resize(1);
        rightconn.resize(1);
            
        leftnodes.clear();rightnodes.clear();

        for(e=1;e<=t_celldata.BulkElmtsNum;e++){
            for(j=1;j<=t_celldata.NodesNumPerBulkElmt;j++){
                i1=(e-1)*t_celldata.MeshOrder+1;
                i2=i1+1;
                i3=i2+1;
                i4=i3+1;

                t_celldata.MeshCell_Total[e-1].Dim=1;
                t_celldata.MeshCell_Local[e-1].NodesNumPerElmt=4;
                t_celldata.MeshCell_Total[e-1].VTKCellType=4;
                t_celldata.MeshCell_Total[e-1].ElmtConn.clear();
                t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i1);
                t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i2);
                t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i3);
                t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i4);

                t_celldata.MeshCell_Total[e-1].ElmtNodeCoords.resize(3);
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

                t_celldata.MeshCell_Total[e-1].ElmtNodeCoords0=t_celldata.MeshCell_Total[e-1].ElmtNodeCoords;

               
                if(i1==1){
                    // for left node
                    leftconn[1-1].Dim=0;
                    leftconn[1-1].NodesNumPerElmt=1;
                    leftconn[1-1].VTKCellType=1;

                    leftconn[1-1].ElmtConn.clear();
                    leftconn[1-1].ElmtConn.push_back(i1);

                    leftnodes.push_back(i1);
                }
                if(i4==t_celldata.NodesNum){
                    // for left node
                    rightconn[1-1].Dim=0;
                    rightconn[1-1].NodesNumPerElmt=1;
                    rightconn[1-1].VTKCellType=1;

                    rightconn[1-1].ElmtConn.clear();
                    rightconn[1-1].ElmtConn.push_back(i4);

                    rightnodes.push_back(i4);
                }
            }
        }// end-of-element-generation-loop
    
        // remove the duplicate node id, for nodal type set, you don't need these duplicate node ids!
        sort(leftnodes.begin(),leftnodes.end());
        leftnodes.erase(unique(leftnodes.begin(),leftnodes.end()),leftnodes.end());
        // for rightnodes
        sort(rightnodes.begin(),rightnodes.end());
        rightnodes.erase(unique(rightnodes.begin(),rightnodes.end()),rightnodes.end());

        // setup the physical group information
        t_celldata.PhyGroupNum_Global=1+2;
        t_celldata.PhyDimVector_Global.resize(1+2,0);
        t_celldata.PhyIDVector_Global.resize(1+2,0);
        t_celldata.PhyNameVector_Global.resize(1+2);
        t_celldata.PhyGroupElmtsNumVector_Global.resize(1+2,0);

        // for physical dim vector
        t_celldata.PhyDimVector_Global[0]=1;
        t_celldata.PhyDimVector_Global[1]=0;
        t_celldata.PhyDimVector_Global[2]=0;

        // for physical id vector
        t_celldata.PhyIDVector_Global[0]=0;
        t_celldata.PhyIDVector_Global[1]=1;
        t_celldata.PhyIDVector_Global[2]=2;

        // for physical name vector
        t_celldata.PhyNameVector_Global[0]="alldomain";
        t_celldata.PhyNameVector_Global[1]="left";
        t_celldata.PhyNameVector_Global[2]="right";
       
        // for phy group elmts num vector
        t_celldata.PhyGroupElmtsNumVector_Global[0]=static_cast<int>(t_celldata.BulkElmtsNum);
        t_celldata.PhyGroupElmtsNumVector_Global[1]=static_cast<int>(leftconn.size());
        t_celldata.PhyGroupElmtsNumVector_Global[2]=static_cast<int>(rightconn.size());

        /**
         * setup id<---->name map
        */
        // id--->name map
        t_celldata.PhyID2NameMap_Global[0]="alldomain";
        t_celldata.PhyID2NameMap_Global[1]="left";
        t_celldata.PhyID2NameMap_Global[2]="right";
        // name--->id map
        t_celldata.PhyName2IDMap_Global["alldomain"]=0;
        t_celldata.PhyName2IDMap_Global["left"]=1;
        t_celldata.PhyName2IDMap_Global["right"]=2;

        /**
         * Setup the global mapping, this should only be nonzero on master rank !!!
        */
        t_celldata.PhyName2MeshCellVectorMap_Global["alldomain"]=t_celldata.MeshCell_Total;
        t_celldata.PhyName2MeshCellVectorMap_Global["left"]=leftconn;
        t_celldata.PhyName2MeshCellVectorMap_Global["right"]=rightconn;

        t_celldata.PhyID2MeshCellVectorMap_Global[0]=t_celldata.MeshCell_Total;
        t_celldata.PhyID2MeshCellVectorMap_Global[1]=leftconn;
        t_celldata.PhyID2MeshCellVectorMap_Global[2]=rightconn;

        /**
         * Setup the nodal physical info group
        */
        t_celldata.NodalPhyGroupNum_Global=2;
        t_celldata.NodalPhyIDVector_Global.resize(2);
        t_celldata.NodalPhyNameVector_Global.resize(2);
        t_celldata.NodalPhyGroupNodesNumVector_Global.resize(2);

        t_celldata.NodalPhyIDVector_Global[0]=10001;
        t_celldata.NodalPhyIDVector_Global[1]=10002;

        t_celldata.NodalPhyNameVector_Global[0]="leftnodes";
        t_celldata.NodalPhyNameVector_Global[1]="rightnodes";

        t_celldata.NodalPhyGroupNodesNumVector_Global[0]=static_cast<int>(leftnodes.size());
        t_celldata.NodalPhyGroupNodesNumVector_Global[1]=static_cast<int>(rightnodes.size());
        
        t_celldata.NodalPhyID2NameMap_Global[10001]="leftnodes";
        t_celldata.NodalPhyID2NameMap_Global[10002]="rightnodes";
        
        t_celldata.NodalPhyName2IDMap_Global["leftnodes"]=10001;
        t_celldata.NodalPhyName2IDMap_Global["rightnodes"]=10002;

        t_celldata.NodalPhyName2NodeIDVecMap_Global["leftnodes"]=leftnodes;
        t_celldata.NodalPhyName2NodeIDVecMap_Global["rightnodes"]=rightnodes;

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
                MPIDataBus::sendPhyName2NodeIDVecMapToOthers("leftnodes",nodeids,1000*cpuid+140,cpuid);
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
                MPIDataBus::sendPhyName2NodeIDVecMapToOthers("rightnodes",nodeids,1000*cpuid+160,cpuid);
            }
        }// end-of-cpuid-loop

    }// end-of-if(rank==0)
    else{
        // now we distribute the global mesh into different ranks
        // setup the physical group information
        t_celldata.PhyGroupNum_Global=1+2;
        t_celldata.PhyDimVector_Global.resize(1+2,0);
        t_celldata.PhyIDVector_Global.resize(1+2,0);
        t_celldata.PhyNameVector_Global.resize(1+2);
        t_celldata.PhyGroupElmtsNumVector_Global.resize(1+2,0);

        // for physical dim vector
        t_celldata.PhyDimVector_Global[0]=1;
        t_celldata.PhyDimVector_Global[1]=0;
        t_celldata.PhyDimVector_Global[2]=0;

        // for physical id vector
        t_celldata.PhyIDVector_Global[0]=0;
        t_celldata.PhyIDVector_Global[1]=1;
        t_celldata.PhyIDVector_Global[2]=2;

        // for physical name vector
        t_celldata.PhyNameVector_Global[0]="alldomain";
        t_celldata.PhyNameVector_Global[1]="left";
        t_celldata.PhyNameVector_Global[2]="right";

        /**
         * setup id<---->name map
        */
        // id--->name map
        t_celldata.PhyID2NameMap_Global[0]="alldomain";
        t_celldata.PhyID2NameMap_Global[1]="left";
        t_celldata.PhyID2NameMap_Global[2]="right";
        // name--->id map
        t_celldata.PhyName2IDMap_Global["alldomain"]=0;
        t_celldata.PhyName2IDMap_Global["left"]=1;
        t_celldata.PhyName2IDMap_Global["right"]=2;

        /**
         * Setup the nodal physical info group
        */
        t_celldata.NodalPhyGroupNum_Global=2;
        t_celldata.NodalPhyIDVector_Global.resize(2);
        t_celldata.NodalPhyNameVector_Global.resize(2);
        t_celldata.NodalPhyGroupNodesNumVector_Global.resize(2);

        t_celldata.NodalPhyIDVector_Global[0]=10001;
        t_celldata.NodalPhyIDVector_Global[1]=10002;

        t_celldata.NodalPhyNameVector_Global[0]="leftnodes";
        t_celldata.NodalPhyNameVector_Global[1]="rightnodes";
        
        t_celldata.NodalPhyID2NameMap_Global[10001]="leftnodes";
        t_celldata.NodalPhyID2NameMap_Global[10002]="rightnodes";
        
        t_celldata.NodalPhyName2IDMap_Global["leftnodes"]=10001;
        t_celldata.NodalPhyName2IDMap_Global["rightnodes"]=10002;
        
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
        
        //*** for nodal physical group
        // for leftnodes
        MPIDataBus::receivePhyName2NodeIDVecMapFromMaster(t_celldata.NodalPhyName2NodeIDVecMap_Local,1000*rank+140);
        // for rightnodes
        MPIDataBus::receivePhyName2NodeIDVecMapFromMaster(t_celldata.NodalPhyName2NodeIDVecMap_Local,1000*rank+160);
    }

    return true;
}