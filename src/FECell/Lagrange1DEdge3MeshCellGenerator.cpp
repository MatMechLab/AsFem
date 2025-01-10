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

#include "FECell/Lagrange1DEdge3MeshCellGenerator.h"
#include "FECell/SingleMeshCell.h"

Lagrange1DEdge3MeshCellGenerator::Lagrange1DEdge3MeshCellGenerator(){
    m_mesh_generated=false;
}
Lagrange1DEdge3MeshCellGenerator::~Lagrange1DEdge3MeshCellGenerator(){
    m_mesh_generated=false;
}
//*********************************************************
bool Lagrange1DEdge3MeshCellGenerator::generateFECell(FECellData &t_celldata){
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
        
        // generate the mesh cell for edge mesh
        t_celldata.MeshOrder=2;
        t_celldata.BulkMeshTypeName="edge3";

        t_celldata.BulkElmtsNum=t_celldata.Nx;
        t_celldata.LineElmtsNum=0;
        t_celldata.SurfElmtsNum=0;
        t_celldata.ElmtsNum=t_celldata.BulkElmtsNum
                           +t_celldata.SurfElmtsNum
                           +t_celldata.LineElmtsNum;
                               
        t_celldata.NodesNum=t_celldata.Nx*t_celldata.MeshOrder+1;
        t_celldata.NodesNumPerBulkElmt=3;
        t_celldata.NodesNumPerSurfElmt=0;
        t_celldata.NodesNumPerLineElmt=0;

        double dx,dy,dz;
        dx=(t_celldata.Xmax-t_celldata.Xmin)/(t_celldata.NodesNum-1);
        dy=(t_celldata.Ymax-t_celldata.Ymin)/(t_celldata.NodesNum-1);
        dz=(t_celldata.Zmax-t_celldata.Zmin)/(t_celldata.NodesNum-1);
            
        t_celldata.BulkElmtVTKCellType=4;
        t_celldata.BulkElmtMeshType=MeshType::EDGE3;

        t_celldata.LineElmtVTKCellType=0;
        t_celldata.LineElmtMeshType=MeshType::NULLTYPE;
            
        t_celldata.SurfElmtVTKCellType=0;
        t_celldata.SurfElmtMeshType=MeshType::NULLTYPE;

        leftnodes.clear();rightnodes.clear();
        t_celldata.NodeCoords_Global.resize(t_celldata.NodesNum*3,0.0);
        for(int i=0;i<t_celldata.NodesNum;i++){
            t_celldata.NodeCoords_Global[i*3+1-1]=t_celldata.Xmin+i*dx;
            t_celldata.NodeCoords_Global[i*3+2-1]=t_celldata.Ymin+i*dy;
            t_celldata.NodeCoords_Global[i*3+3-1]=t_celldata.Zmin+i*dz;
        }// end-of-node-generation

        // for the connectivity information of bulk elements
        t_celldata.MeshCell_Total.resize(t_celldata.BulkElmtsNum);
        leftconn.resize(1);
        rightconn.resize(1);
            
        leftnodes.clear();rightnodes.clear();

        t_celldata.PhyID2BulkFECellIDMap_Global.clear();
        t_celldata.PhyName2BulkFECellIDMap_Global.clear();

        int i1,i2,i3;
        for(int e=1;e<=t_celldata.BulkElmtsNum;e++){
            t_celldata.PhyName2BulkFECellIDMap_Global["alldomain"].push_back(e);
            t_celldata.PhyID2BulkFECellIDMap_Global[0].push_back(e);

            i1=(e-1)*t_celldata.MeshOrder+1;
            i2=i1+1;
            i3=i2+1;

            t_celldata.MeshCell_Total[e-1].Dim=1;
            t_celldata.MeshCell_Total[e-1].NodesNumPerElmt=3;
            t_celldata.MeshCell_Total[e-1].VTKCellType=4;
            t_celldata.MeshCell_Total[e-1].ElmtConn.clear();
            t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i1);
            t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i2);
            t_celldata.MeshCell_Total[e-1].ElmtConn.push_back(i3);

            t_celldata.MeshCell_Total[e-1].ElmtNodeCoords.resize(3);
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

            t_celldata.MeshCell_Total[e-1].ElmtNodeCoords0=t_celldata.MeshCell_Total[e-1].ElmtNodeCoords;

               
            if(i1==1){
                // for left node
                leftconn[1-1].Dim=0;
                leftconn[1-1].NodesNumPerElmt=1;
                leftconn[1-1].VTKCellType=1;

                leftconn[1-1].ElmtConn.clear();
                leftconn[1-1].ElmtConn.push_back(i1);

                leftconn[1-1].ElmtNodeCoords0.resize(1);
                leftconn[1-1].ElmtNodeCoords0(1,1)=t_celldata.NodeCoords_Global[(i1-1)*3+1-1];
                leftconn[1-1].ElmtNodeCoords0(1,2)=t_celldata.NodeCoords_Global[(i1-1)*3+2-1];
                leftconn[1-1].ElmtNodeCoords0(1,3)=t_celldata.NodeCoords_Global[(i1-1)*3+3-1];
                leftconn[1-1].ElmtNodeCoords=leftconn[1-1].ElmtNodeCoords0;

                leftnodes.push_back(i1);
            }
            if(i3==t_celldata.NodesNum) {
                // for left node
                rightconn[1-1].Dim=0;
                rightconn[1-1].NodesNumPerElmt=1;
                rightconn[1-1].VTKCellType=1;

                rightconn[1-1].ElmtConn.clear();
                rightconn[1-1].ElmtConn.push_back(i3);

                rightconn[1-1].ElmtNodeCoords0.resize(1);
                rightconn[1-1].ElmtNodeCoords0(1,1)=t_celldata.NodeCoords_Global[(i3-1)*3+1-1];
                rightconn[1-1].ElmtNodeCoords0(1,2)=t_celldata.NodeCoords_Global[(i3-1)*3+2-1];
                rightconn[1-1].ElmtNodeCoords0(1,3)=t_celldata.NodeCoords_Global[(i3-1)*3+3-1];
                rightconn[1-1].ElmtNodeCoords=rightconn[1-1].ElmtNodeCoords0;

                rightnodes.push_back(i3);
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

        // for mesh cell partition info
        t_celldata.BulkCellPartionInfo_Global.resize(t_celldata.BulkElmtsNum,0);

    }// end-of-if(rank==0)
    else{
        t_celldata.MeshOrder=2;
        t_celldata.BulkMeshTypeName="edge3";

        t_celldata.BulkElmtsNum=t_celldata.Nx;
        t_celldata.LineElmtsNum=0;
        t_celldata.SurfElmtsNum=0;
        t_celldata.ElmtsNum=t_celldata.BulkElmtsNum
                           +t_celldata.SurfElmtsNum
                           +t_celldata.LineElmtsNum;

        t_celldata.NodesNum=t_celldata.Nx*t_celldata.MeshOrder+1;
        t_celldata.NodesNumPerBulkElmt=3;
        t_celldata.NodesNumPerSurfElmt=0;
        t_celldata.NodesNumPerLineElmt=0;
            
        t_celldata.BulkElmtVTKCellType=4;
        t_celldata.BulkElmtMeshType=MeshType::EDGE3;

        t_celldata.LineElmtVTKCellType=0;
        t_celldata.LineElmtMeshType=MeshType::NULLTYPE;
            
        t_celldata.SurfElmtVTKCellType=0;
        t_celldata.SurfElmtMeshType=MeshType::NULLTYPE;
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
    }

    return true;
}