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
//+++ Date   : 2024.01.16
//+++ Purpose: Implement the function of send/receive mesh cell
//+++          among mpi ranks
//+++          Send   : master ----> other ranks
//+++          Receive: others <----- master
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MPIUtils/MPIDataBus.h"

void MPIDataBus::sendMeshCellToOthers(const vector<SingleMeshCell> &meshcellvec,const int &tag,const int &cpuid){
    int basetag;
    basetag=tag;

    int datasize;

    MPI_Request request;

    // send out the length of the mesh cell vector
    datasize=static_cast<int>(meshcellvec.size());
    MPI_Isend(&datasize,1,MPI_INT,cpuid,basetag+1,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);

    for(const auto &cell:meshcellvec){
        // send Dim
        MPI_Isend(&cell.Dim,1,MPI_INT,cpuid,basetag+2,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // send NodesNumPerElmt
        MPI_Isend(&cell.NodesNumPerElmt,1,MPI_INT,cpuid,basetag+3,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // send ElmtConn
        datasize=static_cast<int>(cell.ElmtConn.size());
        MPI_Isend(&datasize,1,MPI_INT,cpuid,basetag+4,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        MPI_Isend(cell.ElmtConn.data(),datasize,MPI_INT,cpuid,basetag+5,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // send ElmtNodeCoords0
        datasize=static_cast<int>(cell.ElmtNodeCoords0.getSize());
        MPI_Isend(&datasize,1,MPI_INT,cpuid,basetag+6,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        MPI_Isend(cell.ElmtNodeCoords0.getCopy().data(),datasize*3,MPI_INT,cpuid,basetag+7,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // send VTKCellType
        MPI_Isend(&cell.VTKCellType,1,MPI_INT,cpuid,basetag+8,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }

}
//********************************************************
void MPIDataBus::receiveMeshCellFromMaster(vector<SingleMeshCell> &meshcellvec,const int &tag){
    int basetag;
    basetag=tag;

    int datasize;

    MPI_Request request;

    // receive the length of the mesh cell vector
    MPI_Irecv(&datasize,1,MPI_INT,0,basetag+1,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);

    meshcellvec.resize(datasize);
    for(auto &cell:meshcellvec){
        // receive Dim 
        MPI_Irecv(&cell.Dim,1,MPI_INT,0,basetag+2,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // receive NodesNumPerElmt
        MPI_Irecv(&cell.NodesNumPerElmt,1,MPI_INT,0,basetag+3,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // receive ElmtConn
        MPI_Irecv(&datasize,1,MPI_INT,0,basetag+4,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        cell.ElmtConn.resize(datasize,0);
        MPI_Irecv(cell.ElmtConn.data(),datasize,MPI_INT,0,basetag+5,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // receive ElmtNodeCoords0
        MPI_Irecv(&datasize,1,MPI_INT,0,basetag+6,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        cell.ElmtNodeCoords0.resize(datasize);
        MPI_Irecv(cell.ElmtNodeCoords0.getData(),datasize*3,MPI_INT,0,basetag+7,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        cell.ElmtNodeCoords=cell.ElmtNodeCoords0;// make a copy

        // receive VTKCellType
        MPI_Irecv(&cell.VTKCellType,1,MPI_INT,0,basetag+8,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
}

//***********************************************************
void MPIDataBus::sendPhyName2MeshCellMapToOthers(const string &phyname,const vector<SingleMeshCell> &meshcellvec,const int &tag,const int &cpuid){
    int basetag;
    basetag=tag;

    int datasize;

    MPI_Request request;

    // send out the physical name to each rank
    datasize=static_cast<int>(phyname.size());
    MPI_Isend(&datasize,1,MPI_INT,cpuid,basetag+1,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    // send the char
    MPI_Isend(phyname.c_str(),datasize,MPI_CHAR,cpuid,basetag+2,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);

    // send out the length of the mesh cell vector
    datasize=static_cast<int>(meshcellvec.size());
    MPI_Isend(&datasize,1,MPI_INT,cpuid,basetag+3,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);

    for(const auto &cell:meshcellvec){
        // send Dim
        MPI_Isend(&cell.Dim,1,MPI_INT,cpuid,basetag+4,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // send NodesNumPerElmt
        MPI_Isend(&cell.NodesNumPerElmt,1,MPI_INT,cpuid,basetag+5,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // send ElmtConn
        datasize=static_cast<int>(cell.ElmtConn.size());
        MPI_Isend(&datasize,1,MPI_INT,cpuid,basetag+6,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        MPI_Isend(cell.ElmtConn.data(),datasize,MPI_INT,cpuid,basetag+7,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // send ElmtNodeCoords0
        datasize=static_cast<int>(cell.ElmtNodeCoords0.getSize());
        MPI_Isend(&datasize,1,MPI_INT,cpuid,basetag+8,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        MPI_Isend(cell.ElmtNodeCoords0.getCopy().data(),datasize*3,MPI_INT,cpuid,basetag+9,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // send VTKCellType
        MPI_Isend(&cell.VTKCellType,1,MPI_INT,cpuid,basetag+10,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
}
//*************************************************************************
void MPIDataBus::receivePhyName2MeshCellMapFromMaster(map<string,vector<SingleMeshCell>> &localmap,const int &tag){
    int basetag;
    basetag=tag;

    string phyname;
    char buff[108];

    int datasize;

    MPI_Request request;
    vector<SingleMeshCell> meshcellvec;

    // receive the physical name 
    MPI_Irecv(&datasize,1,MPI_INT,0,basetag+1,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    // receive the char of phy name
    MPI_Irecv(buff,datasize,MPI_CHAR,0,basetag+2,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    phyname.clear();
    for(int i=0;i<datasize;i++) phyname.push_back(buff[i]);

    // receive the length of the mesh cell vector
    MPI_Irecv(&datasize,1,MPI_INT,0,basetag+3,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);

    meshcellvec.resize(datasize);
    for(auto &cell:meshcellvec){
        // receive Dim
        MPI_Irecv(&cell.Dim,1,MPI_INT,0,basetag+4,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // receive NodesNumPerElmt
        MPI_Irecv(&cell.NodesNumPerElmt,1,MPI_INT,0,basetag+5,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // receive ElmtConn
        datasize=0;
        MPI_Irecv(&datasize,1,MPI_INT,0,basetag+6,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        cell.ElmtConn.resize(datasize,0);
        MPI_Irecv(cell.ElmtConn.data(),datasize,MPI_INT,0,basetag+7,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // receive ElmtNodeCoords0
        datasize=0;
        MPI_Irecv(&datasize,1,MPI_INT,0,basetag+8,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        cell.ElmtNodeCoords0.resize(datasize);
        MPI_Irecv(cell.ElmtNodeCoords0.getData(),datasize*3,MPI_INT,0,basetag+9,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        cell.ElmtNodeCoords=cell.ElmtNodeCoords0;

        // receive VTKCellType
        MPI_Irecv(&cell.VTKCellType,1,MPI_INT,0,basetag+10,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }

    localmap[phyname]=meshcellvec;
    meshcellvec.clear();
}

//**************************************************************
void MPIDataBus::sendPhyID2MeshCellMapToOthers(const int &phyid,const vector<SingleMeshCell> &meshcellvec,const int &tag,const int &cpuid){
    int basetag;
    basetag=tag;

    int datasize;

    MPI_Request request;

    // send out the physical name to each rank
    MPI_Isend(&phyid,1,MPI_INT,cpuid,basetag+1,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);

    // send out the length of the mesh cell vector
    datasize=static_cast<int>(meshcellvec.size());
    MPI_Isend(&datasize,1,MPI_INT,cpuid,basetag+2,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);

    for(const auto &cell:meshcellvec){
        // send Dim
        MPI_Isend(&cell.Dim,1,MPI_INT,cpuid,basetag+3,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // send NodesNumPerElmt
        MPI_Isend(&cell.NodesNumPerElmt,1,MPI_INT,cpuid,basetag+4,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // send ElmtConn
        datasize=static_cast<int>(cell.ElmtConn.size());
        MPI_Isend(&datasize,1,MPI_INT,cpuid,basetag+5,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        MPI_Isend(cell.ElmtConn.data(),datasize,MPI_INT,cpuid,basetag+6,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // send ElmtNodeCoords0
        datasize=static_cast<int>(cell.ElmtNodeCoords0.getSize());
        MPI_Isend(&datasize,1,MPI_INT,cpuid,basetag+7,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        MPI_Isend(cell.ElmtNodeCoords0.getCopy().data(),datasize*3,MPI_INT,cpuid,basetag+8,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // send VTKCellType
        MPI_Isend(&cell.VTKCellType,1,MPI_INT,cpuid,basetag+9,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
}
//***************************************************************
void MPIDataBus::receivePhyID2MeshCellMapFromMaster(map<int,vector<SingleMeshCell>> &localmap,const int &tag){
    int basetag;
    basetag=tag;

    int datasize;
    int phyid;

    MPI_Request request;
    vector<SingleMeshCell> meshcellvec;

    // receive the physical name to each rank
    MPI_Irecv(&phyid,1,MPI_INT,0,basetag+1,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);

    // receive the length of the mesh cell vector
    MPI_Irecv(&datasize,1,MPI_INT,0,basetag+2,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);

    meshcellvec.resize(datasize);
    for(auto &cell:meshcellvec){
        // receive Dim
        MPI_Irecv(&cell.Dim,1,MPI_INT,0,basetag+3,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // receive NodesNumPerElmt
        MPI_Irecv(&cell.NodesNumPerElmt,1,MPI_INT,0,basetag+4,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // receive ElmtConn
        datasize=0;
        MPI_Irecv(&datasize,1,MPI_INT,0,basetag+5,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        cell.ElmtConn.resize(datasize,0);
        MPI_Irecv(cell.ElmtConn.data(),datasize,MPI_INT,0,basetag+6,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // receive ElmtNodeCoords0
        datasize=0;
        MPI_Irecv(&datasize,1,MPI_INT,0,basetag+7,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        cell.ElmtNodeCoords0.resize(datasize);
        MPI_Irecv(cell.ElmtNodeCoords0.getData(),datasize*3,MPI_INT,0,basetag+8,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        cell.ElmtNodeCoords=cell.ElmtNodeCoords0;

        // receive VTKCellType
        MPI_Irecv(&cell.VTKCellType,1,MPI_INT,0,basetag+9,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }

    localmap[phyid]=meshcellvec;
    meshcellvec.clear();
}

//******************************************************************
void MPIDataBus::sendPhyName2NodeIDVecMapToOthers(const string &phyname,const vector<int> &nodesid,const int &tag,const int &cpuid){
    int basetag;
    basetag=tag;

    int datasize;

    MPI_Request request;

    // send out the physical name to each rank
    datasize=static_cast<int>(phyname.size());
    MPI_Isend(&datasize,1,MPI_INT,cpuid,basetag+1,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    // send the char of the phy name
    MPI_Isend(phyname.c_str(),datasize,MPI_CHAR,cpuid,basetag+2,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);

    // send out the length of the node id vector
    datasize=static_cast<int>(nodesid.size());
    MPI_Isend(&datasize,1,MPI_INT,cpuid,basetag+3,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);

    // send out the node id vector
    MPI_Isend(nodesid.data(),datasize,MPI_INT,cpuid,basetag+4,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
}
//*********************************+
void MPIDataBus::receivePhyName2NodeIDVecMapFromMaster(map<string,vector<int>> &localmap,const int &tag){
    int basetag;
    basetag=tag;

    string phyname;
    char buff[108];

    int datasize;

    MPI_Request request;
    vector<int> nodeidvec;

    // receive the physical name length from master rank
    MPI_Irecv(&datasize,1,MPI_INT,0,basetag+1,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    // receive the phy name
    MPI_Irecv(buff,datasize,MPI_CHAR,0,basetag+2,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    phyname.clear();
    for(int i=0;i<datasize;i++) phyname.push_back(buff[i]);

    // receive the length of the nodeid vector
    MPI_Irecv(&datasize,1,MPI_INT,0,basetag+3,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    
    nodeidvec.resize(datasize,0);
    // receive the nodeid vector
    MPI_Irecv(nodeidvec.data(),datasize,MPI_INT,0,basetag+4,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);

    localmap[phyname]=nodeidvec;
    nodeidvec.clear();
}