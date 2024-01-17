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

void MPIDataBus::sendMeshCell2Others(const vector<SingleMeshCell> &meshcellvec,const int &tag,const int &cpuid){
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
        MPI_Isend(cell.ElmtNodeCoords0.getCopy().data(),datasize,MPI_INT,cpuid,basetag+7,MPI_COMM_WORLD,&request);
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

    // send out the length of the mesh cell vector
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

        // send ElmtNodeCoords0
        MPI_Irecv(&datasize,1,MPI_INT,0,basetag+6,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        cell.ElmtNodeCoords0.resize(datasize);
        MPI_Irecv(cell.ElmtNodeCoords0.getData(),datasize,MPI_INT,0,basetag+7,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        cell.ElmtNodeCoords=cell.ElmtNodeCoords0;// make a copy

        // send VTKCellType
        MPI_Irecv(&cell.VTKCellType,1,MPI_INT,0,basetag+8,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
}