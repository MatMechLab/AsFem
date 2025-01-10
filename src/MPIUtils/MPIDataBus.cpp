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

void MPIDataBus::printRankStatus(){
    int rank,size;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    printf("\033[1;91m[Deeeebuuugggg]  =====> in rank-%5d of size=%5d\n",rank,size);
}

void MPIDataBus::sendIntegerToOthers(const int &Val,const int &Tag,const int &Cpuid){
    MPI_Request request;
    MPI_Isend(&Val,1,MPI_INT,Cpuid,Tag,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
}
void MPIDataBus::receiveIntegerFromMaster(int &Val,const int &Tag){
    MPI_Request request;
    MPI_Irecv(&Val,1,MPI_INT,0,Tag,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
}
//********************************************************
void MPIDataBus::sendStringToOthers(const string &Txt,const int &Tag,const int &Cpuid){
    MPI_Request request;
    int datasize;
    datasize=static_cast<int>(Txt.size());
    MPI_Isend(&datasize,1,MPI_INT,Cpuid,Tag,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    //
    MPI_Isend(Txt.data(),datasize,MPI_CHAR,Cpuid,Tag+1,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
}
void MPIDataBus::receiveStringFromMaster(string &Txt,const int &Tag){
    MPI_Request request;
    int datasize;
    MPI_Irecv(&datasize,1,MPI_INT,0,Tag,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    //
    Txt.resize(datasize);
    MPI_Irecv(Txt.data(),datasize,MPI_CHAR,0,Tag+1,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
}
//********************************************************
void MPIDataBus::sendMeshTypeToOthers(const MeshType &t_MeshType,const int &Tag,const int &Cpuid) {
    MPI_Request request;
    int val;
    val=static_cast<int>(t_MeshType);
    MPI_Isend(&val,1,MPI_INT,Cpuid,Tag,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
}
void MPIDataBus::receiveMeshTypeFromMater(MeshType &t_MeshType,const int &Tag) {
    MPI_Request request;
    int val;
    MPI_Irecv(&val,1,MPI_INT,0,Tag,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    t_MeshType=static_cast<MeshType>(val);
}
//********************************************************
void MPIDataBus::sentIntegerVectorToOthers(const vector<int> &Vec,const int &Tag,const int &Cpuid){
    MPI_Request request;
    int datasize;
    datasize=static_cast<int>(Vec.size());
    MPI_Isend(&datasize,1,MPI_INT,Cpuid,Tag,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    if (datasize>0) {
        MPI_Isend(Vec.data(),datasize,MPI_INT,Cpuid,Tag+1,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
}
void MPIDataBus::receiveIntegerVectorFromMaster(vector<int> &Vec,const int &Tag){
    MPI_Request request;
    int vecsize;
    MPI_Irecv(&vecsize,1,MPI_INT,0,Tag,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    Vec.resize(vecsize,0);
    if (vecsize>0) {
        MPI_Irecv(Vec.data(),vecsize,MPI_INT,0,Tag+1,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
}
//********************************************************
void MPIDataBus::sentVectorOfIntegerVectorToOthers(const vector<vector<int>> &Vec,const int &Tag,const int &Cpuid){
    MPI_Request request;
    int datasize,vecsize;
    vecsize=static_cast<int>(Vec.size());
    MPI_Isend(&vecsize,1,MPI_INT,Cpuid,Tag,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    for(int i=0;i<vecsize;i++){
        datasize=static_cast<int>(Vec[i].size());
        MPI_Isend(&datasize,1,MPI_INT,Cpuid,Tag+1,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        MPI_Isend(Vec[i].data(),datasize,MPI_INT,Cpuid,Tag+2,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
}
void MPIDataBus::receiveVectorOfIntegerVectorFromMaster(vector<vector<int>> &Vec,const int &Tag){
    MPI_Request request;
    int vecsize,datasize;
    MPI_Irecv(&vecsize,1,MPI_INT,0,Tag,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    Vec.resize(vecsize,vector<int>(0));
    for(int i=0;i<vecsize;i++){
        MPI_Irecv(&datasize,1,MPI_INT,0,Tag+1,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        Vec[i].resize(datasize,0);
        MPI_Irecv(Vec[i].data(),datasize,MPI_INT,0,Tag+2,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
}
//********************************************************
void MPIDataBus::sendAlignedVectorOfIntegerVectorToOthers(const vector<vector<int>> &Vec,const int &Tag,const int &Cpuid){
    /**
     * The input Vec must have the same number of cols, otherwise it dosen\'t work!!!
     */
    MPI_Request request;
    int vecsize,width;
    vector<int> tempvec;
    vecsize=static_cast<int>(Vec.size());
    MPI_Isend(&vecsize,1,MPI_INT,Cpuid,Tag,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    //
    if (vecsize>0) {
        width=static_cast<int>(Vec[0].size());
        tempvec.clear();
        for(int i=0;i<vecsize;i++){
            for(int j=0;j<width;j++){
                tempvec.push_back(Vec[i][j]);
            }
        }
        MPI_Isend(&width,1,MPI_INT,Cpuid,Tag+1,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        MPI_Isend(tempvec.data(),vecsize*width,MPI_INT,Cpuid,Tag+2,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
    tempvec=vector<int>(0);// release memory
}
void MPIDataBus::receiveAlignedVectorOfIntegerVectorFromMaster(vector<vector<int>> &Vec,const int &Tag){
    MPI_Request request;
    int vecsize,width;
    vector<int> tempvec;
    MPI_Irecv(&vecsize,1,MPI_INT,0,Tag,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    //
    if (vecsize>0) {
        MPI_Irecv(&width,1,MPI_INT,0,Tag+1,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        tempvec.resize(vecsize*width,0);
        MPI_Irecv(tempvec.data(),vecsize*width,MPI_INT,0,Tag+2,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        Vec.resize(vecsize,vector<int>(width,0));
        int k;
        k=0;
        for(int i=0;i<vecsize;i++){
            for(int j=0;j<width;j++){
                Vec[i][j]=tempvec[k];
                k+=1;
            }
        }
    }
    tempvec=vector<int>(0);
}
//********************************************************
void MPIDataBus::sentStringVectorToOthers(const vector<string> &Vec,const int &Tag,const int &Cpuid){
    MPI_Request request;
    int datasize,vecsize;
    vecsize=static_cast<int>(Vec.size());
    MPI_Isend(&vecsize,1,MPI_INT,Cpuid,Tag,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    for(int i=0;i<vecsize;i++){
        datasize=static_cast<int>(Vec[i].size());
        MPI_Isend(&datasize,1,MPI_INT,Cpuid,Tag+1,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        MPI_Isend(Vec[i].data(),datasize,MPI_CHAR,Cpuid,Tag+2,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
}
void MPIDataBus::receiveStringVectorFromMaster(vector<string> &Vec,const int &Tag){
    MPI_Request request;
    int datasize,vecsize;
    MPI_Irecv(&vecsize,1,MPI_INT,0,Tag,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    Vec.resize(vecsize);
    for(int i=0;i<vecsize;i++){
        MPI_Irecv(&datasize,1,MPI_INT,0,Tag+1,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        Vec[i].resize(datasize);
        MPI_Irecv(Vec[i].data(),datasize,MPI_CHAR,0,Tag+2,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
}
//********************************************************
void MPIDataBus::sendMeshCellToOthers(const vector<SingleMeshCell> &Meshcellvec,const int &Tag,const int &Cpuid){
    int Basetag;
    Basetag=Tag;

    int datasize;

    MPI_Request request;

    // send out the length of the mesh cell vector
    datasize=static_cast<int>(Meshcellvec.size());
    MPI_Isend(&datasize,1,MPI_INT,Cpuid,Basetag+1,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);

    for(const auto &cell:Meshcellvec){
        // send Dim
        MPI_Isend(&cell.Dim,1,MPI_INT,Cpuid,Basetag+2,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // send NodesNumPerElmt
        MPI_Isend(&cell.NodesNumPerElmt,1,MPI_INT,Cpuid,Basetag+3,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // send ElmtConn
        datasize=static_cast<int>(cell.ElmtConn.size());
        MPI_Isend(&datasize,1,MPI_INT,Cpuid,Basetag+4,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        MPI_Isend(cell.ElmtConn.data(),datasize,MPI_INT,Cpuid,Basetag+5,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // send ElmtNodeCoords0
        datasize=static_cast<int>(cell.ElmtNodeCoords0.getSize());
        MPI_Isend(&datasize,1,MPI_INT,Cpuid,Basetag+6,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        MPI_Isend(cell.ElmtNodeCoords0.getCopy().data(),datasize*3,MPI_DOUBLE,Cpuid,Basetag+7,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // send VTKCellType
        MPI_Isend(&cell.VTKCellType,1,MPI_INT,Cpuid,Basetag+8,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }

}
//********************************************************
void MPIDataBus::receiveMeshCellFromMaster(vector<SingleMeshCell> &Meshcellvec,const int &Tag){
    int Basetag;
    Basetag=Tag;

    int datasize;

    MPI_Request request;

    // receive the length of the mesh cell vector
    MPI_Irecv(&datasize,1,MPI_INT,0,Basetag+1,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);

    Meshcellvec.resize(datasize);
    for(auto &cell:Meshcellvec){
        // receive Dim 
        MPI_Irecv(&cell.Dim,1,MPI_INT,0,Basetag+2,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // receive NodesNumPerElmt
        MPI_Irecv(&cell.NodesNumPerElmt,1,MPI_INT,0,Basetag+3,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // receive ElmtConn
        MPI_Irecv(&datasize,1,MPI_INT,0,Basetag+4,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        cell.ElmtConn.resize(datasize,0);
        MPI_Irecv(cell.ElmtConn.data(),datasize,MPI_INT,0,Basetag+5,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // receive ElmtNodeCoords0
        MPI_Irecv(&datasize,1,MPI_INT,0,Basetag+6,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        cell.ElmtNodeCoords0.resize(datasize);
        MPI_Irecv(cell.ElmtNodeCoords0.getData(),datasize*3,MPI_DOUBLE,0,Basetag+7,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        cell.ElmtNodeCoords=cell.ElmtNodeCoords0;// make a copy

        // receive VTKCellType
        MPI_Irecv(&cell.VTKCellType,1,MPI_INT,0,Basetag+8,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
}

//***********************************************************
void MPIDataBus::sendPhyName2MeshCellMapToOthers(const string &Phyname,const vector<SingleMeshCell> &Meshcellvec,const int &Tag,const int &Cpuid){
    int Basetag;
    Basetag=Tag;

    int datasize;

    MPI_Request request;

    // send out the physical name to each rank
    datasize=static_cast<int>(Phyname.size());
    MPI_Isend(&datasize,1,MPI_INT,Cpuid,Basetag+1,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    // send the char
    MPI_Isend(Phyname.c_str(),datasize,MPI_CHAR,Cpuid,Basetag+2,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);

    // send out the length of the mesh cell vector
    datasize=static_cast<int>(Meshcellvec.size());
    MPI_Isend(&datasize,1,MPI_INT,Cpuid,Basetag+3,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);

    for(const auto &cell:Meshcellvec){
        // send Dim
        MPI_Isend(&cell.Dim,1,MPI_INT,Cpuid,Basetag+4,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // send NodesNumPerElmt
        MPI_Isend(&cell.NodesNumPerElmt,1,MPI_INT,Cpuid,Basetag+5,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // send ElmtConn
        datasize=static_cast<int>(cell.ElmtConn.size());
        MPI_Isend(&datasize,1,MPI_INT,Cpuid,Basetag+6,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        MPI_Isend(cell.ElmtConn.data(),datasize,MPI_INT,Cpuid,Basetag+7,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // send ElmtNodeCoords0
        datasize=static_cast<int>(cell.ElmtNodeCoords0.getSize());
        MPI_Isend(&datasize,1,MPI_INT,Cpuid,Basetag+8,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        MPI_Isend(cell.ElmtNodeCoords0.getCopy().data(),datasize*3,MPI_DOUBLE,Cpuid,Basetag+9,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // send VTKCellType
        MPI_Isend(&cell.VTKCellType,1,MPI_INT,Cpuid,Basetag+10,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
}
//*************************************************************************
void MPIDataBus::receivePhyName2MeshCellMapFromMaster(map<string,vector<SingleMeshCell>> &Localmap,const int &Tag){
    int Basetag;
    Basetag=Tag;

    string Phyname;
    char buff[108];

    int datasize;

    MPI_Request request;
    vector<SingleMeshCell> Meshcellvec;

    // receive the size of physical name 
    MPI_Irecv(&datasize,1,MPI_INT,0,Basetag+1,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    // receive the char of phy name
    MPI_Irecv(buff,datasize,MPI_CHAR,0,Basetag+2,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    Phyname.clear();
    for(int i=0;i<datasize;i++) Phyname.push_back(buff[i]);

    // receive the length of the mesh cell vector
    MPI_Irecv(&datasize,1,MPI_INT,0,Basetag+3,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);

    Meshcellvec.resize(datasize);
    for(auto &cell:Meshcellvec){
        // receive Dim
        MPI_Irecv(&cell.Dim,1,MPI_INT,0,Basetag+4,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // receive NodesNumPerElmt
        MPI_Irecv(&cell.NodesNumPerElmt,1,MPI_INT,0,Basetag+5,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // receive ElmtConn
        datasize=0;
        MPI_Irecv(&datasize,1,MPI_INT,0,Basetag+6,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        cell.ElmtConn.resize(datasize,0);
        MPI_Irecv(cell.ElmtConn.data(),datasize,MPI_INT,0,Basetag+7,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // receive ElmtNodeCoords0
        datasize=0;
        MPI_Irecv(&datasize,1,MPI_INT,0,Basetag+8,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        cell.ElmtNodeCoords0.resize(datasize);
        MPI_Irecv(cell.ElmtNodeCoords0.getData(),datasize*3,MPI_DOUBLE,0,Basetag+9,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        cell.ElmtNodeCoords=cell.ElmtNodeCoords0;

        // receive VTKCellType
        MPI_Irecv(&cell.VTKCellType,1,MPI_INT,0,Basetag+10,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }

    Localmap[Phyname]=Meshcellvec;
    Meshcellvec.clear();
}

//**************************************************************
void MPIDataBus::sendPhyID2MeshCellMapToOthers(const int &Phyid,const vector<SingleMeshCell> &Meshcellvec,const int &Tag,const int &Cpuid){
    int Basetag;
    Basetag=Tag;

    int datasize;

    MPI_Request request;

    // send out the physical name to each rank
    MPI_Isend(&Phyid,1,MPI_INT,Cpuid,Basetag+1,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);

    // send out the length of the mesh cell vector
    datasize=static_cast<int>(Meshcellvec.size());
    MPI_Isend(&datasize,1,MPI_INT,Cpuid,Basetag+2,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);


    for(const auto &cell:Meshcellvec){
        // send Dim
        MPI_Isend(&cell.Dim,1,MPI_INT,Cpuid,Basetag+3,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // send NodesNumPerElmt
        MPI_Isend(&cell.NodesNumPerElmt,1,MPI_INT,Cpuid,Basetag+4,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // send ElmtConn
        datasize=static_cast<int>(cell.ElmtConn.size());
        MPI_Isend(&datasize,1,MPI_INT,Cpuid,Basetag+5,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        MPI_Isend(cell.ElmtConn.data(),datasize,MPI_INT,Cpuid,Basetag+6,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // send ElmtNodeCoords0
        datasize=static_cast<int>(cell.ElmtNodeCoords0.getSize());
        MPI_Isend(&datasize,1,MPI_INT,Cpuid,Basetag+7,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        MPI_Isend(cell.ElmtNodeCoords0.getCopy().data(),datasize*3,MPI_DOUBLE,Cpuid,Basetag+8,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // send VTKCellType
        MPI_Isend(&cell.VTKCellType,1,MPI_INT,Cpuid,Basetag+9,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
}
//***************************************************************
void MPIDataBus::receivePhyID2MeshCellMapFromMaster(map<int,vector<SingleMeshCell>> &Localmap,const int &Tag){
    int Basetag;
    Basetag=Tag;

    int datasize;
    int phyid;

    MPI_Request request;
    vector<SingleMeshCell> Meshcellvec;

    // receive the physical name to each rank
    MPI_Irecv(&phyid,1,MPI_INT,0,Basetag+1,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);

    // receive the length of the mesh cell vector
    MPI_Irecv(&datasize,1,MPI_INT,0,Basetag+2,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);

    Meshcellvec.resize(datasize);
    for(auto &cell:Meshcellvec){
        // receive Dim
        MPI_Irecv(&cell.Dim,1,MPI_INT,0,Basetag+3,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // receive NodesNumPerElmt
        MPI_Irecv(&cell.NodesNumPerElmt,1,MPI_INT,0,Basetag+4,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // receive ElmtConn
        datasize=0;
        MPI_Irecv(&datasize,1,MPI_INT,0,Basetag+5,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        cell.ElmtConn.resize(datasize,0);
        MPI_Irecv(cell.ElmtConn.data(),datasize,MPI_INT,0,Basetag+6,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // receive ElmtNodeCoords0
        datasize=0;
        MPI_Irecv(&datasize,1,MPI_INT,0,Basetag+7,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        cell.ElmtNodeCoords0.resize(datasize);
        MPI_Irecv(cell.ElmtNodeCoords0.getData(),datasize*3,MPI_DOUBLE,0,Basetag+8,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        cell.ElmtNodeCoords=cell.ElmtNodeCoords0;

        // receive VTKCellType
        MPI_Irecv(&cell.VTKCellType,1,MPI_INT,0,Basetag+9,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }

    Localmap[phyid]=Meshcellvec;
    Meshcellvec.clear();
}

//******************************************************************
void MPIDataBus::sendPhyName2NodeIDVecMapToOthers(const string &Phyname,const vector<int> &Nodesid,const int &Tag,const int &Cpuid){
    int Basetag;
    Basetag=Tag;

    int datasize;

    MPI_Request request;

    // send out the physical name to each rank
    datasize=static_cast<int>(Phyname.size());
    MPI_Isend(&datasize,1,MPI_INT,Cpuid,Basetag+1,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    // send the char of the phy name
    MPI_Isend(Phyname.c_str(),datasize,MPI_CHAR,Cpuid,Basetag+2,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);

    // send out the length of the node id vector
    datasize=static_cast<int>(Nodesid.size());
    MPI_Isend(&datasize,1,MPI_INT,Cpuid,Basetag+3,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    if (datasize>0) {
        // send out the node id vector
        MPI_Isend(Nodesid.data(),datasize,MPI_INT,Cpuid,Basetag+4,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
}
//*********************************+
void MPIDataBus::receivePhyName2NodeIDVecMapFromMaster(map<string,vector<int>> &Localmap,const int &Tag){
    int Basetag;
    Basetag=Tag;

    string phyname;
    char buff[108];

    int datasize;

    MPI_Request request;
    vector<int> nodeidvec;

    // receive the physical name length from master rank
    MPI_Irecv(&datasize,1,MPI_INT,0,Basetag+1,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    // receive the phy name
    MPI_Irecv(buff,datasize,MPI_CHAR,0,Basetag+2,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    phyname.clear();
    for(int i=0;i<datasize;i++) phyname.push_back(buff[i]);

    // receive the length of the nodeid vector
    MPI_Irecv(&datasize,1,MPI_INT,0,Basetag+3,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    
    nodeidvec.resize(datasize,0);
    if (datasize>0) {
        // receive the nodeid vector
        MPI_Irecv(nodeidvec.data(),datasize,MPI_INT,0,Basetag+4,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }

    Localmap[phyname]=nodeidvec;
    nodeidvec.clear();
}