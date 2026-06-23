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

namespace{
    /**
     * Pack a vector of mesh cells into three flat buffers so that a whole cell
     * vector can be transferred with a handful of MPI messages instead of one
     * message per field per cell. Packing into named buffers also keeps the
     * send buffers alive until MPI_Wait completes.
     * @param cells  the mesh cell vector to be packed
     * @param header per-cell metadata, 5 ints per cell: Dim, NodesNumPerElmt, VTKCellType, ElmtConn size, nodes number
     * @param conn   the concatenation of every cell's ElmtConn vector
     * @param coords the concatenation of every cell's undeformed nodal coordinates (3 doubles per node)
     */
    void packMeshCells(const vector<SingleMeshCell> &cells,
                       vector<int> &header,vector<int> &conn,vector<double> &coords){
        header.clear();conn.clear();coords.clear();
        size_t conntotal=0,coordtotal=0;
        for(const auto &cell:cells){
            conntotal +=cell.ElmtConn.size();
            coordtotal+=static_cast<size_t>(cell.ElmtNodeCoords0.getSize())*3;
        }
        header.reserve(cells.size()*5);
        conn.reserve(conntotal);
        coords.reserve(coordtotal);
        for(const auto &cell:cells){
            const int connsize=static_cast<int>(cell.ElmtConn.size());
            const int nodesnum=cell.ElmtNodeCoords0.getSize();
            header.push_back(cell.Dim);
            header.push_back(cell.NodesNumPerElmt);
            header.push_back(cell.VTKCellType);
            header.push_back(connsize);
            header.push_back(nodesnum);
            for(int k=0;k<connsize;k++) conn.push_back(cell.ElmtConn[k]);
            const double *coordptr=cell.ElmtNodeCoords0.getDataPtr();
            for(int k=0;k<nodesnum*3;k++) coords.push_back(coordptr[k]);
        }
    }
    /**
     * Rebuild a mesh cell vector from the three flat buffers produced by packMeshCells.
     * @param cells  the mesh cell vector to be filled
     * @param ncells the number of mesh cells
     * @param header per-cell metadata (see packMeshCells)
     * @param conn   the flattened connectivity buffer
     * @param coords the flattened nodal coordinates buffer
     */
    void unpackMeshCells(vector<SingleMeshCell> &cells,const int &ncells,
                         const vector<int> &header,const vector<int> &conn,const vector<double> &coords){
        cells.clear();
        cells.resize(ncells);
        int hoff=0,coff=0,xoff=0;
        for(int i=0;i<ncells;i++){
            SingleMeshCell &cell=cells[i];
            cell.Dim            =header[hoff+0];
            cell.NodesNumPerElmt=header[hoff+1];
            cell.VTKCellType    =header[hoff+2];
            const int connsize  =header[hoff+3];
            const int nodesnum  =header[hoff+4];
            hoff+=5;
            cell.ElmtConn.resize(connsize,0);
            for(int k=0;k<connsize;k++) cell.ElmtConn[k]=conn[coff++];
            cell.ElmtNodeCoords0.resize(nodesnum);
            double *coordptr=cell.ElmtNodeCoords0.getData();
            for(int k=0;k<nodesnum*3;k++) coordptr[k]=coords[xoff++];
            cell.ElmtNodeCoords=cell.ElmtNodeCoords0;// deformed copy starts from the undeformed one
        }
    }
}

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

    MPI_Request request;

    /**
     * Pack the whole cell vector into flat buffers and send them in bulk. This
     * replaces the previous one-message-per-field-per-cell scheme (O(ncells)
     * messages) with a fixed number of messages, and keeps the coordinate send
     * buffer alive until MPI_Wait completes.
     */
    vector<int>    header;
    vector<int>    conn;
    vector<double> coords;
    packMeshCells(Meshcellvec,header,conn,coords);

    int ncells,headersize,connsize,coordssize;
    ncells    =static_cast<int>(Meshcellvec.size());
    headersize=static_cast<int>(header.size());
    connsize  =static_cast<int>(conn.size());
    coordssize=static_cast<int>(coords.size());

    // send out the length of the mesh cell vector
    MPI_Isend(&ncells,1,MPI_INT,Cpuid,Basetag+1,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);

    // send out the per-cell metadata
    MPI_Isend(&headersize,1,MPI_INT,Cpuid,Basetag+2,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    if(headersize>0){
        MPI_Isend(header.data(),headersize,MPI_INT,Cpuid,Basetag+3,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }

    // send out the flattened connectivity
    MPI_Isend(&connsize,1,MPI_INT,Cpuid,Basetag+4,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    if(connsize>0){
        MPI_Isend(conn.data(),connsize,MPI_INT,Cpuid,Basetag+5,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }

    // send out the flattened nodal coordinates
    MPI_Isend(&coordssize,1,MPI_INT,Cpuid,Basetag+6,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    if(coordssize>0){
        MPI_Isend(coords.data(),coordssize,MPI_DOUBLE,Cpuid,Basetag+7,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
}
//********************************************************
void MPIDataBus::receiveMeshCellFromMaster(vector<SingleMeshCell> &Meshcellvec,const int &Tag){
    int Basetag;
    Basetag=Tag;

    MPI_Request request;

    int ncells,headersize,connsize,coordssize;

    // receive the length of the mesh cell vector
    MPI_Irecv(&ncells,1,MPI_INT,0,Basetag+1,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);

    // receive the per-cell metadata
    MPI_Irecv(&headersize,1,MPI_INT,0,Basetag+2,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    vector<int> header(headersize,0);
    if(headersize>0){
        MPI_Irecv(header.data(),headersize,MPI_INT,0,Basetag+3,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }

    // receive the flattened connectivity
    MPI_Irecv(&connsize,1,MPI_INT,0,Basetag+4,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    vector<int> conn(connsize,0);
    if(connsize>0){
        MPI_Irecv(conn.data(),connsize,MPI_INT,0,Basetag+5,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }

    // receive the flattened nodal coordinates
    MPI_Irecv(&coordssize,1,MPI_INT,0,Basetag+6,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    vector<double> coords(coordssize,0.0);
    if(coordssize>0){
        MPI_Irecv(coords.data(),coordssize,MPI_DOUBLE,0,Basetag+7,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }

    unpackMeshCells(Meshcellvec,ncells,header,conn,coords);
}

//***********************************************************
void MPIDataBus::sendPhyName2MeshCellMapToOthers(const string &Phyname,const vector<SingleMeshCell> &Meshcellvec,const int &Tag,const int &Cpuid){
    int Basetag;
    Basetag=Tag;

    MPI_Request request;

    // send out the physical name to each rank
    int namesize;
    namesize=static_cast<int>(Phyname.size());
    MPI_Isend(&namesize,1,MPI_INT,Cpuid,Basetag+1,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    if(namesize>0){
        MPI_Isend(Phyname.c_str(),namesize,MPI_CHAR,Cpuid,Basetag+2,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }

    // pack and send the mesh cells in bulk (see sendMeshCellToOthers)
    vector<int>    header;
    vector<int>    conn;
    vector<double> coords;
    packMeshCells(Meshcellvec,header,conn,coords);

    int ncells,headersize,connsize,coordssize;
    ncells    =static_cast<int>(Meshcellvec.size());
    headersize=static_cast<int>(header.size());
    connsize  =static_cast<int>(conn.size());
    coordssize=static_cast<int>(coords.size());

    MPI_Isend(&ncells,1,MPI_INT,Cpuid,Basetag+3,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    MPI_Isend(&headersize,1,MPI_INT,Cpuid,Basetag+4,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    if(headersize>0){
        MPI_Isend(header.data(),headersize,MPI_INT,Cpuid,Basetag+5,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
    MPI_Isend(&connsize,1,MPI_INT,Cpuid,Basetag+6,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    if(connsize>0){
        MPI_Isend(conn.data(),connsize,MPI_INT,Cpuid,Basetag+7,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
    MPI_Isend(&coordssize,1,MPI_INT,Cpuid,Basetag+8,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    if(coordssize>0){
        MPI_Isend(coords.data(),coordssize,MPI_DOUBLE,Cpuid,Basetag+9,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
}
//*************************************************************************
void MPIDataBus::receivePhyName2MeshCellMapFromMaster(map<string,vector<SingleMeshCell>> &Localmap,const int &Tag){
    int Basetag;
    Basetag=Tag;

    MPI_Request request;

    // receive the physical name into a properly sized string (no fixed-size buffer)
    int namesize;
    MPI_Irecv(&namesize,1,MPI_INT,0,Basetag+1,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    string Phyname;
    Phyname.resize(namesize);
    if(namesize>0){
        MPI_Irecv(Phyname.data(),namesize,MPI_CHAR,0,Basetag+2,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }

    // receive the packed mesh cells (see receiveMeshCellFromMaster)
    int ncells,headersize,connsize,coordssize;
    MPI_Irecv(&ncells,1,MPI_INT,0,Basetag+3,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    MPI_Irecv(&headersize,1,MPI_INT,0,Basetag+4,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    vector<int> header(headersize,0);
    if(headersize>0){
        MPI_Irecv(header.data(),headersize,MPI_INT,0,Basetag+5,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
    MPI_Irecv(&connsize,1,MPI_INT,0,Basetag+6,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    vector<int> conn(connsize,0);
    if(connsize>0){
        MPI_Irecv(conn.data(),connsize,MPI_INT,0,Basetag+7,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
    MPI_Irecv(&coordssize,1,MPI_INT,0,Basetag+8,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    vector<double> coords(coordssize,0.0);
    if(coordssize>0){
        MPI_Irecv(coords.data(),coordssize,MPI_DOUBLE,0,Basetag+9,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }

    vector<SingleMeshCell> Meshcellvec;
    unpackMeshCells(Meshcellvec,ncells,header,conn,coords);
    Localmap[Phyname]=Meshcellvec;
}

//**************************************************************
void MPIDataBus::sendPhyID2MeshCellMapToOthers(const int &Phyid,const vector<SingleMeshCell> &Meshcellvec,const int &Tag,const int &Cpuid){
    int Basetag;
    Basetag=Tag;

    MPI_Request request;

    // send out the physical id to each rank
    MPI_Isend(&Phyid,1,MPI_INT,Cpuid,Basetag+1,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);

    // pack and send the mesh cells in bulk (see sendMeshCellToOthers)
    vector<int>    header;
    vector<int>    conn;
    vector<double> coords;
    packMeshCells(Meshcellvec,header,conn,coords);

    int ncells,headersize,connsize,coordssize;
    ncells    =static_cast<int>(Meshcellvec.size());
    headersize=static_cast<int>(header.size());
    connsize  =static_cast<int>(conn.size());
    coordssize=static_cast<int>(coords.size());

    MPI_Isend(&ncells,1,MPI_INT,Cpuid,Basetag+2,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    MPI_Isend(&headersize,1,MPI_INT,Cpuid,Basetag+3,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    if(headersize>0){
        MPI_Isend(header.data(),headersize,MPI_INT,Cpuid,Basetag+4,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
    MPI_Isend(&connsize,1,MPI_INT,Cpuid,Basetag+5,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    if(connsize>0){
        MPI_Isend(conn.data(),connsize,MPI_INT,Cpuid,Basetag+6,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
    MPI_Isend(&coordssize,1,MPI_INT,Cpuid,Basetag+7,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    if(coordssize>0){
        MPI_Isend(coords.data(),coordssize,MPI_DOUBLE,Cpuid,Basetag+8,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
}
//***************************************************************
void MPIDataBus::receivePhyID2MeshCellMapFromMaster(map<int,vector<SingleMeshCell>> &Localmap,const int &Tag){
    int Basetag;
    Basetag=Tag;

    int phyid;

    MPI_Request request;

    // receive the physical id from master rank
    MPI_Irecv(&phyid,1,MPI_INT,0,Basetag+1,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);

    // receive the packed mesh cells (see receiveMeshCellFromMaster)
    int ncells,headersize,connsize,coordssize;
    MPI_Irecv(&ncells,1,MPI_INT,0,Basetag+2,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    MPI_Irecv(&headersize,1,MPI_INT,0,Basetag+3,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    vector<int> header(headersize,0);
    if(headersize>0){
        MPI_Irecv(header.data(),headersize,MPI_INT,0,Basetag+4,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
    MPI_Irecv(&connsize,1,MPI_INT,0,Basetag+5,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    vector<int> conn(connsize,0);
    if(connsize>0){
        MPI_Irecv(conn.data(),connsize,MPI_INT,0,Basetag+6,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }
    MPI_Irecv(&coordssize,1,MPI_INT,0,Basetag+7,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    vector<double> coords(coordssize,0.0);
    if(coordssize>0){
        MPI_Irecv(coords.data(),coordssize,MPI_DOUBLE,0,Basetag+8,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
    }

    vector<SingleMeshCell> Meshcellvec;
    unpackMeshCells(Meshcellvec,ncells,header,conn,coords);
    Localmap[phyid]=Meshcellvec;
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

    int datasize;

    MPI_Request request;
    vector<int> nodeidvec;

    // receive the physical name length from master rank
    MPI_Irecv(&datasize,1,MPI_INT,0,Basetag+1,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    // receive the phy name into a properly sized string (no fixed-size buffer)
    phyname.resize(datasize);
    MPI_Irecv(phyname.data(),datasize,MPI_CHAR,0,Basetag+2,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);

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