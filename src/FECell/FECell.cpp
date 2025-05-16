//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author  : Yang Bai
//+++ Date    : 2024.01.19
//+++ Function: finite element cell management class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <fstream>
#include "FECell/FECell.h"
#include "FECell/FECellPartioner.h"

FECell::FECell() {
    m_CellData.MaxDim=1;
    m_CellData.MinDim=0;
    m_CellData.Nx=1;
    m_CellData.Xmin=0.0;
    m_CellData.Xmax=1.0;
    m_CellData.Ymin=0.0;
    m_CellData.Ymax=0.0;
    m_CellData.Zmin=0.0;
    m_CellData.Zmax=0.0;
    m_CellData.BulkElmtMeshType=MeshType::EDGE2;
    m_CellData.NodeCoords_Global.clear();
    //
    m_CellData.NodesNum=0;
    m_CellData.ElmtsNum=0;
    m_CellData.TotalDofsNum=0;
    // for physical group
    m_CellData.PhyGroupNum_Global=0;
    m_CellData.PhyDimVector_Global.clear();
    m_CellData.PhyNameVector_Global.clear();
    m_CellData.PhyIDVector_Global.clear();
    m_CellData.PhyGroupElmtsNumVector_Global.clear();

    m_CellData.NodalPhyGroupNum_Global=0;
}

// setup 1d mesh info
void FECell::setMeshInfo(const int &nx,const double &xmin,const double &xmax,const MeshType &meshtype){
    m_CellData.MaxDim=1;
    m_CellData.MinDim=0;
    m_CellData.Nx=nx;
    m_CellData.Xmin=xmin;
    m_CellData.Xmax=xmax;
    m_CellData.BulkElmtMeshType=meshtype;
}
// setup 2d mesh info
void FECell::setMeshInfo(const int &nx,const int &ny,
                         const double &xmin,const double &xmax,
                         const double &ymin,const double &ymax,const MeshType &meshtype){
    m_CellData.MaxDim=2;
    m_CellData.MinDim=1;
    m_CellData.Nx=nx;
    m_CellData.Ny=ny;
    m_CellData.Xmin=xmin;
    m_CellData.Xmax=xmax;
    m_CellData.Ymin=ymin;
    m_CellData.Ymax=ymax;
    m_CellData.BulkElmtMeshType=meshtype;
}
// setup 3d mesh info
void FECell::setMeshInfo(const int &nx,const int &ny,const int &nz,
                         const double &xmin,const double &xmax,
                         const double &ymin,const double &ymax,
                         const double &zmin,const double &zmax,const MeshType &meshtype){
    m_CellData.MaxDim=3;
    m_CellData.MinDim=2;
    m_CellData.Nx=nx;
    m_CellData.Ny=ny;
    m_CellData.Nz=nz;
    m_CellData.Xmin=xmin;
    m_CellData.Xmax=xmax;
    m_CellData.Ymin=ymin;
    m_CellData.Ymax=ymax;
    m_CellData.Zmin=zmin;
    m_CellData.Zmax=zmax;
    m_CellData.BulkElmtMeshType=meshtype;
}
//****************************************************
void FECell::saveFECell2VTUFile(const string filename)const{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if(rank==0){
        std::ofstream out;
        out.open(filename,std::ios::out);
        if(!out.is_open()){
            MessagePrinter::printErrorTxt("can\'t create/open"+filename+", please make sure you have the write permission");
            MessagePrinter::exitAsFem();
        }
        out << "<?xml version=\"1.0\"?>\n";
        out << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\">\n";
        out << "<UnstructuredGrid>\n";
        out << "<Piece NumberOfPoints=\"" << m_CellData.NodesNum << "\" NumberOfCells=\"" << m_CellData.BulkElmtsNum << "\">\n";
        out << "<Points>\n";
        out << "<DataArray type=\"Float64\" Name=\"nodes\"  NumberOfComponents=\"3\"  format=\"ascii\">\n";

        //*****************************
        // print out node coordinates
        out <<std::scientific << std::setprecision(6);
        for (int i = 0; i < m_CellData.NodesNum; i++){
            out << m_CellData.NodeCoords_Global[i*3+1-1] << " "
                << m_CellData.NodeCoords_Global[i*3+2-1] << " "
                << m_CellData.NodeCoords_Global[i*3+3-1] << "\n";
        }
        out << "</DataArray>\n";
        out << "</Points>\n";

        //***************************************
        //*** For cell information
        //***************************************
        out << "<Cells>\n";
        out << "<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for(const auto &cell:m_CellData.MeshCell_Total){
            for(const auto &id:cell.ElmtConn){
                out<<id-1<<" ";
            }
            out<<"\n";
        }
        out << "</DataArray>\n";

        //***************************************
        //*** For offset
        //***************************************
        out << "<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        int offset = 0;
        for(const auto &cell:m_CellData.MeshCell_Total){
            offset+=cell.NodesNumPerElmt;
            out << offset << "\n";
        }
        out << "</DataArray>\n";

        //***************************************
        //*** For vtk cell type
        //***************************************
        out << "<DataArray type=\"Int32\" Name=\"types\"  NumberOfComponents=\"1\"  format=\"ascii\">\n";
        for(const auto &cell:m_CellData.MeshCell_Total){
            out<<cell.VTKCellType<<"\n";
        }
        out << "</DataArray>\n";
        out << "</Cells>\n";

        //***************************************
        //*** End of output
        //***************************************
        out << "</Piece>\n";
        out << "</UnstructuredGrid>\n";
        out << "</VTKFile>" << endl;

        out.close();
    }// end-of-master-rank-if
}
//****************************************************
void FECell::saveFECellPartionInfo2VTUFile(const string filename)const {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if(rank==0){
        std::ofstream out;
        out.open(filename,std::ios::out);
        if(!out.is_open()){
            MessagePrinter::printErrorTxt("can\'t create/open"+filename+", please make sure you have the write permission");
            MessagePrinter::exitAsFem();
        }
        out << "<?xml version=\"1.0\"?>\n";
        out << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\">\n";
        out << "<UnstructuredGrid>\n";
        out << "<Piece NumberOfPoints=\"" << m_CellData.NodesNum << "\" NumberOfCells=\"" << m_CellData.BulkElmtsNum << "\">\n";

        // for cell data
        out << "<CellData>\n";
        // for ranks information
        out << "<DataArray type=\"Int32\" Name=\"cpu\" format=\"ascii\">\n";
        for (int e=1;e<=m_CellData.BulkElmtsNum;e++) {
            out << m_CellData.BulkCellPartionInfo_Global[e-1] << "\n";
        }
        out<<"</DataArray>\n";
        out << "</CellData>\n";


        // for points data
        out << "<Points>\n";
        out << "<DataArray type=\"Float64\" Name=\"nodes\"  NumberOfComponents=\"3\"  format=\"ascii\">\n";
        //*****************************
        // print out node coordinates
        out <<std::scientific << std::setprecision(6);
        for (int i = 0; i < m_CellData.NodesNum; i++){
            out << m_CellData.NodeCoords_Global[i*3+1-1] << " "
                << m_CellData.NodeCoords_Global[i*3+2-1] << " "
                << m_CellData.NodeCoords_Global[i*3+3-1] << "\n";
        }
        out << "</DataArray>\n";
        out << "</Points>\n";

        //***************************************
        //*** For cell information
        //***************************************
        out << "<Cells>\n";
        out << "<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for(const auto &cell:m_CellData.MeshCell_Total){
            for(const auto &id:cell.ElmtConn){
                out<<id-1<<" ";
            }
            out<<"\n";
        }
        out << "</DataArray>\n";

        //***************************************
        //*** For offset
        //***************************************
        out << "<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        int offset = 0;
        for(const auto &cell:m_CellData.MeshCell_Total){
            offset+=cell.NodesNumPerElmt;
            out << offset << "\n";
        }
        out << "</DataArray>\n";

        //***************************************
        //*** For vtk cell type
        //***************************************
        out << "<DataArray type=\"Int32\" Name=\"types\"  NumberOfComponents=\"1\"  format=\"ascii\">\n";
        for(const auto &cell:m_CellData.MeshCell_Total){
            out<<cell.VTKCellType<<"\n";
        }
        out << "</DataArray>\n";
        out << "</Cells>\n";

        //*** For cell data
        out << "<CellData>\n";
        //***************************************
        //*** For cell partition info
        //***************************************
        out << "<DataArray type=\"Int32\" Name=\"partition\"  NumberOfComponents=\"1\"  format=\"ascii\">\n";
        for(int e = 1; e <= m_CellData.BulkElmtsNum; e++){
            out << m_CellData.BulkCellPartionInfo_Global[e-1] << "\n";
        }
        out << "</DataArray>\n";
        out << "</CellData>\n";

        //***************************************
        //*** End of output
        //***************************************
        out << "</Piece>\n";
        out << "</UnstructuredGrid>\n";
        out << "</VTKFile>" << endl;

        out.close();
    }// end-of-master-rank-if
}
//****************************************************
void FECell::printSummaryInfo()const{
    MessagePrinter::printStars();
    MessagePrinter::printNormalTxt("Summary information of FECell system");
    MessagePrinter::printDashLine();
    string txt;
    char buff[65];
    snprintf(buff,65,"  total nodes=%5d, total elmts=%5d, total dofs=%6d",m_CellData.NodesNum,m_CellData.ElmtsNum,m_CellData.TotalDofsNum);
    txt=string(buff);
    MessagePrinter::printNormalTxt(txt);

    snprintf(buff,65,"  bulk elmts=%6d, line elmts=%6d, surf elmts=%6d",m_CellData.BulkElmtsNum,m_CellData.LineElmtsNum,m_CellData.SurfElmtsNum);
    txt=string(buff);
    MessagePrinter::printNormalTxt(txt);

    snprintf(buff,65,"  bulk elmt nodes=%2d, nodal dofs=%6d",m_CellData.NodesNumPerBulkElmt,m_CellData.MaxDofsPerNode);
    txt=string(buff);
    MessagePrinter::printNormalTxt(txt);

    MessagePrinter::printNormalTxt("  cell distribution method:"+MeshDistributionMethod);

    snprintf(buff,65,"  elemental physical groups=%4d, nodal physical groups=%4d",m_CellData.PhyGroupNum_Global,m_CellData.NodalPhyGroupNum_Global);
    txt=string(buff);
    MessagePrinter::printNormalTxt(txt);

    MessagePrinter::printDashLine();
    MessagePrinter::printNormalTxt("  physical ID<------>physical name<------>dim<---------->elmts num");
    for(int i=0;i<m_CellData.PhyGroupNum_Global;i++){
        snprintf(buff,65,"%8d    %20s         %2d          %8d",m_CellData.PhyIDVector_Global[i],m_CellData.PhyNameVector_Global[i].c_str(),m_CellData.PhyDimVector_Global[i],m_CellData.PhyGroupElmtsNumVector_Global[i]);
        txt=string(buff);
        MessagePrinter::printNormalTxt(txt);
    }

    MessagePrinter::printDashLine();
    MessagePrinter::printNormalTxt("  nodal physical ID<------>physical name<---------->nodes num");
    for(int i=0;i<m_CellData.NodalPhyGroupNum_Global;i++){
        snprintf(buff,65,"  %8d        %20s           %8d",m_CellData.NodalPhyIDVector_Global[i],m_CellData.NodalPhyNameVector_Global[i].c_str(),m_CellData.NodalPhyGroupNodesNumVector_Global[i]);
        txt=string(buff);
        MessagePrinter::printNormalTxt(txt);
    }

    int rank,size;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if(size==1){
        snprintf(buff,65,"  FECell is distributed in %5d cpu",size);
    }
    else{
        snprintf(buff,65,"  FECell is distributed in %5d cpus",size);
    }
    txt=string(buff);
    MessagePrinter::printDashLine();
    MessagePrinter::printNormalTxt(txt);

    if(rank!=0){
        MPI_Request request;
        int datasize;
        datasize=static_cast<int>(m_CellData.MeshCell_Local.size());
        MPI_Isend(&datasize,1,MPI_INT,0,1,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        datasize=static_cast<int>(m_CellData.PhyName2MeshCellVectorMap_Local.size());
        MPI_Isend(&datasize,1,MPI_INT,0,2,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        for(const auto &it:m_CellData.PhyName2MeshCellVectorMap_Local){
            datasize=static_cast<int>(it.first.size());
            MPI_Isend(&datasize,1,MPI_INT,0,3,MPI_COMM_WORLD,&request);
            MPI_Wait(&request,MPI_STATUS_IGNORE);

            MPI_Isend(it.first.c_str(),datasize,MPI_CHAR,0,4,MPI_COMM_WORLD,&request);
            MPI_Wait(&request,MPI_STATUS_IGNORE);

            datasize=static_cast<int>(it.second.size());
            MPI_Isend(&datasize,1,MPI_INT,0,5,MPI_COMM_WORLD,&request);
            MPI_Wait(&request,MPI_STATUS_IGNORE);
        }
        // for nodal physical info
        datasize=0;
        for(const auto &it:m_CellData.NodalPhyName2NodeIDVecMap_Local) datasize+=static_cast<int>(it.second.size());
        MPI_Isend(&datasize,1,MPI_INT,0,6,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        datasize=static_cast<int>(m_CellData.NodalPhyName2NodeIDVecMap_Local.size());
        MPI_Isend(&datasize,1,MPI_INT,0,7,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        for(const auto &it:m_CellData.NodalPhyName2NodeIDVecMap_Local){
            datasize=static_cast<int>(it.first.size());
            MPI_Isend(&datasize,1,MPI_INT,0,8,MPI_COMM_WORLD,&request);
            MPI_Wait(&request,MPI_STATUS_IGNORE);

            MPI_Isend(it.first.c_str(),datasize,MPI_CHAR,0,9,MPI_COMM_WORLD,&request);
            MPI_Wait(&request,MPI_STATUS_IGNORE);

            datasize=static_cast<int>(it.second.size());
            MPI_Isend(&datasize,1,MPI_INT,0,10,MPI_COMM_WORLD,&request);
            MPI_Wait(&request,MPI_STATUS_IGNORE);
        }
        datasize=static_cast<int>(m_CellData.NodeIDs_Local.size());
        MPI_Isend(&datasize,1,MPI_INT,0,11,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        //
        if (datasize>0) {
            MPI_Isend(&m_CellData.NodeIDs_Local[0],1,MPI_INT,0,12,MPI_COMM_WORLD,&request);
            MPI_Wait(&request,MPI_STATUS_IGNORE);
            //
            MPI_Isend(&m_CellData.NodeIDs_Local[datasize-1],1,MPI_INT,0,13,MPI_COMM_WORLD,&request);
            MPI_Wait(&request,MPI_STATUS_IGNORE);
        }
    } // end-of-other-ranks
    else{
        int cpuid,datasize;
        cpuid=0;
        MessagePrinter::printDashLine();
        snprintf(buff,65,"  For rank-%3d:",cpuid);
        txt=string(buff);
        MessagePrinter::printNormalTxt(txt);

        snprintf(buff,65,"    Local bulk elmts num=%6d",static_cast<int>(m_CellData.MeshCell_Local.size()));
        txt=string(buff);
        MessagePrinter::printNormalTxt(txt);
        for(const auto &it:m_CellData.PhyName2MeshCellVectorMap_Local){
            snprintf(buff,65,"      elmtset[%20s]====> elmts num=%6d",it.first.c_str(),static_cast<int>(it.second.size()));
            txt=string(buff);
            MessagePrinter::printNormalTxt(txt);
        }

        datasize=0;
        for(const auto &it:m_CellData.NodalPhyName2NodeIDVecMap_Local) datasize+=static_cast<int>(it.second.size());

        snprintf(buff,65,"    Local nodes num=%6d",datasize);
        txt=string(buff);
        MessagePrinter::printNormalTxt(txt);
        for(const auto &it:m_CellData.NodalPhyName2NodeIDVecMap_Local){
            snprintf(buff,65,"      nodeset[%20s]====> nodes num=%6d",it.first.c_str(),static_cast<int>(it.second.size()));
            txt=string(buff);
            MessagePrinter::printNormalTxt(txt);
        }

        datasize=static_cast<int>(m_CellData.NodeIDs_Local.size());
        if (datasize > 0) {
            snprintf(buff,65,"    Local nodes ids=%6d(%6d~%6d)",datasize,m_CellData.NodeIDs_Local[0],m_CellData.NodeIDs_Local[datasize-1]);
        }
        else {
            snprintf(buff,65,"    Local nodes ids=%6d(%6d~%6d)",datasize,0,0);
        }
        txt=string(buff);
        MessagePrinter::printNormalTxt(txt);

        MPI_Request req;
        int phynum,nodesnum;
        string phyname;

        for(cpuid=1;cpuid<size;cpuid++){
            MPI_Irecv(&datasize,1,MPI_INT,cpuid,1,MPI_COMM_WORLD,&req);
            MPI_Wait(&req,MPI_STATUS_IGNORE);

            MessagePrinter::printDashLine();
            snprintf(buff,65,"  For rank-%3d:",cpuid);
            txt=string(buff);
            MessagePrinter::printNormalTxt(txt);
            
            snprintf(buff,65,"    Local bulk elmts num=%6d",datasize);
            txt=string(buff);
            MessagePrinter::printNormalTxt(txt);

            MPI_Irecv(&phynum,1,MPI_INT,cpuid,2,MPI_COMM_WORLD,&req);
            MPI_Wait(&req,MPI_STATUS_IGNORE);

            for(int i=0;i<phynum;i++){
                MPI_Irecv(&datasize,1,MPI_INT,cpuid,3,MPI_COMM_WORLD,&req);
                MPI_Wait(&req,MPI_STATUS_IGNORE);

                phyname.clear();
                MPI_Irecv(buff,datasize,MPI_CHAR,cpuid,4,MPI_COMM_WORLD,&req);
                MPI_Wait(&req,MPI_STATUS_IGNORE);
                for(int j=0;j<datasize;j++) phyname.push_back(buff[j]);

                MPI_Irecv(&datasize,1,MPI_INT,cpuid,5,MPI_COMM_WORLD,&req);
                MPI_Wait(&req,MPI_STATUS_IGNORE);

                snprintf(buff,65,"      elmtset[%20s]====> elmts num=%6d",phyname.c_str(),datasize);
                txt=string(buff);
                MessagePrinter::printNormalTxt(txt);
            }

            // for nodal phy info
            nodesnum=0;
            MPI_Irecv(&nodesnum,1,MPI_INT,cpuid,6,MPI_COMM_WORLD,&req);
            MPI_Wait(&req,MPI_STATUS_IGNORE);

            snprintf(buff,65,"    Local nodes num=%6d",nodesnum);
            txt=string(buff);
            MessagePrinter::printNormalTxt(txt);

            MPI_Irecv(&phynum,1,MPI_INT,cpuid,7,MPI_COMM_WORLD,&req);
            MPI_Wait(&req,MPI_STATUS_IGNORE);

            for(int i=0;i<phynum;i++){
                MPI_Irecv(&datasize,1,MPI_INT,cpuid,8,MPI_COMM_WORLD,&req);
                MPI_Wait(&req,MPI_STATUS_IGNORE);

                phyname.clear();
                MPI_Irecv(buff,datasize,MPI_CHAR,cpuid,9,MPI_COMM_WORLD,&req);
                MPI_Wait(&req,MPI_STATUS_IGNORE);
                for(int j=0;j<datasize;j++) phyname.push_back(buff[j]);

                MPI_Irecv(&datasize,1,MPI_INT,cpuid,10,MPI_COMM_WORLD,&req);
                MPI_Wait(&req,MPI_STATUS_IGNORE);

                snprintf(buff,65,"      nodeset[%20s]====> nodes num=%6d",phyname.c_str(),datasize);
                txt=string(buff);
                MessagePrinter::printNormalTxt(txt);
            }

            MPI_Irecv(&datasize,1,MPI_INT,cpuid,11,MPI_COMM_WORLD,&req);
            MPI_Wait(&req,MPI_STATUS_IGNORE);
            int min,max;
            min=-1;max=-1;
            if (datasize>0) {
                MPI_Irecv(&min,1,MPI_INT,cpuid,12,MPI_COMM_WORLD,&req);
                MPI_Wait(&req,MPI_STATUS_IGNORE);
                MPI_Irecv(&max,1,MPI_INT,cpuid,13,MPI_COMM_WORLD,&req);
                MPI_Wait(&req,MPI_STATUS_IGNORE);
            }
            else {
                min=0;max=0;
            }
            snprintf(buff,65,"    Local nodes ids=%6d(%6d~%6d)",datasize,min,max);
            txt=string(buff);
            MessagePrinter::printNormalTxt(txt);
        }
    }
    MessagePrinter::printStars();
}

void FECell::printMeshInfo()const{
    MessagePrinter::printStars();
    MessagePrinter::printNormalTxt("Summary information of FECell system");
    MessagePrinter::printDashLine();
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if(rank==0){
        char buff1[6+17],buff2[7];
        string str;
        for(int e=1;e<=m_CellData.BulkElmtsNum;e++){
            str.clear();
            snprintf(buff1,6+17,"%5d-th element:",e);
            str+=buff1;
            for(int i=1;i<=m_CellData.NodesNumPerBulkElmt;i++){
                snprintf(buff2,7,"%5d ",m_CellData.MeshCell_Total[e-1].ElmtConn[i-1]);
                str+=buff2;
            }
            MessagePrinter::printNormalTxt(str);
        }
    }
    MessagePrinter::printDashLine();
    MessagePrinter::printStars();
}

void FECell::releaseMemory(){

}

void FECell::distributeMesh(){
    FECellPartioner Part;
    Part.partFECell(MeshDistributionMethod,m_CellData);
    MPI_Barrier(MPI_COMM_WORLD);
}