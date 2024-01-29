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

FECell::FECell(){}

// setup 1d mesh info
void FECell::setMeshInfo(const int &nx,const double &xmin,const double &xmax,const MeshType &meshtype){
    m_CellData.Nx=nx;
    m_CellData.Xmin=xmin;
    m_CellData.Xmax=xmax;
    m_CellData.BulkElmtMeshType=meshtype;
}
// setup 2d mesh info
void FECell::setMeshInfo(const int &nx,const int &ny,
                         const double &xmin,const double &xmax,
                         const double &ymin,const double &ymax,const MeshType &meshtype){
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
void FECell::saveFECell2VTUFile()const{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if(rank==0){
        std::ofstream out;
        out.open("mesh.vtu",std::ios::out);
        if(!out.is_open()){
            MessagePrinter::printErrorTxt("can\'t create/open mesh.vtu, please make sure you have the write permission");
            MessagePrinter::exitAsFem();
        }
        out << "<?xml version=\"1.0\"?>\n";
        out << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\">\n";
        out << "<UnstructuredGrid>\n";
        out << "<Piece NumberOfPoints=\"" << m_CellData.NodesNum << "\" NumberOfCells=\"" << m_CellData.BulkElmtsNum << "\">\n";
        out << "<Points>\n";
        out << "<DataArray type=\"Float64\" Name=\"nodes\"  NumberOfComponents=\"3\"  format=\"ascii\">\n";

        int i;
        //*****************************
        // print out node coordinates
        out <<std::scientific << std::setprecision(6);
        for (i = 1; i <= m_CellData.NodesNum; i++){
            out << m_CellData.NodeCoords_Global[(i-1)*3+1-1] << " ";
            out << m_CellData.NodeCoords_Global[(i-1)*3+2-1] << " ";
            out << m_CellData.NodeCoords_Global[(i-1)*3+3-1] << "\n";
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
    }
}
void FECell::printSummaryInfo()const{
    MessagePrinter::printStars();
    MessagePrinter::printNormalTxt("Summary information of FECell system");
    MessagePrinter::printDashLine();
    string txt;
    char buff[65];
    snprintf(buff,65,"Total nodes=%5d, total elmts=%5d, total dofs=%6d",m_CellData.NodesNum,m_CellData.ElmtsNum,m_CellData.TotalDofsNum);
    txt=string(buff);
    MessagePrinter::printNormalTxt(txt);

    snprintf(buff,65,"Bulk elmts=%6d, line elmts=%6d, surf elmts=%6d",m_CellData.BulkElmtsNum,m_CellData.LineElmtsNum,m_CellData.SurfElmtsNum);
    txt=string(buff);
    MessagePrinter::printNormalTxt(txt);

    snprintf(buff,65,"Bulk elmt nodes=%2d, nodal dofs=%6d",m_CellData.NodesNumPerBulkElmt,m_CellData.MaxDofsPerNode);
    txt=string(buff);
    MessagePrinter::printNormalTxt(txt);

    snprintf(buff,65,"Elemental physical groups=%4d, nodal physical groups=%4d",m_CellData.PhyGroupNum_Global,m_CellData.NodalPhyGroupNum_Global);
    txt=string(buff);
    MessagePrinter::printNormalTxt(txt);

    MessagePrinter::printDashLine();
    MessagePrinter::printNormalTxt("Physical ID<------>Physical name<------>Dim<---------->Elmts num");
    for(int i=0;i<m_CellData.PhyGroupNum_Global;i++){
        snprintf(buff,65,"%8d    %20s         %2d          %8d",m_CellData.PhyIDVector_Global[i],m_CellData.PhyNameVector_Global[i].c_str(),m_CellData.PhyDimVector_Global[i],m_CellData.PhyGroupElmtsNumVector_Global[i]);
        txt=string(buff);
        MessagePrinter::printNormalTxt(txt);
    }

    MessagePrinter::printDashLine();
    MessagePrinter::printNormalTxt("Nodal physical ID<------>Physical name<---------->Nodes num");
    for(int i=0;i<m_CellData.NodalPhyGroupNum_Global;i++){
        snprintf(buff,65,"  %8d        %20s           %8d",m_CellData.NodalPhyIDVector_Global[i],m_CellData.NodalPhyNameVector_Global[i].c_str(),m_CellData.NodalPhyGroupNodesNumVector_Global[i]);
        txt=string(buff);
        MessagePrinter::printNormalTxt(txt);
    }

    int rank,size;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if(size==1){
        snprintf(buff,65,"FECell is distributed in %5d cpu",size);
    }
    else{
        snprintf(buff,65,"FECell is distributed in %5d cpus",size);
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
    }
    else{
        int cpuid,datasize;
        char buff[65];
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

        }
    }
    MessagePrinter::printStars();
}