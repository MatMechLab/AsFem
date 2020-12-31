//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.06.30
//+++ Purpose: Here we implement the save function for our lagrange
//+++          mesh. Currently, the mesh can only be saved as vtu
//+++          file!
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Mesh/LagrangeMesh.h"

void LagrangeMesh::SaveLagrangeMesh(string inputfilename) const{
    PetscMPIInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    if(rank==0){
        char buff[110];
        string str;
        string _MeshFileName;
        if(inputfilename.size()<2){
            _MeshFileName="mesh.vtu";
        }
        else{
            int i=inputfilename.find_last_of(".");
            _MeshFileName=inputfilename.substr(0,i)+"_mesh"+".vtu";
        }
        ofstream meshout;
        meshout.open(_MeshFileName,ios::out);
        if(!meshout.is_open()){
            snprintf(buff,110,"can\'t write mesh to vtu file(=%28s), please make sure you have the write permission",_MeshFileName.c_str());
            str=buff;
            MessagePrinter::PrintErrorTxt(str);
            MessagePrinter::AsFem_Exit();
        }
        int i,j,e;

        //****************************************
        //*** print out header information
        //****************************************
        meshout<<"<?xml version=\"1.0\"?>\n";
        meshout<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
        meshout<<"<UnstructuredGrid>\n";
        meshout<<"<Piece NumberOfPoints=\""<<GetBulkMeshNodesNum()
               <<"\" NumberOfCells=\""<<GetBulkMeshBulkElmtsNum()<<"\">\n";
    
    
        meshout<<"<Points>\n";
        meshout<<"<DataArray type=\"Float64\" Name=\"nodes\"  NumberOfComponents=\"3\"  format=\"ascii\">\n";
        //*****************************
        // print out node coordinates
        meshout<<scientific<<setprecision(6);
        for(i=1;i<=GetBulkMeshNodesNum();++i){
            meshout<<GetBulkMeshIthNodeJthCoord(i,1)<<" ";
            meshout<<GetBulkMeshIthNodeJthCoord(i,2)<<" ";
            meshout<<GetBulkMeshIthNodeJthCoord(i,3)<<"\n";
        }
        meshout<<"</DataArray>\n";
        meshout<<"</Points>\n";
        //***************************************
        //*** For cell information
        //***************************************
        meshout<<"<Cells>\n";
        meshout<<"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for(e=1;e<=GetBulkMeshBulkElmtsNum();++e){
            for(j=1;j<=GetBulkMeshIthBulkElmtNodesNum(e);++j){
                meshout<<GetBulkMeshIthBulkElmtJthNodeID(e,j)-1<<" ";
            }
            meshout<<"\n";
        }
        meshout<<"</DataArray>\n";
        // for offset
        meshout<<"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        int offset=0;
        for(e=1;e<=GetBulkMeshBulkElmtsNum();++e){
            offset+=GetBulkMeshIthBulkElmtNodesNum(e);
            meshout<<offset<<"\n";
        }
        meshout<<"</DataArray>\n";
        // For connectivity
        meshout<<"<DataArray type=\"Int32\" Name=\"types\"  NumberOfComponents=\"1\"  format=\"ascii\">\n";
        for(e=1;e<=GetBulkMeshBulkElmtsNum();++e){
            meshout<<GetBulkMeshIthBulkElmtVTKCellType(e)<<"\n";
        }
        meshout<<"</DataArray>\n";
        meshout<<"</Cells>\n";

    
        meshout<<"</Piece>\n";
        meshout<<"</UnstructuredGrid>\n";
        meshout<<"</VTKFile>"<<endl;
        meshout.close();
    }
}