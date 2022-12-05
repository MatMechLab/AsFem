//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2022.05.07
//+++ Purpose: save mesh by using different file format
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Mesh/BulkMesh.h"

void BulkMesh::saveBulkMesh2VTU(const string &inputfilename)const{
    // here the inputfilename must be "xxxx.json", which should has the .json extension
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
            _MeshFileName=inputfilename.substr(0,i)+"-mesh"+".vtu";
        }
        std::ofstream meshout;
        meshout.open(_MeshFileName,std::ios::out);
        if(!meshout.is_open()){
            snprintf(buff,110,"can\'t write mesh to vtu file(=%28s), please make sure you have the write permission",_MeshFileName.c_str());
            str=buff;
            MessagePrinter::printErrorTxt(str);
            MessagePrinter::exitAsFem();
        }
        int i,j,e;

        //****************************************
        //*** print out header information
        //****************************************
        meshout<<"<?xml version=\"1.0\"?>\n";
        meshout<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
        meshout<<"<UnstructuredGrid>\n";
        meshout<<"<Piece NumberOfPoints=\""<<getBulkMeshNodesNum()
               <<"\" NumberOfCells=\""<<getBulkMeshBulkElmtsNum()<<"\">\n";
    
    
        meshout<<"<Points>\n";
        meshout<<"<DataArray type=\"Float64\" Name=\"nodes\"  NumberOfComponents=\"3\"  format=\"ascii\">\n";
        //*****************************
        // print out node coordinates
        meshout<<std::scientific<<std::setprecision(6);
        for(i=1;i<=getBulkMeshNodesNum();++i){
            meshout<<getBulkMeshIthNodeJthCoord0(i,1)<<" ";
            meshout<<getBulkMeshIthNodeJthCoord0(i,2)<<" ";
            meshout<<getBulkMeshIthNodeJthCoord0(i,3)<<"\n";
        }
        meshout<<"</DataArray>\n";
        meshout<<"</Points>\n";
        //***************************************
        //*** For cell information
        //***************************************
        meshout<<"<Cells>\n";
        meshout<<"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for(e=1;e<=getBulkMeshBulkElmtsNum();++e){
            for(j=1;j<=getBulkMeshNodesNumPerBulkElmt();++j){
                meshout<<getBulkMeshIthBulkElmtJthNodeID(e,j)-1<<" ";
            }
            meshout<<"\n";
        }
        meshout<<"</DataArray>\n";
        // for offset
        meshout<<"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        int offset=0;
        for(e=1;e<=getBulkMeshBulkElmtsNum();++e){
            offset+=getBulkMeshNodesNumPerBulkElmt();
            meshout<<offset<<"\n";
        }
        meshout<<"</DataArray>\n";
        // For connectivity
        meshout<<"<DataArray type=\"Int32\" Name=\"types\"  NumberOfComponents=\"1\"  format=\"ascii\">\n";
        for(e=1;e<=getBulkMeshBulkElmtsNum();++e){
            meshout<<getBulkMeshBulkElmtVTKCellType()<<"\n";
        }
        meshout<<"</DataArray>\n";
        meshout<<"</Cells>\n";

    
        meshout<<"</Piece>\n";
        meshout<<"</UnstructuredGrid>\n";
        meshout<<"</VTKFile>"<<endl;
        meshout.close();
    }
}