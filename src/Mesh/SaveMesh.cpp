//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "Mesh/Mesh.h"

void Mesh::SaveMesh(){
    PetscMPIInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    if(rank==0){
        // only master rank can have the access to file write!
        ofstream meshout;
        meshout.open(_MeshFileName,ios::out);
        if(!meshout.is_open()){
            PetscPrintf(PETSC_COMM_WORLD,"*** Error: can\'t write mesh to vtu file(=%27s)!***\n",_MeshFileName.c_str());
            PetscPrintf(PETSC_COMM_WORLD,"***        please make sure you have write permission!                ***\n");
            Msg_AsFem_Exit();
        }
        int i,j,e;

        //****************************************
        //*** print out header information
        //****************************************
        meshout<<"<?xml version=\"1.0\"?>\n";
        meshout<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
        meshout<<"<UnstructuredGrid>\n";
        meshout<<"<Piece NumberOfPoints=\""<<GetNodesNum()<<"\" NumberOfCells=\""<<GetBulkElmtsNum()<<"\">\n";
    
    
        meshout<<"<Points>\n";
        meshout<<"<DataArray type=\"Float64\" Name=\"nodes\"  NumberOfComponents=\"3\"  format=\"ascii\">\n";
        //*****************************
        // print out node coordinates
        meshout<<scientific<<setprecision(6);
        for(i=1;i<=GetNodesNum();++i){
            meshout<<GetIthNodeJthCoord(i,1)<<" ";
            meshout<<GetIthNodeJthCoord(i,2)<<" ";
            meshout<<GetIthNodeJthCoord(i,3)<<"\n";
        }
        meshout<<"</DataArray>\n";
        meshout<<"</Points>\n";
        //***************************************
        //*** For cell information
        //***************************************
        meshout<<"<Cells>\n";
        meshout<<"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for(e=1;e<=GetBulkElmtsNum();++e){
            for(j=1;j<=GetIthBulkElmtNodesNum(e);++j){
                meshout<<GetIthBulkElmtJthConn(e,j)-1<<" ";
            }
            meshout<<"\n";
        }
        meshout<<"</DataArray>\n";
        // for offset
        meshout<<"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        int offset=0;
        for(e=1;e<=GetBulkElmtsNum();++e){
            offset+=GetIthBulkElmtNodesNum(e);
            meshout<<offset<<"\n";
        }
        meshout<<"</DataArray>\n";
        // For connectivity
        meshout<<"<DataArray type=\"Int32\" Name=\"types\"  NumberOfComponents=\"1\"  format=\"ascii\">\n";
        for(e=1;e<=GetBulkElmtsNum();++e){
            meshout<<GetIthBulkElmtVTKCellType(e)<<"\n";
        }
        meshout<<"</DataArray>\n";
        meshout<<"</Cells>\n";

    
        meshout<<"</Piece>\n";
        meshout<<"</UnstructuredGrid>\n";
        meshout<<"</VTKFile>"<<endl;
        meshout.close();
    }
}