//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.06.26
//+++ Purpose: implement a general mesh io class, which can import
//+++          mesh file from gmsh, netgen, etc...
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#include "Mesh/MeshIO.h"

MeshIO::MeshIO(){
    Gmsh2IO::_MeshFileName.clear();
    Gmsh4IO::_MeshFileName.clear();
    
    Gmsh2IO::_HasSetMeshFileName=false;
    Gmsh4IO::_HasSetMeshFileName=false;

    _MeshIOType=MeshIOType::GMSH2;
}

bool MeshIO::ReadMeshFromFile(Mesh &mesh){
    switch (_MeshIOType)
    {
    case MeshIOType::GMSH2:
        // cout<<"using gmsh2"<<endl;
        return Gmsh2IO::ReadMeshFromFile(mesh);
        break;
    case MeshIOType::GMSH4:
        // cout<<"using gmsh4"<<endl;
        return Gmsh4IO::ReadMeshFromFile(mesh);
        break;
    case MeshIOType::NETGEN:
        MessagePrinter::PrintErrorTxt("Netgen mesh is not supported yet!");
        return false;
    case MeshIOType::ABAQUS:
        return AbaqusIO::ReadMeshFromFile(mesh);
    default:
        return false;
    }
}
void MeshIO::SetMeshFileName(string filename){
    if(IsGmsh2MeshFile(filename)){
        Gmsh2IO::SetMeshFileName(filename);
        _MeshIOType=MeshIOType::GMSH2;
    }
    else if(IsGmsh4MeshFile(filename)){
        Gmsh4IO::SetMeshFileName(filename);
        _MeshIOType=MeshIOType::GMSH4;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
    }
    else if(IsNetgenMeshFile(filename)){
        _MeshIOType=MeshIOType::NETGEN;
    }
    else if(IsAbaqusMeshFile(filename)){
        AbaqusIO::SetMeshFileName(filename);
        _MeshIOType=MeshIOType::ABAQUS;
    }
    else{
        MessagePrinter::PrintErrorTxt("can\'t read mesh file(="+filename+")!!!, the mesh file type is not supported yet");
        MessagePrinter::AsFem_Exit();
    }
}
string MeshIO::GetMeshFileName()const{
    switch (_MeshIOType)
    {
    case MeshIOType::GMSH2:
        return Gmsh2IO::GetMeshFileName();
        break;
    case MeshIOType::GMSH4:
        return Gmsh4IO::GetMeshFileName();
        break;
    case MeshIOType::ABAQUS:
        return AbaqusIO::GetMeshFileName();
        break;
    default:
        return "";
    }
}