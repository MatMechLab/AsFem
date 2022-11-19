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
//+++ Date   : 2022.08.26
//+++ Purpose: Implement the general mesh file import function
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Mesh/MeshFileImporter.h"

MeshFileImporter::MeshFileImporter(){

}

bool MeshFileImporter::importMsh2Mesh(const string meshfile,MeshData &t_meshdata){
    return Msh2FileImporter::importMeshFile(meshfile,t_meshdata);
}


bool MeshFileImporter::importMsh4Mesh(const string meshfile,MeshData &t_meshdata){
    return Msh4FileImporter::importMeshFile(meshfile,t_meshdata);
}


bool MeshFileImporter::importGmsh2Mesh(const string meshfile,MeshData &t_meshdata){
    return Gmsh2FileImporter::importMeshFile(meshfile,t_meshdata);
}