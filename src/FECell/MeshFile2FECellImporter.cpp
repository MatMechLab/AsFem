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
//+++ Date   : 2024.08.01
//+++ Purpose: import mesh file and generate the fecell data structure
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FECell/MeshFile2FECellImporter.h"


MeshFile2FECellImporter::MeshFile2FECellImporter(){}

bool MeshFile2FECellImporter::importMeshFile2FECell(const string &filetypename,const string &filename,FECellData &t_celldata){
    if(filetypename.find("msh4")!=string::npos){
        return Msh4File2FECellImporter::importMeshFile(filename,t_celldata);
    }
    else if(filetypename.find("msh2")!=string::npos && filetypename.find("gmsh2")==string::npos){
        return Msh2File2FECellImporter::importMeshFile(filename,t_celldata);
    }
    else if(filetypename.find("gmsh2")!=string::npos){
        return Gmsh2File2FECellImporter::importMeshFile(filename,t_celldata);
    }
    else{
        MessagePrinter::printErrorTxt("Unsupported mesh file import option, currently, only msh2, msh4, and gmsh2 file are supported");
        MessagePrinter::exitAsFem();
    }
    return false;
}