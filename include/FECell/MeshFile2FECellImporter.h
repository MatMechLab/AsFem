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

#pragma once

#include "FECell/Msh2File2FECellImporter.h"
#include "FECell/Msh4File2FECellImporter.h"
#include "FECell/Gmsh2File2FECellImporter.h"


class MeshFile2FECellImporter:public Msh2File2FECellImporter,
                              public Msh4File2FECellImporter,
                              public Gmsh2File2FECellImporter{
public:
    /**
     * Constructor
     */
    MeshFile2FECellImporter();

    /**
     * read the mesh file and generate the fecell data structure
     * @param filetypename the type name of the mesh file
     * @param filename mesh file name
     * @param t_celldata the fe cell data structure
     */
    bool importMeshFile2FECell(const string &filetypename,const string &filename,FECellData &t_celldata);

};