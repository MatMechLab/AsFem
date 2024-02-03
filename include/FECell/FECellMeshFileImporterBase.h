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
//+++ Date   : 2022.08.26
//+++ Purpose: Define the abstract class for FE cell mesh import function.
//+++          Other meshio, i.e. gmsh and netgen should inherit
//+++          from this one
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdio>
#include <vector>
#include <algorithm>


#include "Utils/MessagePrinter.h"

#include "FECell/FECellData.h"
#include "Mesh/MeshType.h"

using std::transform;
using std::remove;
using std::ifstream;
using std::ios;
using std::istringstream;
using std::make_pair;

/**
 * This is the abstract class for fe cell mesh file reader, if one wants to have the specific mesh file importer, one 
 * should inherit this class and offers the implementation details
 */
class FECellMeshFileImporterBase{
protected:
    /**
     * This function implement the mesh file import.
     * @param filename the mesh file's string name
     * @param t_celldata the fe cell mesh data structure
     */
    virtual bool importMeshFile(const string &filename,FECellData &t_celldata)=0;
    
};