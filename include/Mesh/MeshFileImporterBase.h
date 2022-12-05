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
//+++ Purpose: Define the abstract class for mesh import function.
//+++          Other meshio, i.e. gmsh and netgen should inherit
//+++          from this
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

#include "Mesh/MeshData.h"
#include "Mesh/MeshType.h"

using std::transform;
using std::remove;
using std::ifstream;
using std::ios;
using std::istringstream;
using std::make_pair;

/**
 * This is the abstract class for mesh file reader, if one wants to have the specific mesh file import, one 
 * should inherit this class and offers the implementation details
 */
class MeshFileImporterBase{
protected:
    /**
     * This function implement the mesh file import.
     * @param filename the mesh file's string name
     * @param t_meshdata the mesh data structure
     */
    virtual bool importMeshFile(const string &filename,MeshData &t_meshdata)=0;
    
};
