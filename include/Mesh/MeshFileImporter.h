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

#pragma once

#include "Mesh/Msh2FileImporter.h"
#include "Mesh/Msh4FileImporter.h"
#include "Mesh/Gmsh2FileImporter.h"

/**
 * This class offers the mesh file import function for the supported mesh file format
 */
class MeshFileImporter:public Msh2FileImporter,
                       public Msh4FileImporter,
                       public Gmsh2FileImporter{
public:
    /**
     * constructor
     */
    MeshFileImporter();

    /**
     * import the mesh from mesh file (msh2 format)
     * @param meshfile the string name for the mesh file
     * @param t_meshdata the mesh data structure
     */
    bool importMsh2Mesh(const string meshfile,MeshData &t_meshdata);

    /**
     * import the mesh from mesh file (msh4 format)
     * @param meshfile the string name for the mesh file
     * @param t_meshdata the mesh data structure
     */
    bool importMsh4Mesh(const string meshfile,MeshData &t_meshdata);

    /**
     * import the mesh from gmsh2 file (msh2 format from netgen)
     * @param meshfile the string name for the mesh file
     * @param t_meshdata the mesh data structure
     */
    bool importGmsh2Mesh(const string meshfile,MeshData &t_meshdata);

};