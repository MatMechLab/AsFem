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
//+++ Purpose: Implement the msh file (version-2) import function.
//+++          This mesh file must be the *.msh in version-2.
//+++          For version-4, please use Msh4FileImporter.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "Mesh/MeshFileImporterBase.h"
#include "Mesh/MshFileUtils.h"

/**
 * This class implements the msh file (version-2) import function.
 */
class Msh2FileImporter:public MeshFileImporterBase{
protected:
    /**
     * This function implement the mesh file import.
     * @param filename the mesh file's string name
     * @param t_meshdata the mesh data structure
     */
    virtual bool importMeshFile(const string &filename,MeshData &t_meshdata) override;

private:
    /**
     * get the maximum dimension of the bulk mesh from msh file
     * @param filename the msh file name
     */
    int getMaxMeshDim(const string &filename)const;

};