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
//+++ Date   : 2024.02.03
//+++ Purpose: Implement the msh file (version-2) import function.
//+++          This mesh file must be the *.msh in version-2.
//+++          For version-4, please use Msh4FileImporter.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "FECell/MeshFile2FECellImporterBase.h"
#include "FECell/MshFileUtils.h"

/**
 * This class implements the msh file (version-2) import function.
 */
class Msh2File2FECellImporter:public MeshFile2FECellImporterBase{
protected:
    /**
     * This function implement the mesh file import.
     * @param filename the mesh file's string name
     * @param t_celldata the fe cell mesh data structure
     */
    virtual bool importMeshFile(const string &filename,FECellData &t_celldata) override;

private:
    /**
     * get the maximum dimension of the bulk mesh from msh file
     * @param filename the msh file name
     */
    int getMaxMeshDim(const string &filename)const;

};