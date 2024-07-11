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
//+++ Date   : 2024.07.11
//+++ Purpose: Implement the msh file (version-4) import function.
//+++          This mesh file must be the *.msh in version-4.
//+++          For version-2, please use Msh2FileImporter.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "FECell/MeshFile2FECellImporterBase.h"
#include "FECell/MshFileUtils.h"

/**
 * This class implements the msh file (version-2) import function.
 */
class Msh4File2FECellImporter:public MeshFile2FECellImporterBase{
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

    /**
     * get the maximum dimension of the bulk mesh from msh file
     * @param entityDim the dimension of current entity
     * @param entityTag the tag of current entity
     */
    int getPhysicalIDViaEntityTag(const int &entityDim,const int &entityTag)const;

private:
    vector<int> m_PointsEntityPhyIDs,m_CurvesEntityPhyIDS,m_SurfaceEntityPhyIDs,m_VolumesEntityPhyIDs;

};