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
//+++ Date   : 2024.01.19
//+++ Purpose: the FE cell generator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "FECell/Lagrange3DMeshCellGenerator.h"


/**
 * This class offers the FE mesh cell generator for default mesh in AsFem
*/
class FECellGenerator:public Lagrange3DMeshCellGenerator{
public:
    /**
     * Constructor
    */
    FECellGenerator();

    /**
     * create the finite element mesh cell and store them in meshdata(distributed)
     * @param meshtype the mesh type used for mesh generation
     * @param meshdata the meshdata structure
    */
    bool createFEMeshCell(const MeshType &meshtype,FECellData &meshdata);
};