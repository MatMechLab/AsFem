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
//+++ Purpose: the 1D edge lagrange FE cell generator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "FECell/FECellGeneratorBase.h"


/**
 * the cell generator for 1d lagrange mesh
 */
class Lagrange1DEdge3MeshCellGenerator:public FECellGeneratorBase{
public:
    /**
     * constructor
     */
    Lagrange1DEdge3MeshCellGenerator();
    /**
     * deconstructor
     */
    ~Lagrange1DEdge3MeshCellGenerator();

    /**
     * function for the details of 3d lagrange mesh generation, if everything works fine, it will return true.
     * @param t_meshtype the type of mesh one want to use for the mesh generation
     * @param t_celldata the fe cell data structure, which should be updated within each cell generator!
     */
    virtual bool generateFECell(FECellData &t_celldata) override;
private:
    bool m_mesh_generated=false;

};