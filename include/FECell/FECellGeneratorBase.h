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
//+++ Date   : 2024.01.06
//+++ Purpose: the base class for FE cell generator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <algorithm>

#include "mpi.h"
#include "FECell/FECellData.h"

using std::sort;
using std::make_pair;

/**
 * this is the base class for fe cell(mesh) generator, all the generator should inherit this base class
 */
class FECellGeneratorBase{
public:
    /**
     * constructor
     */
    FECellGeneratorBase(){}

    /**
     * virtual function for the details of different cell(mesh) generation, the child class should
     * offer the implementations
     * @param t_meshtype the type of mesh one want to use for the mesh generation
     * @param t_celldata the fe cell data structure, which should be updated within each cell generator!
     */
    virtual bool generateFECell(const MeshType &t_meshtype,FECellData &t_celldata)=0;

};