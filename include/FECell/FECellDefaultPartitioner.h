//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author  : Yang Bai
//+++ Date    : 2024.08.01
//+++ Function: the built-in fe cell partitioner
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "FECell/FECellPartitionerBase.h"


class FECellDefaultPartitioner:public FECellPartitionerBase{
public:
    /**
     * partition the fe cell based on mpi rank
     * @param t_CellData the fecell data 
     */
    virtual void partitionFECell(FECellData &t_CellData) override;
};