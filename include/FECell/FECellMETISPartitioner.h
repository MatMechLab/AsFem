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
//+++ Date    : 2025.01.22
//+++ Function: the fe cell partitioner based on METIS
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "FECell/FECellPartitionerBase.h"


class FECellMETISPartitioner:public FECellPartitionerBase{
public:
    /**
     * partition the fe cell based on mpi rank
     * @param t_CellData the fecell data
     */
    virtual void partitionFECell(FECellData &t_CellData) override;
};