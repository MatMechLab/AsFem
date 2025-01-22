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
//+++ Date    : 2024.08.04
//+++ Function: the fe cell partitioner used by AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "FECell/FECellDefaultPartitioner.h"
#include "FECell/FECellMETISPartitioner.h"

class FECellPartioner:public FECellDefaultPartitioner,
                      public FECellMETISPartitioner{
public:
    /**
     * the fe cell partition
     * @param PartTypeName the name of the partition algorithm
     * @param t_CellData the fe cell data
     */
    void partFECell(const string &PartTypeName,FECellData &t_CellData);
};