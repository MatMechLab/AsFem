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
//+++ Function: base class for fe cell partion based on mpi
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once


#include "Utils/MessagePrinter.h"
#include "MPIUtils/MPIDataBus.h"
#include "FECell/FECellData.h"


class FECellPartitionerBase{
public:
    /**
     * partition the fe cell based on mpi rank
     * @param t_celldata the fecell data 
     */
    virtual void partitionFECell(FECellData &t_celldata)=0;
};