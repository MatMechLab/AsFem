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

#include "FECell/FECellPartioner.h"

void FECellPartioner::partFECell(const string &PartTypeName,FECellData &t_CellData){
    if(PartTypeName.find("asfem")!=string::npos){
        FECellDefaultPartitioner::partitionFECell(t_CellData);
    }
    else{
        MessagePrinter::printErrorTxt("Unsupported partition type name("+PartTypeName+"), please check you input file");
        MessagePrinter::exitAsFem();
    }
}