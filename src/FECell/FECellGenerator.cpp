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

#include "FECell/FECellGenerator.h"

FECellGenerator::FECellGenerator(){

}
//********************************************
bool FECellGenerator::createFEMeshCell(const MeshType &meshtype,FECellData &meshdata){
    switch (meshtype)
    {
    case MeshType::HEX8:
        return Lagrange3DHex8MeshCellGenerator::generateFECell(meshdata);
        break;
    case MeshType::HEX20:
        return Lagrange3DHex20MeshCellGenerator::generateFECell(meshdata);
        break;
    case MeshType::HEX27:
        return Lagrange3DHex27MeshCellGenerator::generateFECell(meshdata);
        break;
    default:
        MessagePrinter::printErrorTxt("Unsupported meshtype for FEMeshCell generation");
        MessagePrinter::exitAsFem();
        break;
    }
    return true;
}
