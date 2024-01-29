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
    if(meshtype==MeshType::HEX8||
       meshtype==MeshType::HEX20||
       meshtype==MeshType::HEX27){
        return Lagrange3DMeshCellGenerator::generateFECell(meshtype,meshdata);
    }
    else{
        return false;
    }
    return true;
}
