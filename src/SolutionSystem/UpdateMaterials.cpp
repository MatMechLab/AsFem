//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.10.31
//+++ Purpose: update the material properties calculated from MateSystem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "SolutionSystem/SolutionSystem.h"

void SolutionSystem::UpdateMaterials(){
    for(int i=0;i<static_cast<int>(_ScalarMaterials.size());i++){
        _ScalarMaterialsOld[i]=_ScalarMaterials[i];
    }
    for(int i=0;i<static_cast<int>(_VectorMaterials.size());i++){
        _ScalarMaterialsOld[i]=_ScalarMaterials[i];
    }
    for(int i=0;i<static_cast<int>(_Rank2TensorMaterials.size());i++){
        _Rank2TensorMaterialsOld[i]=_Rank2TensorMaterials[i];
    }
    for(int i=0;i<static_cast<int>(_Rank4TensorMaterials.size());i++){
        _Rank4TensorMaterialsOld[i]=_Rank4TensorMaterials[i];
    }
}
