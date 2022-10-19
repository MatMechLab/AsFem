//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2022.08.10
//+++ Purpose: Implement the random IC, the given dofs will be
//+++          assigned with a random value
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ICSystem/RandomIC.h"

void RandomIC::computeInitialValue(const nlohmann::json &t_params,
                                   const double &icvalue,
                                   const int &dim,
                                   const int &dofs,
                                   const Vector3d &nodecoords,
                                   VectorXd &localU){
    if(t_params.size()||dim||nodecoords(1)){}


    for(int i=1;i<=dofs;i++){
        localU(i)=icvalue;
    }
}