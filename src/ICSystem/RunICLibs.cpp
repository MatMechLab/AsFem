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
//+++ Date   : 2020.07.10
//+++ Purpose: Call the different initial conditions
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ICSystem/ICSystem.h"

void ICSystem::runICLibs(const ICType &t_type,
                         const nlohmann::json &t_params,
                         const double &icvalue,
                         const int &dim,
                         const int &dofs,
                         const Vector3d &nodecoords,
                         VectorXd &localU){
    switch (t_type)
    {
    case ICType::CONSTIC:
        ConstantIC::computeInitialValue(t_params,icvalue,dim,dofs,nodecoords,localU);
        break;
    case ICType::RANDOMIC:
        RandomIC::computeInitialValue(t_params,icvalue,dim,dofs,nodecoords,localU);
        break;
    case ICType::CIRCLEIC:
        CircleIC::computeInitialValue(t_params,icvalue,dim,dofs,nodecoords,localU);
        break;
    default:
        MessagePrinter::printErrorTxt("Unsupported initial condition type in RunICLibs.cpp, please check your code or your input file");
        MessagePrinter::exitAsFem();
        break;
    }
}