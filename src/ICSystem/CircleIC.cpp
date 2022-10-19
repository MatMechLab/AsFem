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
//+++ Purpose: Implement the circle IC, the given dofs will be
//+++          assigned with two values(in and out circle)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ICSystem/CircleIC.h"

void CircleIC::computeInitialValue(const nlohmann::json &t_params,
                                   const double &icvalue,
                                   const int &dim,
                                   const int &dofs,
                                   const Vector3d &nodecoords,
                                   VectorXd &localU){
    if(icvalue||dim||nodecoords(1)){}

    if(dim!=2){
        MessagePrinter::printErrorTxt("circle ic must be applied in 2d case");
        MessagePrinter::exitAsFem();
    }

    x0=JsonUtils::getValue(t_params,"x0");
    y0=JsonUtils::getValue(t_params,"y0");
    r=JsonUtils::getValue(t_params,"radius");
    innerval=JsonUtils::getValue(t_params,"inside-value");
    outval=JsonUtils::getValue(t_params,"outside-value");

    dist=sqrt((nodecoords(1)-x0)*(nodecoords(1)-x0)+(nodecoords(2)-y0)*(nodecoords(2)-y0));
    if(dist<=r){
        value=innerval;
    }
    else{
        value=outval;
    }

    for(int i=1;i<=dofs;i++){
        localU(i)=value;
    }
}