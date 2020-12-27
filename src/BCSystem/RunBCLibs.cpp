//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.12.27
//+++ Purpose: here we can apply the different types of boundary
//+++          condition, and, once more, we only need to focus on
//+++          the calculation on each gauss point !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/BCSystem.h"

void BCSystem::RunBCLibs(const BCType bctype,const FECalcType &calctype,
                const Vector3d &normals,const double &gpU,const Vector3d &gpGradU,
                const double &bcvalue,
                const double &test,const double &trial,
                const Vector3d &grad_test,const Vector3d &grad_trial,
                double &localK,double &localR){
    switch (bctype)
    {
    case BCType::USER1BC:
        User1BC(calctype,normals,gpU,gpGradU,bcvalue,test,trial,grad_test,grad_trial,localK,localR);
        break;
    default:
        MessagePrinter::PrintErrorTxt("unsupported boundary condition type");
        MessagePrinter::AsFem_Exit();
        break;
    }
}