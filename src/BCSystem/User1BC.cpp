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
//+++ Date   : 2020.12.27
//+++ Purpose: here we can apply the user-defined-boundary-condition-1
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/BCSystem.h"

void BCSystem::User1BC(const FECalcType &calctype,
                const Vector3d &normals,const double &gpU,const Vector3d &gpGradU,
                const double &bcvalue,
                const double &test,const double &trial,
                const Vector3d &grad_test,const Vector3d &grad_trial,
                double &localK,double &localR){
    if(gpU||grad_test(1)||trial||grad_trial(1)){}
    if(calctype==FECalcType::ComputeResidual){
        localR=(gpGradU*normals-bcvalue)*test;
    }
    else if(calctype==FECalcType::ComputeJacobian){
        localK=(grad_trial*normals)*test;
    }
}