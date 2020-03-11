//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "FESystem/FESystem.h"

void FESystem::AssembleLocalToGlobal(const int &isw,const int &ndofs,vector<int> &elDofs,
                            vector<double> &localK,vector<double> &localR,
                            Mat &AMATRIX,Vec &RHS){
    if(isw==3||isw==6){
        VecSetValues(RHS,ndofs,elDofs.data(),localR.data(),ADD_VALUES);
        if(isw==6){
            MatSetValues(AMATRIX,ndofs,elDofs.data(),ndofs,elDofs.data(),localK.data(),ADD_VALUES);
        }
    }
}