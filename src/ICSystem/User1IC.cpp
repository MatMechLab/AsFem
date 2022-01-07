//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.12.30
//+++ Purpose: Apply the user-defined IC to U vector
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ICSystem/ICSystem.h"

void ICSystem::User1IC(const int &DofIndex, const vector<double> &Parameters, const vector<string> &DomainList,
                       const Mesh &mesh, const DofHandler &dofHandler, Vec &U){
    if(DofIndex||Parameters.size()||DomainList.size()||
    mesh.GetDim()||dofHandler.GetTotalDofsNum()){}
    int i;
    VecGetSize(U,&i);
    if(i){return;}
}
