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
//+++ Purpose: Define the boundary condition system in AsFem
//+++          Here we can apply:
//+++               1) constant value ic
//+++               2) random value ic
//+++               3) other type or user defined ic
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ICSystem/ICSystem.h"


void ICSystem::ApplyIC(const Mesh &mesh,const DofHandler &dofHandler,Vec &U){
    vector<string> domainlist;
    vector<double> parameters;
    int DofIndex;
    for(auto it:_ICBlockList){
        domainlist=it._DomainNameList;
        parameters=it._Parameters;
        DofIndex=it._DofID;
        RunICLibs(it._ICType,DofIndex,parameters,domainlist,mesh,dofHandler,U);
    }
    VecAssemblyBegin(U);
    VecAssemblyEnd(U);
}