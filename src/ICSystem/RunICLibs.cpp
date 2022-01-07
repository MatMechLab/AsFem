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
//+++ Purpose: Apply different types of ICs in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ICSystem/ICSystem.h"

void ICSystem::RunICLibs(const ICType &ictype,const int &DofIndex,const vector<double> &Parameters,const vector<string> &DomainList,
                         const Mesh &mesh,const DofHandler &dofHandler,Vec &U){
    switch (ictype){
        case ICType::NULLIC:
            break;
        case ICType::CONSTIC:
            ApplyConstantIC(DofIndex,Parameters,DomainList,mesh,dofHandler,U);
            break;
        case ICType::RANDOMIC:
            ApplyRandomIC(DofIndex,Parameters,DomainList,mesh,dofHandler,U);
            break;
        case ICType::CIRCLEIC:
            ApplyCircleIC(DofIndex,Parameters,DomainList,mesh,dofHandler,U);
            break;
        case ICType::SMOOTHCIRCLEIC:
            ApplySmoothCircleIC(DofIndex,Parameters,DomainList,mesh,dofHandler,U);
            break;
        case ICType::RECTANGLEIC:
            ApplyRectangleIC(DofIndex,Parameters,DomainList,mesh,dofHandler,U);
            break;
        case ICType::CUBICIC:
            ApplyCubicIC(DofIndex,Parameters,DomainList,mesh,dofHandler,U);
            break;
        case ICType::SPHERICALIC:
            ApplySphereIC(DofIndex,Parameters,DomainList,mesh,dofHandler,U);
            break;
        case ICType::USER1IC:
            User1IC(DofIndex,Parameters,DomainList,mesh,dofHandler,U);
            break;
        default:
            MessagePrinter::PrintErrorTxt("unsupported initial condition type in ICSystem, please check your input file");
            MessagePrinter::AsFem_Exit();
            break;
    }
}
