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
//+++ Date   : 2020.11.29
//+++ Purpose: we list all of our elements here for different models
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/BulkElmtSystem.h"

void BulkElmtSystem::RunBulkElmtLibs(const FECalcType &calctype,const ElmtType &elmtytype,
                         const double (&ctan)[2],
                         const LocalElmtInfo &elmtinfo,
                         const LocalElmtSolution &soln,
                         const LocalShapeFun &shp,
                         const Materials &Mate,const Materials &MateOld,
                         ScalarMateType &gpProj,
                         MatrixXd &localK,VectorXd &localR){
    switch (elmtytype){
        case ElmtType::LAPLACEELMT:
            break;
        case ElmtType::POISSONELMT:
            PoissonElmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::TIMEDERIVELMT:
            break;
        case ElmtType::DIFFUSIONELMT:
            DiffusionElmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::CAHNHILLIARDELMT:
            CahnHilliardElmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::MECHANICSELMT:
            MechanicsElmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::MIEHEFRACELMT:
            MieheFractureElmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        default:
            MessagePrinter::PrintErrorTxt("unsupported element type in ElmtSystem, please check your code or your input file");
            MessagePrinter::AsFem_Exit();
            break;
    }
}
