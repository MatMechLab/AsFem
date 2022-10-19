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
//+++ Date   : 2020.12.27
//+++ Purpose: here we can apply the different types of boundary
//+++          condition (integrated type)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/BCSystem.h"

void BCSystem::runBCLibs(const FECalcType &calctype,
                         const BCType &bctype,
                         const double &bcvalue,
                         const double (&ctan)[3],
                         const nlohmann::json &json,
                         const Vector3d &normal,
                         const LocalElmtInfo &elmtinfo,
                         const LocalElmtSolution &elmtsoln,
                         const LocalShapeFun &shp,
                         MatrixXd &localK,
                         VectorXd &localR){
    switch (bctype)
    {
    case BCType::NEUMANNBC:
        NeumannBC::computeBCValue(calctype,bcvalue,json,elmtinfo,elmtsoln,normal,shp,ctan,localK,localR);
        break;
    case BCType::PRESSUREBC:
        PressureBC::computeBCValue(calctype,bcvalue,json,elmtinfo,elmtsoln,normal,shp,ctan,localK,localR);
        break;
    case BCType::TRACTIONBC:
        TractionBC::computeBCValue(calctype,bcvalue,json,elmtinfo,elmtsoln,normal,shp,ctan,localK,localR);
        break;
    default:
        MessagePrinter::printErrorTxt("unsupported boundary condition type in RunBCLibs.cpp, please check your input file or your code");
        MessagePrinter::exitAsFem();
        break;
    }
}