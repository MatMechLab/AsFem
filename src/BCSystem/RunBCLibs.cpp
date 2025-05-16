//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
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

void BCSystem::runBCLibs(const FECalcType &CalcType,
                         const BCType &t_BCType,
                         const double &BCValue,
                         const double (&Ctan)[3],
                         const nlohmann::json &Params,
                         const Vector3d &Normal,
                         const LocalElmtInfo &ElmtInfo,
                         const LocalElmtSolution &ElmtSoln,
                         const LocalShapeFun &Shp,
                         MatrixXd &LocalK,
                         VectorXd &LocalR){
    switch (t_BCType)
    {
    case BCType::NEUMANNBC:
        NeumannBC::computeBCValue(CalcType,BCValue,Params,ElmtInfo,ElmtSoln,Normal,Shp,Ctan,LocalK,LocalR);
        break;
    case BCType::PRESSUREBC:
        PressureBC::computeBCValue(CalcType,BCValue,Params,ElmtInfo,ElmtSoln,Normal,Shp,Ctan,LocalK,LocalR);
        break;
    case BCType::TRACTIONBC:
        TractionBC::computeBCValue(CalcType,BCValue,Params,ElmtInfo,ElmtSoln,Normal,Shp,Ctan,LocalK,LocalR);
        break;
    default:
        MessagePrinter::printErrorTxt("unsupported boundary condition type in RunBCLibs.cpp, please check your input file or your code");
        MessagePrinter::exitAsFem();
        break;
    }
}