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
//+++ Date   : 2020.12.27
//+++ Purpose: here we can apply the different types of boundary
//+++          condition, and, once more, we only need to focus on
//+++          the calculation on each gauss point !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/BCSystem.h"

void BCSystem::RunBCLibs(const FECalcType &calctype,const BCType &bctype,
        const double &bcvalue,const vector<double> &parameters,
        const Vector3d &normals,const double (&ctan)[3],
        const LocalElmtInfo &elmtinfo,
        const LocalElmtSolution &soln,
        const LocalShapeFun &shp,
        VectorXd &localR,
        MatrixXd &localK){
    switch (bctype){
        case BCType::NEUMANNBC:
            NeumannBC::ComputeBCValue(calctype,bcvalue,parameters,elmtinfo,soln,normals,shp,ctan,localK,localR);
            break;
        case BCType::USER1BC:
            User1BC::ComputeBCValue(calctype,bcvalue,parameters,elmtinfo,soln,normals,shp,ctan,localK,localR);
            break;
        case BCType::USER2BC:
            User2BC::ComputeBCValue(calctype,bcvalue,parameters,elmtinfo,soln,normals,shp,ctan,localK,localR);
            break;
        case BCType::USER3BC:
            User3BC::ComputeBCValue(calctype,bcvalue,parameters,elmtinfo,soln,normals,shp,ctan,localK,localR);
            break;
        case BCType::USER4BC:
            User4BC::ComputeBCValue(calctype,bcvalue,parameters,elmtinfo,soln,normals,shp,ctan,localK,localR);
            break;
        case BCType::USER5BC:
            User5BC::ComputeBCValue(calctype,bcvalue,parameters,elmtinfo,soln,normals,shp,ctan,localK,localR);
            break;
        case BCType::USER6BC:
            User6BC::ComputeBCValue(calctype,bcvalue,parameters,elmtinfo,soln,normals,shp,ctan,localK,localR);
            break;
        case BCType::USER7BC:
            User7BC::ComputeBCValue(calctype,bcvalue,parameters,elmtinfo,soln,normals,shp,ctan,localK,localR);
            break;
        case BCType::USER8BC:
            User8BC::ComputeBCValue(calctype,bcvalue,parameters,elmtinfo,soln,normals,shp,ctan,localK,localR);
            break;
        case BCType::USER9BC:
            User9BC::ComputeBCValue(calctype,bcvalue,parameters,elmtinfo,soln,normals,shp,ctan,localK,localR);
            break;
        case BCType::USER10BC:
            User10BC::ComputeBCValue(calctype,bcvalue,parameters,elmtinfo,soln,normals,shp,ctan,localK,localR);
            break;
   default:
        MessagePrinter::PrintErrorTxt("unsupported boundary condition type in RunBCLibs, please check your input file or your code");
        MessagePrinter::AsFem_Exit();
        break;
    }
}
