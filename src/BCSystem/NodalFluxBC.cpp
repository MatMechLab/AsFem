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
//+++ Date   : 2021.10.06
//+++ Purpose: implement the nodal flux boundary condition 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/NodalFluxBC.h"

void NodalFluxBC::ComputeBCValue(const FECalcType &calctype, const double &bcvalue, const vector<double> &params, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &soln,MatrixXd &localK, VectorXd &localR){
    if(calctype==FECalcType::ComputeResidual){
        ComputeResidual(bcvalue,params,elmtinfo,soln,localR);
    }
    else if(calctype==FECalcType::ComputeJacobian){
        ComputeJacobian(bcvalue,params,elmtinfo,soln,localK);
    }
}

void NodalFluxBC::ComputeResidual(const double &bcvalue, const vector<double> &params, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &soln, VectorXd &localR){
    // get rid of unused warnings
    if(bcvalue||params.size()||elmtinfo.dt||soln.gpU[0]){}

    if(params.size()<4){
        MessagePrinter::PrintErrorTxt("for nodal force boundary condition, your params must have at least 4 parameters, i.e. fx fy fz");
        MessagePrinter::AsFem_Exit();
    }
    for(int i=1;i<=elmtinfo.nDofs;i++){
        localR(i)=params[i-1];
    }
}

void NodalFluxBC::ComputeJacobian(const double &bcvalue, const vector<double> &params, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &soln, MatrixXd &localK){
    // get rid of unused warning
    if(bcvalue||params.size()||elmtinfo.dt||soln.gpU[0]){}
    
    if(params.size()<3){
        MessagePrinter::PrintErrorTxt("for nodal force boundary condition, your params must have at least 3 parameters, i.e. fx fy fz");
        MessagePrinter::AsFem_Exit();
    }
    localK.setZero();
}


