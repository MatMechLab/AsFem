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
//+++ Purpose: implement the Neumann type boundary condition
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/NeumannBC.h"

void NeumannBC::ComputeBCValue(const FECalcType &calctype, const double &bcvalue, const vector<double> &params, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &soln, const Vector3d &normal, const LocalShapeFun &shp,const double (&ctan)[3], MatrixXd &localK, VectorXd &localR){
    if(calctype==FECalcType::ComputeResidual){
        ComputeResidual(bcvalue,params,elmtinfo,soln,normal,shp,localR);
    }
    else if(calctype==FECalcType::ComputeJacobian){
        ComputeJacobian(bcvalue,params,elmtinfo,soln,normal,shp,ctan,localK);
    }
    else{
        MessagePrinter::PrintErrorTxt("Unsupported calculation type in NeumannBC class, please check your code");
        MessagePrinter::AsFem_Exit();
    }
}


void NeumannBC::ComputeResidual(const double &bcvalue, const vector<double> &params, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &soln, const Vector3d &normal, const LocalShapeFun &shp, VectorXd &localR){
    // get rid of unused warnings
    if(params.size()||elmtinfo.dt||soln.gpU[0]||normal(1)){}

    localR(1)=bcvalue*shp.test;
}

void NeumannBC::ComputeJacobian(const double &bcvalue, const vector<double> &params, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &soln, const Vector3d &normal, const LocalShapeFun &shp,const double (&ctan)[3], MatrixXd &localK){
    // get rid of unused warnings
    if(bcvalue||params.size()||elmtinfo.dt||soln.gpU[0]||normal(1)||shp.test||ctan[0]){}
    localK(1,1)=0.0; // since bcvalue is a constant, its derivative should be zero!

}

