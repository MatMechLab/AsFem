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
//+++ Date   : 2021.10.17
//+++ Purpose: implement the user-defined-integrated type boundary
//+++          condition, if one wants to use dirichlet type bc,
//+++          please call UserXXDirichletBC !!!
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/User1BC.h"

void User1BC::ComputeBCValue(const FECalcType &calctype, const double &bcvalue, const vector<double> &params, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &soln, const Vector3d &normal, const LocalShapeFun &shp,const double (&ctan)[3], MatrixXd &localK, VectorXd &localR){
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


void User1BC::ComputeResidual(const double &bcvalue, const vector<double> &params, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &soln, const Vector3d &normal, const LocalShapeFun &shp, VectorXd &localR){
    // get rid of unused warnings
    if(params.size()||elmtinfo.dt||soln.gpU[0]||normal(1)){}
    localR.setZero();
    localR(1)=bcvalue*shp.test;// here the bcvalue is your "flux"!
}

void User1BC::ComputeJacobian(const double &bcvalue, const vector<double> &params, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &soln, const Vector3d &normal, const LocalShapeFun &shp,const double (&ctan)[3], MatrixXd &localK){
    // get rid of unused warnings
    if(bcvalue||params.size()||elmtinfo.dt||soln.gpU[0]||normal(1)||shp.test||ctan[0]){}
    localK.setZero();
    localK(1,1)=0.0; // since bcvalue is a constant, its derivative should be zero!
}

