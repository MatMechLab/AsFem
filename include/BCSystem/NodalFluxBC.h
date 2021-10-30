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
//+++ Date   : 2021.10.06
//+++ Purpose: implement the nodal flux boundary condition 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "BCSystem/NodalBCBase.h"


class NodalFluxBC:public NodalBCBase{
public:
    void ComputeBCValue(const FECalcType &calctype, const double &bcvalue, const vector<double> &params, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &soln, MatrixXd &localK, VectorXd &localR) override;

private:
    void ComputeResidual(const double &bcvalue, const vector<double> &params, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &soln, VectorXd &localR) override;

    void ComputeJacobian(const double &bcvalue, const vector<double> &params, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &soln, MatrixXd &localK) override;

};
