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
//+++ Date   : 2021.04.09
//+++ Purpose: here we define the base class for the FEM calculation
//+++          i.e., the residual and jacobian calculation
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <vector>

#include "Utils/Vector3d.h"
#include "Utils/MathFuns.h"
#include "Utils/VectorXd.h"
#include "Utils/MatrixXd.h"

#include "Utils/RankTwoTensor.h"
#include "Utils/RankFourTensor.h"

#include "FESystem/FECalcType.h"
#include "MateSystem/Materials.h"

#include "Utils/MessagePrinter.h"

class BulkElmtBase{
public:
    virtual void ComputeAll(const FECalcType &calctype,const int &nDim,const int &nNodes,const int &nDofs,
                            const double &t,const double &dt,const double (&ctan)[2],
                            const Vector3d &gpCoords,
                            const vector<double> &gpU,const vector<double> &gpUold,
                            const vector<double> &gpV,const vector<double> &gpVold,
                            const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradUold,
                            const vector<Vector3d> &gpGradV,const vector<Vector3d> &gpGradVold,
                            const double &test,const double &trial,
                            const Vector3d &grad_test,const Vector3d &grad_trial,
                            const Materials &Mate,const Materials &MateOld,
                            map<string,double> &gpProj,
                            MatrixXd &localK,VectorXd &localR)=0;

protected:
    virtual void ComputeResidual(const int &nDim,const int &nNodes,const int &nDofs,
                                 const double &t,const double &dt,const Vector3d &gpCoords,
                                 const vector<double> &gpU,const vector<double> &gpUold,
                                 const vector<double> &gpV,const vector<double> &gpVold,
                                 const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradUold,
                                 const vector<Vector3d> &gpGradV,const vector<Vector3d> &gpGradVold,
                                 const double &test,const Vector3d &grad_test,
                                 const Materials &Mate,const Materials &MateOld,
                                 VectorXd &localR)=0;

    virtual void ComputeJacobian(const int &nDim,const int &nNodes,const int &nDofs,
                                 const double &t,const double &dt,const double (&ctan)[2],
                                 const Vector3d &gpCoords,
                                 const vector<double> &gpU,const vector<double> &gpUold,
                                 const vector<double> &gpV,const vector<double> &gpVold,
                                 const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradUold,
                                 const vector<Vector3d> &gpGradV,const vector<Vector3d> &gpGradVold,
                                 const double &test,const double &trial,
                                 const Vector3d &grad_test,const Vector3d &grad_trial,
                                 const Materials &Mate,const Materials &MateOld,
                                 MatrixXd &localK)=0;

    virtual void ComputeProjection(const int &nDim,const int &nNodes,const int &nDofs,
                                   const double &t,const double &dt,const double (&ctan)[2],
                                   const Vector3d &gpCoords,
                                   const vector<double> &gpU,const vector<double> &gpUold,
                                   const vector<double> &gpV,const vector<double> &gpVold,
                                   const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradUold,
                                   const vector<Vector3d> &gpGradV,const vector<Vector3d> &gpGradVold,
                                   const double &test,const Vector3d &grad_test,
                                   const Materials &Mate,const Materials &MateOld,
                                   map<string,double> &gpProj)=0;

};
