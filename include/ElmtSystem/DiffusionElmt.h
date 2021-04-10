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
//+++ Date   : 2021.04.10
//+++ Purpose: implement the residual and jacobian for general
//+++          diffusion equation:
//+++          dc/dt=div(D*grad(c))
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "ElmtSystem/BulkElmtBase.h"

class DiffusionElmt:public BulkElmtBase{
public:
    virtual void ComputeAll(const FECalcType &calctype, const int &nDim, const int &nNodes,
                            const int &nDofs, const double &t, const double &dt, const double (&ctan)[2],
                            const Vector3d &gpCoords, const vector<double> &gpU,
                            const vector<double> &gpUold, const vector<double> &gpV,
                            const vector<double> &gpVold, const vector<Vector3d> &gpGradU,
                            const vector<Vector3d> &gpGradUold, const vector<Vector3d> &gpGradV,
                            const vector<Vector3d> &gpGradVold, const double &test, const double &trial,
                            const Vector3d &grad_test, const Vector3d &grad_trial, const Materials &Mate,
                            const Materials &MateOld, map<string, double> &gpProj, MatrixXd &localK,
                            VectorXd &localR) override;

private:
    virtual void ComputeResidual(const int &nDim, const int &nNodes, const int &nDofs, const double &t,
                                 const double &dt, const Vector3d &gpCoords, const vector<double> &gpU,
                                 const vector<double> &gpUold, const vector<double> &gpV,
                                 const vector<double> &gpVold, const vector<Vector3d> &gpGradU,
                                 const vector<Vector3d> &gpGradUold, const vector<Vector3d> &gpGradV,
                                 const vector<Vector3d> &gpGradVold, const double &test,
                                 const Vector3d &grad_test, const Materials &Mate,
                                 const Materials &MateOld, VectorXd &localR) override;

    virtual void ComputeJacobian(const int &nDim, const int &nNodes, const int &nDofs, const double &t,
                                 const double &dt, const double (&ctan)[2], const Vector3d &gpCoords,
                                 const vector<double> &gpU, const vector<double> &gpUold,
                                 const vector<double> &gpV, const vector<double> &gpVold,
                                 const vector<Vector3d> &gpGradU, const vector<Vector3d> &gpGradUold,
                                 const vector<Vector3d> &gpGradV, const vector<Vector3d> &gpGradVold,
                                 const double &test, const double &trial, const Vector3d &grad_test,
                                 const Vector3d &grad_trial, const Materials &Mate,
                                 const Materials &MateOld, MatrixXd &localK) override;

    virtual void ComputeProjection(const int &nDim, const int &nNodes, const int &nDofs, const double &t,
                                   const double &dt, const double (&ctan)[2], const Vector3d &gpCoords,
                                   const vector<double> &gpU, const vector<double> &gpUold,
                                   const vector<double> &gpV, const vector<double> &gpVold,
                                   const vector<Vector3d> &gpGradU, const vector<Vector3d> &gpGradUold,
                                   const vector<Vector3d> &gpGradV, const vector<Vector3d> &gpGradVold,
                                   const double &test, const Vector3d &grad_test, const Materials &Mate,
                                   const Materials &MateOld, map<string, double> &gpProj) override;

};