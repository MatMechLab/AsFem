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
                                     const int &nDim,const int &nNodes,const int &nDofs,
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
                                     MatrixXd &localK,VectorXd &localR){
    switch (elmtytype){
    case ElmtType::LAPLACEELMT:
        break;
    case ElmtType::POISSONELMT:
        break;
    case ElmtType::TIMEDERIVELMT:
        break;
    case ElmtType::DIFFUSIONELMT:
        break;
    case ElmtType::CAHNHILLIARDELMT:
        CahnHilliardElmt::ComputeAll(calctype,nDim,nNodes,nDofs,t,dt,ctan,
                                     gpCoords,gpU,gpUold,gpV,gpVold,
                                     gpGradU,gpGradUold,gpGradV,gpGradVold,
                                     test,trial,grad_test,grad_trial,
                                     Mate,MateOld,gpProj,localK,localR);
        break;
    case ElmtType::MECHANICSELMT:
        MechanicsElmt::ComputeAll(calctype,nDim,nNodes,nDofs,t,dt,ctan,
                                  gpCoords,gpU,gpUold,gpV,gpVold,
                                  gpGradU,gpGradUold,gpGradV,gpGradVold,
                                  test,trial,grad_test,grad_trial,
                                  Mate,MateOld,gpProj,localK,localR);
        break;
    case ElmtType::MIEHEFRACELMT:
        break;
    default:
        MessagePrinter::PrintErrorTxt("unsupported element type in ElmtSystem, please check your code or your input file");
        MessagePrinter::AsFem_Exit();
        break;
    }
}