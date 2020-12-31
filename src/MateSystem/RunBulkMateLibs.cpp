//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.11.30
//+++ Purpose: call all the materials list in MateSystem, include
//+++          built-in materials and UMAT
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/BulkMateSystem.h"

void BulkMateSystem::RunBulkMateLibs(const MateType &imate,const int &mateindex,const int &nDim,
                    const double &t,const double &dt,
                    const Vector3d &gpCoord,
                    const vector<double> &gpU,const vector<double> &gpV,
                    const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                    vector<double> &gpHist,const vector<double> &gpHistOld){
    switch (imate){
    case MateType::NULLMATE:
        return;
        break;
    case MateType::CONSTPOISSONMATE:
        ConstPoissonMaterial(nDim,t,dt,_BulkMateBlockList[mateindex-1]._Parameters,
            gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
        break;
    case MateType::CONSTDIFFUSIONMATE:
        ConstDiffusionMaterial(nDim,t,dt,_BulkMateBlockList[mateindex-1]._Parameters,
                               gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
        break;
    case MateType::CAHNHILLIARDMATE:
        CahnHilliardMaterial(nDim,t,dt,_BulkMateBlockList[mateindex-1]._Parameters,
                             gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
        break;
    case MateType::LINEARELASTICMATE:
        LinearElasticMaterial(nDim,t,dt,_BulkMateBlockList[mateindex-1]._Parameters,
                              gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
        break;
    default:
        MessagePrinter::PrintErrorTxt("unsupported material type in RunBulkMateLibs of MateSystem, please check either your code or your input file");
        MessagePrinter::AsFem_Exit();
        break;
    }
}