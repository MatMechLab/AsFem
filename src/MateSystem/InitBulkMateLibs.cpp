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
//+++ Purpose: Init all material properties, set the intial status
//+++          of history variables
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/BulkMateSystem.h"

void BulkMateSystem::InitBulkMateLibs(const MateType &imate,const int &mateindex,const int &nDim,const Vector3d &gpCoord,
                                      const vector<double> &gpU,const vector<double> &gpUdot,
                                      const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradUdot){
    switch (imate){
        case MateType::NULLMATE:
            break;
        case MateType::CONSTPOISSONMATE:
            ConstPoissonMaterial::InitMaterialProperties(nDim,gpCoord,_BulkMateBlockList[mateindex-1]._Parameters,
                                                         gpU,gpUdot,gpGradU,gpGradUdot,_Materials);
            break;
        case MateType::CONSTDIFFUSIONMATE:
            ConstDiffusionMaterial::InitMaterialProperties(nDim,gpCoord,_BulkMateBlockList[mateindex-1]._Parameters,
                                                           gpU,gpUdot,gpGradU,gpGradUdot,_Materials);
            break;
        case MateType::DOUBLEWELLFREENERGYMATE:
            DoubleWellFreeEnergyMaterial::InitMaterialProperties(nDim,gpCoord,_BulkMateBlockList[mateindex-1]._Parameters,
                                                                 gpU,gpUdot,gpGradU,gpGradUdot,_Materials);
            break;
        case MateType::LINEARELASTICMATE:
            LinearElasticMaterial::InitMaterialProperties(nDim,gpCoord,_BulkMateBlockList[mateindex-1]._Parameters,
                                                          gpU,gpUdot,gpGradU,gpGradUdot,_Materials);
            break;
        case MateType::INCREMENTSMALLSTRAINMATE:
            IncrementSmallStrainMaterial::InitMaterialProperties(nDim,gpCoord,_BulkMateBlockList[mateindex-1]._Parameters,
                                                                 gpU,gpUdot,gpGradU,gpGradUdot,_Materials);
            break;
        case MateType::NEOHOOKEANMATE:
            NeoHookeanMaterial::InitMaterialProperties(nDim,gpCoord,_BulkMateBlockList[mateindex-1]._Parameters,
                                                       gpU,gpUdot,gpGradU,gpGradUdot,_Materials);
            break;
        case MateType::MIEHEFRACTUREMATE:
            MieheFractureMaterial::InitMaterialProperties(nDim,gpCoord,_BulkMateBlockList[mateindex-1]._Parameters,
                                                          gpU,gpUdot,gpGradU,gpGradUdot,_Materials);
            break;
        default:
            MessagePrinter::PrintErrorTxt("unsupported material type in RunBulkMateLibs of MateSystem, please check either your code or your input file");
            MessagePrinter::AsFem_Exit();
            break;
    }
}