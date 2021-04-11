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
//+++ Date   : 2020.11.30
//+++ Purpose: call all the materials list in MateSystem, include
//+++          built-in materials and UMAT
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/BulkMateSystem.h"

void BulkMateSystem::RunBulkMateLibs(const MateType &imate,const int &mateindex,const int &nDim,
                                     const double &t,const double &dt,const Vector3d &gpCoord,
                                     const vector<double> &gpU,const vector<double> &gpUOld,
                                     const vector<double> &gpUdot,const vector<double> &gpUdotOld,
                                     const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradUOld,
                                     const vector<Vector3d> &gpGradUdot,const vector<Vector3d> &gpGradUdotOld){
    switch (imate){
        case MateType::NULLMATE:
            break;
        case MateType::CONSTPOISSONMATE:
            ConstPoissonMaterial::ComputeMaterialProperties(t,dt,nDim,gpCoord,_BulkMateBlockList[mateindex-1]._Parameters,
                                                        gpU,gpUOld,gpUdot,gpUdotOld,
                                                        gpGradU,gpGradUOld,gpGradUdot,gpGradUdotOld,
                                                        _MaterialsOld,_Materials);
            break;
        case MateType::CONSTDIFFUSIONMATE:
            ConstDiffusionMaterial::ComputeMaterialProperties(t,dt,nDim,gpCoord,_BulkMateBlockList[mateindex-1]._Parameters,
                                                          gpU,gpUOld,gpUdot,gpUdotOld,
                                                          gpGradU,gpGradUOld,gpGradUdot,gpGradUdotOld,
                                                          _MaterialsOld,_Materials);
            break;
        case MateType::DOUBLEWELLFREENERGYMATE:
            DoubleWellFreeEnergyMaterial::ComputeMaterialProperties(t,dt,nDim,gpCoord,_BulkMateBlockList[mateindex-1]._Parameters,
                                                                gpU,gpUOld,gpUdot,gpUdotOld,
                                                                gpGradU,gpGradUOld,gpGradUdot,gpGradUdotOld,
                                                                _MaterialsOld,_Materials);
            break;
        case MateType::LINEARELASTICMATE:
            LinearElasticMaterial::ComputeMaterialProperties(t,dt,nDim,gpCoord,_BulkMateBlockList[mateindex-1]._Parameters,
                                                         gpU,gpUOld,gpUdot,gpUdotOld,
                                                         gpGradU,gpGradUOld,gpGradUdot,gpGradUdotOld,
                                                         _MaterialsOld,_Materials);
            break;
        case MateType::INCREMENTSMALLSTRAINMATE:
            IncrementSmallStrainMaterial::ComputeMaterialProperties(t,dt,nDim,gpCoord,_BulkMateBlockList[mateindex-1]._Parameters,
                                                                gpU,gpUOld,gpUdot,gpUdotOld,
                                                                gpGradU,gpGradUOld,gpGradUdot,gpGradUdotOld,
                                                                _MaterialsOld,_Materials);
            break;
        case MateType::NEOHOOKEANMATE:
            NeoHookeanMaterial::ComputeMaterialProperties(t,dt,nDim,gpCoord,_BulkMateBlockList[mateindex-1]._Parameters,
                                                      gpU,gpUOld,gpUdot,gpUdotOld,
                                                      gpGradU,gpGradUOld,gpGradUdot,gpGradUdotOld,
                                                      _MaterialsOld,_Materials);
            break;
        case MateType::PLASTIC1DMATE:
            Plastic1DMaterial::ComputeMaterialProperties(t,dt,nDim,gpCoord,_BulkMateBlockList[mateindex-1]._Parameters,
                                                         gpU,gpUOld,gpUdot,gpUdotOld,
                                                         gpGradU,gpGradUOld,gpGradUdot,gpGradUdotOld,
                                                         _MaterialsOld,_Materials);
            break;
        case MateType::J2PLASTICITYMATE:
            J2PlasticityMaterial::ComputeMaterialProperties(t,dt,nDim,gpCoord,_BulkMateBlockList[mateindex-1]._Parameters,
                                                            gpU,gpUOld,gpUdot,gpUdotOld,
                                                            gpGradU,gpGradUOld,gpGradUdot,gpGradUdotOld,
                                                            _MaterialsOld,_Materials);
            break;
        case MateType::MIEHEFRACTUREMATE:
            MieheFractureMaterial::ComputeMaterialProperties(t,dt,nDim,gpCoord,_BulkMateBlockList[mateindex-1]._Parameters,
                                                         gpU,gpUOld,gpUdot,gpUdotOld,
                                                         gpGradU,gpGradUOld,gpGradUdot,gpGradUdotOld,
                                                         _MaterialsOld,_Materials);
            break;
        default:
            MessagePrinter::PrintErrorTxt("unsupported material type in RunBulkMateLibs of MateSystem, please check either your code or your input file");
            MessagePrinter::AsFem_Exit();
            break;
    }
}