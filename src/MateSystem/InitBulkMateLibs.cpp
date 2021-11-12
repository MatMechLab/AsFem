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


void BulkMateSystem::InitBulkMateLibs(const MateType &imate, const int &mateindex, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln){
    switch (imate){
        case MateType::NULLMATE:
            break;
        case MateType::CONSTPOISSONMATE:
            ConstPoissonMaterial::InitMaterialProperties(_BulkMateBlockList[mateindex-1]._Parameters,elmtinfo,elmtsoln,_Materials);
            break;
        case MateType::CONSTDIFFUSIONMATE:
            ConstDiffusionMaterial::InitMaterialProperties(_BulkMateBlockList[mateindex-1]._Parameters,elmtinfo,elmtsoln,_Materials);
            break;
        case MateType::DOUBLEWELLFREENERGYMATE:
            DoubleWellFreeEnergyMaterial::InitMaterialProperties(_BulkMateBlockList[mateindex-1]._Parameters,elmtinfo,elmtsoln,_Materials);
            break;
        case MateType::IDEALSOLUTIONFREENERGYMATE:
            IdealSolutionFreeEnergyMaterial::InitMaterialProperties(_BulkMateBlockList[mateindex-1]._Parameters,elmtinfo,elmtsoln,_Materials);
            break;
        case MateType::LINEARELASTICMATE:
            LinearElasticMaterial::InitMaterialProperties(_BulkMateBlockList[mateindex-1]._Parameters,elmtinfo,elmtsoln,_Materials);
            break;
        case MateType::INCREMENTSMALLSTRAINMATE:
            IncrementSmallStrainMaterial::InitMaterialProperties(_BulkMateBlockList[mateindex-1]._Parameters,elmtinfo,elmtsoln,_Materials);
            break;
        case MateType::NEOHOOKEANMATE:
            NeoHookeanMaterial::InitMaterialProperties(_BulkMateBlockList[mateindex-1]._Parameters,elmtinfo,elmtsoln,_Materials);
            break;
        case MateType::SAINTVENANTMATE:
            SaintVenantMaterial::InitMaterialProperties(_BulkMateBlockList[mateindex-1]._Parameters,elmtinfo,elmtsoln,_Materials);
            break;
        case MateType::PLASTIC1DMATE:
            Plastic1DMaterial::InitMaterialProperties(_BulkMateBlockList[mateindex-1]._Parameters,elmtinfo,elmtsoln,_Materials);
            break;
        case MateType::J2PLASTICITYMATE:
            J2PlasticityMaterial::InitMaterialProperties(_BulkMateBlockList[mateindex-1]._Parameters,elmtinfo,elmtsoln,_Materials);
            break;
        case MateType::MIEHEFRACTUREMATE:
            MieheFractureMaterial::InitMaterialProperties(_BulkMateBlockList[mateindex-1]._Parameters,elmtinfo,elmtsoln,_Materials);
            break;
        case MateType::STRESSDECOMPOSITIONMATE:
            StressDecompositionMaterial::InitMaterialProperties(_BulkMateBlockList[mateindex-1]._Parameters,elmtinfo,elmtsoln,_Materials);
            break;
        case MateType::NEOHOOKEANPFFRACTUREMATE:
            NeoHookeanPFFractureMaterial::InitMaterialProperties(_BulkMateBlockList[mateindex-1]._Parameters,elmtinfo,elmtsoln,_Materials);
            break;
        case MateType::KOBAYASHIMATE:
            KobayashiMaterial::InitMaterialProperties(_BulkMateBlockList[mateindex-1]._Parameters,elmtinfo,elmtsoln,_Materials);
            break;
        //**************************************************************
        //*** for User-Defined-Materials (UMAT)
        //**************************************************************
        case MateType::USER1MATE:
            User1Material::InitMaterialProperties(_BulkMateBlockList[mateindex-1]._Parameters,elmtinfo,elmtsoln,_Materials);
            break;
        case MateType::USER2MATE:
            User2Material::InitMaterialProperties(_BulkMateBlockList[mateindex-1]._Parameters,elmtinfo,elmtsoln,_Materials);
            break;
        case MateType::USER3MATE:
            User3Material::InitMaterialProperties(_BulkMateBlockList[mateindex-1]._Parameters,elmtinfo,elmtsoln,_Materials);
            break;
        case MateType::USER4MATE:
            User4Material::InitMaterialProperties(_BulkMateBlockList[mateindex-1]._Parameters,elmtinfo,elmtsoln,_Materials);
            break;
        case MateType::USER5MATE:
            User5Material::InitMaterialProperties(_BulkMateBlockList[mateindex-1]._Parameters,elmtinfo,elmtsoln,_Materials);
            break;
        case MateType::USER6MATE:
            User6Material::InitMaterialProperties(_BulkMateBlockList[mateindex-1]._Parameters,elmtinfo,elmtsoln,_Materials);
            break;
        case MateType::USER7MATE:
            User7Material::InitMaterialProperties(_BulkMateBlockList[mateindex-1]._Parameters,elmtinfo,elmtsoln,_Materials);
            break;
        case MateType::USER8MATE:
            User8Material::InitMaterialProperties(_BulkMateBlockList[mateindex-1]._Parameters,elmtinfo,elmtsoln,_Materials);
            break;
        case MateType::USER9MATE:
            User9Material::InitMaterialProperties(_BulkMateBlockList[mateindex-1]._Parameters,elmtinfo,elmtsoln,_Materials);
            break;
        case MateType::USER10MATE:
            User10Material::InitMaterialProperties(_BulkMateBlockList[mateindex-1]._Parameters,elmtinfo,elmtsoln,_Materials);
            break;
        default:
            MessagePrinter::PrintErrorTxt("unsupported material type in InitBulkMateLibs of the MateSystem, please check either your code or your input file");
            MessagePrinter::AsFem_Exit();
            break;
    }
}
