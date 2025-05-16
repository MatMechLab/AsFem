//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
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


void BulkMateSystem::runBulkMateLibs(const MateType &t_MateType,
                                     const nlohmann::json &Params,
                                     const LocalElmtInfo &ElmtInfo,
                                     const LocalElmtSolution &ElmtSoln){
    switch (t_MateType)
    {
    case MateType::CONSTPOISSONMATE:
        ConstPoissonMaterial::computeMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainerOld,m_MaterialContainer);
        break;
    case MateType::POISSON1DBENCHMARKMATE:
        Poisson1DBenchmarkMaterial::computeMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainerOld,m_MaterialContainer);
        break;
    case MateType::POISSON2DBENCHMARKMATE:
        Poisson2DBenchmarkMaterial::computeMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainerOld,m_MaterialContainer);
        break;
    case MateType::NONLINEARPOISSON2DMATE:
        NonlinearPoisson2DMaterial::computeMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainerOld,m_MaterialContainer);
        break;
    case MateType::NONLINEARPOISSON3DMATE:
        NonlinearPoisson3DMaterial::computeMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainerOld,m_MaterialContainer);
        break;
    case MateType::CONSTDIFFUSIONNMATE:
        ConstDiffusionMaterial::computeMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainerOld,m_MaterialContainer);
        break;
    case MateType::NONLINEARDIFFUSION2DMATE:
        NonlinearDiffusion2DMaterial::computeMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainerOld,m_MaterialContainer);
        break;
    //******************************************
    //*** for free energy materials
    //******************************************
    case MateType::DOUBLEWELLMATE:
        DoubleWellPotentialMaterial::computeMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainerOld,m_MaterialContainer);
        break;
    case MateType::BINARYMIXMATE:
        BinaryMixtureMaterial::computeMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainerOld,m_MaterialContainer);
        break;
    case MateType::KOBAYASHIMATE:
        KobayashiDendriteMaterial::computeMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainerOld,m_MaterialContainer);
        break;
    //******************************************
    //*** for mechanics materials
    //******************************************
    case MateType::LINEARELASTICMATE:
        LinearElasticMaterial::computeMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainerOld,m_MaterialContainer);
        break;
    case MateType::SAINTVENANTMATE:
        SaintVenantMaterial::computeMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainerOld,m_MaterialContainer);
        break;
    case MateType::NEOHOOKEANMATE:
        NeoHookeanMaterial::computeMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainerOld,m_MaterialContainer);
        break;
    case MateType::SMALLSTRAINJ2PLASTICITYMATE:
        SmallStrainJ2PlasticityMaterial::computeMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainerOld,m_MaterialContainer);
        break;
    case MateType::SMALLSTRAINEXPLAWJ2PLASTICITYMATE:
        SmallStrainExpLawJ2PlasticityMaterial::computeMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainerOld,m_MaterialContainer);
        break;
    //******************************************
    //*** for coupled materials
    //******************************************
    case MateType::SMALLSTRAINDIFFUSIONMATE:
        SmallStrainDiffusionMaterial::computeMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainerOld,m_MaterialContainer);
        break;
    case MateType::LINEARELASTICFRACMATE:
        LinearElasticFractureMaterial::computeMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainerOld,m_MaterialContainer);
        break;
    case MateType::NEOHOOKEANPFFRACTUREMATE:
        NeoHookeanPFFractureMaterial::computeMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainerOld,m_MaterialContainer);
        break;
    case MateType::MIEHEFRACTUREMATE:
        MieheFractureMaterial::computeMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainerOld,m_MaterialContainer);
        break;
    case MateType::SMALLSTRAINCAHNHILLIARDMATE:
        SmallStrainCahnHilliardMaterial::computeMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainerOld,m_MaterialContainer);
        break;
    case MateType::SMALLSTRAINDIFFUSIONJ2MATE:
        SmallStrainDiffusionJ2Material::computeMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainerOld,m_MaterialContainer);
        break;
    case MateType::DIFFUSIONACFRACTUREMATE:
        DiffusionACFractureMaterial::computeMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainerOld,m_MaterialContainer);
        break;
    //******************************************
    //*** for User-Defined-Material (UMAT)
    //******************************************
    case MateType::USER1MATE:
        User1Material::computeMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainerOld,m_MaterialContainer);
        break;
    default:
        MessagePrinter::printErrorTxt("Unsupported material type in runBulkMateLibs, please check your code");
        MessagePrinter::exitAsFem();
        break;
    }
}