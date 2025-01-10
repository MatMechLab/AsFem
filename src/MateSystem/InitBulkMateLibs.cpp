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


void BulkMateSystem::initBulkMateLibs(const MateType &t_MateType,
                                      const nlohmann::json &Params,
                                      const LocalElmtInfo &ElmtInfo,
                                      const LocalElmtSolution &ElmtSoln){
    switch (t_MateType)
    {
    case MateType::CONSTPOISSONMATE:
        ConstPoissonMaterial::initMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainer);
        break;
    case MateType::POISSON1DBENCHMARKMATE:
        Poisson1DBenchmarkMaterial::initMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainer);
        break;
    case MateType::POISSON2DBENCHMARKMATE:
        Poisson2DBenchmarkMaterial::initMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainer);
        break;
    case MateType::NONLINEARPOISSON2DMATE:
        NonlinearPoisson2DMaterial::initMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainer);
        break;
    case MateType::NONLINEARPOISSON3DMATE:
        NonlinearPoisson3DMaterial::initMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainer);
        break;
    case MateType::CONSTDIFFUSIONNMATE:
        ConstDiffusionMaterial::initMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainer);
        break;
    case MateType::NONLINEARDIFFUSION2DMATE:
        NonlinearDiffusion2DMaterial::initMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainer);
        break;
    //******************************************
    //*** for free energy materials
    //******************************************
    case MateType::DOUBLEWELLMATE:
        DoubleWellPotentialMaterial::initMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainer);
        break;
    case MateType::BINARYMIXMATE:
        BinaryMixtureMaterial::initMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainer);
        break;
    case MateType::KOBAYASHIMATE:
        KobayashiDendriteMaterial::initMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainer);
        break;
    //******************************************
    //*** for mechanics materials
    //******************************************
    case MateType::LINEARELASTICMATE:
        LinearElasticMaterial::initMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainer);
        break;
    case MateType::SAINTVENANTMATE:
        SaintVenantMaterial::initMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainer);
        break;
    case MateType::NEOHOOKEANMATE:
        NeoHookeanMaterial::initMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainer);
        break;
    case MateType::SMALLSTRAINJ2PLASTICITYMATE:
        SmallStrainJ2PlasticityMaterial::initMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainer);
        break;
    case MateType::SMALLSTRAINEXPLAWJ2PLASTICITYMATE:
        SmallStrainExpLawJ2PlasticityMaterial::initMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainer);
        break;
    //******************************************
    //*** for coupled materials
    //******************************************
    case MateType::SMALLSTRAINDIFFUSIONMATE:
        SmallStrainDiffusionMaterial::initMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainer);
        break;
    case MateType::LINEARELASTICFRACMATE:
        LinearElasticFractureMaterial::initMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainer);
        break;
    case MateType::NEOHOOKEANPFFRACTUREMATE:
        NeoHookeanPFFractureMaterial::initMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainer);
        break;
    case MateType::MIEHEFRACTUREMATE:
        MieheFractureMaterial::initMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainer);
        break;
    case MateType::SMALLSTRAINCAHNHILLIARDMATE:
        SmallStrainCahnHilliardMaterial::initMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainer);
        break;
    case MateType::SMALLSTRAINDIFFUSIONJ2MATE:
        SmallStrainDiffusionJ2Material::initMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainer);
        break;
    case MateType::DIFFUSIONACFRACTUREMATE:
        DiffusionACFractureMaterial::initMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainer);
        break;
    //******************************************
    //*** for UMAT
    //******************************************
    case MateType::USER1MATE:
        User1Material::initMaterialProperties(Params,ElmtInfo,ElmtSoln,m_MaterialContainer);
        break;
    default:
        MessagePrinter::printErrorTxt("Unsupported material type in initBulkMateLibs, please check your code");
        MessagePrinter::exitAsFem();
        break;
    }
}