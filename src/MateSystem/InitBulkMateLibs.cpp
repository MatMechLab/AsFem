//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
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


void BulkMateSystem::initBulkMateLibs(const MateType &t_matetype,const nlohmann::json &t_params,
                                      const LocalElmtInfo &t_elmtinfo,const LocalElmtSolution &t_elmtsoln){
    switch (t_matetype)
    {
    case MateType::CONSTPOISSONMATE:
        ConstPoissonMaterial::initMaterialProperties(t_params,t_elmtinfo,t_elmtsoln,m_materialcontainer);
        break;
    case MateType::POISSON1DBENCHMARKMATE:
        Poisson1DBenchmarkMaterial::initMaterialProperties(t_params,t_elmtinfo,t_elmtsoln,m_materialcontainer);
        break;
    case MateType::POISSON2DBENCHMARKMATE:
        Poisson2DBenchmarkMaterial::initMaterialProperties(t_params,t_elmtinfo,t_elmtsoln,m_materialcontainer);
        break;
    case MateType::NONLINEARPOISSON2DMATE:
        NonlinearPoisson2DMaterial::initMaterialProperties(t_params,t_elmtinfo,t_elmtsoln,m_materialcontainer);
        break;
    case MateType::NONLINEARPOISSON3DMATE:
        NonlinearPoisson3DMaterial::initMaterialProperties(t_params,t_elmtinfo,t_elmtsoln,m_materialcontainer);
        break;
    case MateType::CONSTDIFFUSIONNMATE:
        ConstDiffusionMaterial::initMaterialProperties(t_params,t_elmtinfo,t_elmtsoln,m_materialcontainer);
        break;
    case MateType::NONLINEARDIFFUSION2DMATE:
        NonlinearDiffusion2DMaterial::initMaterialProperties(t_params,t_elmtinfo,t_elmtsoln,m_materialcontainer);
        break;
    //******************************************
    //*** for free energy materials
    //******************************************
    case MateType::DOUBLEWELLMATE:
        DoubleWellPotentialMaterial::initMaterialProperties(t_params,t_elmtinfo,t_elmtsoln,m_materialcontainer);
        break;
    case MateType::BINARYMIXMATE:
        BinaryMixtureMaterial::initMaterialProperties(t_params,t_elmtinfo,t_elmtsoln,m_materialcontainer);
        break;
    case MateType::KOBAYASHIMATE:
        KobayashiDendriteMaterial::initMaterialProperties(t_params,t_elmtinfo,t_elmtsoln,m_materialcontainer);
        break;
    //******************************************
    //*** for mechanics materials
    //******************************************
    case MateType::LINEARELASTICMATE:
        LinearElasticMaterial::initMaterialProperties(t_params,t_elmtinfo,t_elmtsoln,m_materialcontainer);
        break;
    case MateType::SAINTVENANTMATE:
        SaintVenantMaterial::initMaterialProperties(t_params,t_elmtinfo,t_elmtsoln,m_materialcontainer);
        break;
    case MateType::NEOHOOKEANMATE:
        NeoHookeanMaterial::initMaterialProperties(t_params,t_elmtinfo,t_elmtsoln,m_materialcontainer);
        break;
    case MateType::SMALLSTRAINJ2PLASTICITYMATE:
        SmallStrainJ2PlasticityMaterial::initMaterialProperties(t_params,t_elmtinfo,t_elmtsoln,m_materialcontainer);
        break;
    case MateType::SMALLSTRAINEXPLAWJ2PLASTICITYMATE:
        SmallStrainExpLawJ2PlasticityMaterial::initMaterialProperties(t_params,t_elmtinfo,t_elmtsoln,m_materialcontainer);
        break;
    //******************************************
    //*** for coupled materials
    //******************************************
    case MateType::SMALLSTRAINDIFFUSIONMATE:
        SmallStrainDiffusionMaterial::initMaterialProperties(t_params,t_elmtinfo,t_elmtsoln,m_materialcontainer);
        break;
    case MateType::LINEARELASTICFRACMATE:
        LinearElasticFractureMaterial::initMaterialProperties(t_params,t_elmtinfo,t_elmtsoln,m_materialcontainer);
        break;
    case MateType::NEOHOOKEANPFFRACTUREMATE:
        NeoHookeanPFFractureMaterial::initMaterialProperties(t_params,t_elmtinfo,t_elmtsoln,m_materialcontainer);
        break;
    case MateType::MIEHEFRACTUREMATE:
        MieheFractureMaterial::initMaterialProperties(t_params,t_elmtinfo,t_elmtsoln,m_materialcontainer);
        break;
    case MateType::SMALLSTRAINCAHNHILLIARDMATE:
        SmallStrainCahnHilliardMaterial::initMaterialProperties(t_params,t_elmtinfo,t_elmtsoln,m_materialcontainer);
        break;
    case MateType::SMALLSTRAINDIFFUSIONJ2MATE:
        SmallStrainDiffusionJ2Material::initMaterialProperties(t_params,t_elmtinfo,t_elmtsoln,m_materialcontainer);
        break;
    case MateType::DIFFUSIONACFRACTUREMATE:
        DiffusionACFractureMaterial::initMaterialProperties(t_params,t_elmtinfo,t_elmtsoln,m_materialcontainer);
        break;
    default:
        MessagePrinter::printErrorTxt("Unsupported material type in initBulkMateLibs, please check your code");
        MessagePrinter::exitAsFem();
        break;
    }
}