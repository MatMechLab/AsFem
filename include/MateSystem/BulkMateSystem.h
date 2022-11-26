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
//+++ Purpose: Implement the materials system for the bulk element
//+++          in AsFem. It is different from another one, namely
//+++          the interface material system
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#pragma once

/**
 * For AsFem's own headers
 */
#include "MateSystem/MateType.h"
#include "MateSystem/MaterialsContainer.h"

/**
 * For built-in and user-defined materials (UMAT)
 */
#include "MateSystem/ConstPoissonMaterial.h"
#include "MateSystem/Poisson1DBenchmarkMaterial.h"
#include "MateSystem/Poisson2DBenchmarkMaterial.h"
#include "MateSystem/NonlinearPoisson2DMaterial.h"
#include "MateSystem/NonlinearPoisson3DMaterial.h"
#include "MateSystem/ConstDiffusionMaterial.h"
#include "MateSystem/NonlinearDiffusion2DMaterial.h"

#include "MateSystem/DoubleWellPotentialMaterial.h"
#include "MateSystem/BinaryMixtureMaterial.h"
#include "MateSystem/KobayashiDendriteMaterial.h"

#include "MateSystem/LinearElasticMaterial.h"
#include "MateSystem/SaintVenantMaterial.h"
#include "MateSystem/NeoHookeanMaterial.h"
#include "MateSystem/SmallStrainJ2PlasticityMaterial.h"
#include "MateSystem/SmallStrainExpLawJ2PlasticityMaterial.h"

#include "MateSystem/SmallStrainDiffusionMaterial.h"
#include "MateSystem/LinearElasticFractureMaterial.h"
#include "MateSystem/MieheFractureMaterial.h"
#include "MateSystem/SmallStrainCahnHilliardMaterial.h"
#include "MateSystem/NeoHookeanPFFractureMaterial.h"
#include "MateSystem/SmallStrainDiffusionJ2Material.h"
#include "MateSystem/DiffusionACFractureMaterial.h"


/**
 * This class implement the materials calculation for the bulk element in AsFem.
 */
class BulkMateSystem:public ConstPoissonMaterial,
                     public Poisson1DBenchmarkMaterial,
                     public Poisson2DBenchmarkMaterial,
                     public NonlinearPoisson2DMaterial,
                     public NonlinearPoisson3DMaterial,
                     public ConstDiffusionMaterial,
                     public NonlinearDiffusion2DMaterial,
                     // for free energy materials
                     public DoubleWellPotentialMaterial,
                     public BinaryMixtureMaterial,
                     public KobayashiDendriteMaterial,
                     // for mechanics materials
                     public LinearElasticMaterial,
                     public SaintVenantMaterial,
                     public NeoHookeanMaterial,
                     public SmallStrainJ2PlasticityMaterial,
                     public SmallStrainExpLawJ2PlasticityMaterial,
                     // for coupled materials
                     public SmallStrainDiffusionMaterial,
                     public LinearElasticFractureMaterial,
                     public NeoHookeanPFFractureMaterial,
                     public MieheFractureMaterial,
                     public SmallStrainCahnHilliardMaterial,
                     public SmallStrainDiffusionJ2Material,
                     public DiffusionACFractureMaterial{
public:
    /**
     * constructor
     */
    BulkMateSystem();


    /**
     * initialize the material libs for different material model calculation
     * @param t_matetype the type of material calculation
     * @param t_params the parameters read from json file
     * @param t_elmtinfo the local element information
     * @param t_elmtsoln the local element solution
     */
    void initBulkMateLibs(const MateType &t_matetype,const nlohmann::json &t_params,
                          const LocalElmtInfo &t_elmtinfo,const LocalElmtSolution &t_elmtsoln);

    /**
     * run the material libs for different material model calculation
     * @param t_matetype the type of material calculation
     * @param t_params the parameters read from json file
     * @param t_elmtinfo the local element information
     * @param t_elmtsoln the local element solution
     */
    void runBulkMateLibs(const MateType &t_matetype,const nlohmann::json &t_params,
                         const LocalElmtInfo &t_elmtinfo,const LocalElmtSolution &t_elmtsoln);

public:
    MaterialsContainer m_materialcontainer_old;/**< the materials container of previous step */
    MaterialsContainer m_materialcontainer;/**< the material container of current step */
    
};