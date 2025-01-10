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
 * For User-Defined-Material (UMAT)
*/
#include "MateSystem/User1Material.h"

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
                     public DiffusionACFractureMaterial,
                     // for umat
                     public User1Material{
public:
    /**
     * constructor
     */
    BulkMateSystem();


    /**
     * initialize the material libs for different material model calculation
     * @param t_MateType the type of material calculation
     * @param Params the parameters read from json file
     * @param ElmtInfo the local element information
     * @param ElmtSoln the local element solution
     */
    void initBulkMateLibs(const MateType &t_MateType,
                          const nlohmann::json &Params,
                          const LocalElmtInfo &ElmtInfo,
                          const LocalElmtSolution &ElmtSoln);

    /**
     * run the material libs for different material model calculation
     * @param t_MateType the type of material calculation
     * @param Params the parameters read from json file
     * @param ElmtInfo the local element information
     * @param ElmtSoln the local element solution
     */
    void runBulkMateLibs(const MateType &t_MateType,
                         const nlohmann::json &Params,
                         const LocalElmtInfo &ElmtInfo,
                         const LocalElmtSolution &ElmtSoln);

public:
    MaterialsContainer m_MaterialContainerOld;/**< the materials container of previous step */
    MaterialsContainer m_MaterialContainer;/**< the material container of current step */
    
};