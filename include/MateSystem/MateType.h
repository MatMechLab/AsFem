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
//+++ Date   : 2022.05.11
//+++ Purpose: the material type defined in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

/**
 * this enum class defines the supported material type in AsFem
 */
enum class MateType{
    NULLMATE,
    // For material properties used in 'single'/'uncoupled' standard element
    // for poisson material
    CONSTPOISSONMATE,
    POISSON1DBENCHMARKMATE,
    POISSON2DBENCHMARKMATE,
    NONLINEARPOISSON2DMATE,
    NONLINEARPOISSON3DMATE,
    // for diffusion material
    CONSTDIFFUSIONNMATE,
    NONLINEARDIFFUSION2DMATE,
    // for thermal material
    CONSTTHERMALMATE,
    // for mechanics material
    LINEARELASTICMATE,
    SAINTVENANTMATE,
    NEOHOOKEANMATE,
    SMALLSTRAINJ2PLASTICITYMATE,
    SMALLSTRAINEXPLAWJ2PLASTICITYMATE,
    // for free energy material
    DOUBLEWELLMATE,
    BINARYMIXMATE,
    KOBAYASHIMATE,
    // for wave material
    WAVEMATE,
    // For coupled material
    SMALLSTRAINCAHNHILLIARDMATE,
    SMALLSTRAINDIFFUSIONMATE,
    LINEARELASTICFRACMATE,
    NEOHOOKEANPFFRACTUREMATE,
    MIEHEFRACTUREMATE,
    SMALLSTRAINDIFFUSIONJ2MATE,
    DIFFUSIONACFRACTUREMATE,
    // For user-defined-materials(UMAT)
    USER1MATE,
    USER2MATE,
    USER3MATE,
    USER4MATE,
    USER5MATE,
    USER6MATE,
    USER7MATE,
    USER8MATE,
    USER9MATE,
    USER10MATE,
    //
    USER11MATE,
    USER12MATE,
    USER13MATE,
    USER14MATE,
    USER15MATE,
    USER16MATE,
    USER17MATE,
    USER18MATE,
    USER19MATE,
    USER20MATE
};