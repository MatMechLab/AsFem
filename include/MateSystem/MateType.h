//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.07.09
//+++ Purpose: Define the material type in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

/**
 * This enum structure list all the pre-registered materials
 */
enum class MateType{
    NULLMATE,
    LINEARELASTICMATE,
    INCREMENTSMALLSTRAINMATE,
    SAINTVENANTMATE,
    NEOHOOKEANMATE,
    MOONEYRIVLINMATE,
    PLASTIC1DMATE,
    J2PLASTICITYMATE,
    CONSTPOISSONMATE,
    NLPOISSONMATE,
    CONSTDIFFUSIONMATE,
    NLDIFFUSIONMATE,
    THERMALMATE,
    THERMALELASTICMATE,
    CONSTWAVEMATE,
    CAHNHILLIARDMATE,
    LINEARELASTICCHMATE,
    IDEALSOLUTIONFREENERGYMATE,
    DOUBLEWELLFREENERGYMATE,
    ELASTICCAHNHILLIARDMATE,
    LINEARELASTICPFFRACTUREMATE,
    MIEHEFRACTUREMATE,
    STRESSDECOMPOSITIONMATE,
    NEOHOOKEANPFFRACTUREMATE,
    KOBAYASHIMATE,
    DIFFNEOHOOKEANMATE,
    DIFFUSIONFRACTUREMATE,
    WAVEMATE,
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
