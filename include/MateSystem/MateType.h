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
    IDEALSOLUTIONFREENERGYMATE,
    DOUBLEWELLFREENERGYMATE,
    ELASTICCAHNHILLIARDMATE,
    LINEARELASTICPFFRACTUREMATE,
    MIEHEFRACTUREMATE,
    HYPERELASTICPFFRACTUREMATE,
    DENDRITEMATE,
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
