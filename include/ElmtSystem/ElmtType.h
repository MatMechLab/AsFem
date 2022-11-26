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
//+++ Purpose: the element type defined in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

/**
 * This enum class defines the supported element type in AsFem
 */
enum class ElmtType{
    NULLELMT,
    // For 'single'/'uncoupled' standard element
    LAPLACEELMT,
    SCALARBODYSOURCEELMT,
    POISSONELMT,
    DIFFUSIONELMT,
    THERMALELMT,
    MECHANICSELMT,
    DYNAMICMECHANICSELMT,
    ALLENCAHNELMT,
    CAHNHILLIARDELMT,
    WAVEELMT,
    KOBAYASHIELMT,
    // For coupled elements
    STRESSCAHNHILLIARDELMT,
    STRESSDIFFUSIONELMT,
    ALLENCAHNFRACTUREELMT,
    MIEHEFRACTUREELMT,
    DIFFUSIONACFRACTUREELMT,
    // For user-defined-elements(UEL)
    USER1ELMT,
    USER2ELMT,
    USER3ELMT,
    USER4ELMT,
    USER5ELMT,
    USER6ELMT,
    USER7ELMT,
    USER8ELMT,
    USER9ELMT,
    USER10ELMT,
    //
    USER11ELMT,
    USER12ELMT,
    USER13ELMT,
    USER14ELMT,
    USER15ELMT,
    USER16ELMT,
    USER17ELMT,
    USER18ELMT,
    USER19ELMT,
    USER20ELMT
};