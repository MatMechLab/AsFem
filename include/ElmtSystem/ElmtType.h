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
//+++ Date   : 2020.07.01
//+++ Purpose: Define the basic element type in AsFem
//+++          The elements can be classified into two groups:
//+++             1. built-in elements,i.e. mechanics,cahnhilliard...
//+++             2. user-defined-elements(uel),i.e. user1,user2...
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once


enum class ElmtType{
    NULLELMT,
    TIMEDERIVELMT,
    LAPLACEELMT,
    POISSONELMT,
    MECHANICSELMT,
    CAHNHILLIARDELMT,
    MECHCAHNHILLIARDELMT,
    DIFFUSIONELMT,
    WAVEELMT,
    THERMALCONDUCTELMT,
    THERMALMECHANICSELMT,
    DIFFUSIONMECHANICSELMT,
    ALLANCAHNELMT,
    MIEHEFRACELMT,
    KOBAYASHIELMT,
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
