//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
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
    POISSON,
    MECHANIC,
    CAHNHILLIARD,
    MECHCAHNHILLIARD,
    DIFFUSION,
    THERMALCONDUCT,
    THERMALMECHANICS,
    DIFFUSIONMECHANICS,
    ALLANCAHN
};