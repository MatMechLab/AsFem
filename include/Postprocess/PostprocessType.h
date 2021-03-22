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
//+++ Date   : 2021.02.21
//+++ Purpose: define the postprocess type in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

enum class PostprocessType{
    NULLPPS,
    // For node and node set pps
    NODALVALUEPPS,
    NODESETVALUEPPS,
    NODESETMAXVALUEPPS,
    NODESETMINVALUEPPS,
    NODESETMIDVALUEPPS,
    // For sub-elemental(side or boundary) type pps
    AREAPPS,
    SIDEINTEGRALPPS,
    // For elemental pps
    VOLUMEPPS,
    ELEMENTINTEGRALPPS,
    ELEMENTVALUEPPS,
    ELEMENTSETVALUEPPS,
    ELEMENTSETMAXVALUEPPS,
    ELEMENTSETMINVALUEPPS,
    ELEMENTSETMIDVALUEPPS,
    // For reaction force
    PROJVARIABLESIDEINTEGRALPPS,
    // For rank-2 materials
    RANK2MATESIDEINTEGRALPPS
};