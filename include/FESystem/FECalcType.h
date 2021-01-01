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
//+++ Date   : 2020.11.29
//+++ Purpose: Define some commonly used calculation type in FEM
//+++          calculation, i.e. compute residual, compute jacobian
//+++          projection from gauss point to nodal point
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

enum class FECalcType{
    ComputeResidual,
    ComputeJacobian,
    Projection,
    InitHistoryVariable,
    UpdateHistoryVariable
};