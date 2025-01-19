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
//+++ Date   : 2022.05.13
//+++ Purpose: defines the basic calculations of FEM
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

/**
 * the basic fem calculation type
 */
enum class FECalcType{
    COMPUTERESIDUAL,
    COMPUTEJACOBIAN,
    COMPUTERESIDUALANDJACOBIAN,
    INITMATERIAL,
    UPDATEMATERIAL,
    UPDATEU
};