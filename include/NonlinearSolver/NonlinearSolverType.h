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
//+++ Date   : 2020.07.12
//+++ Purpose: define the nonlinear solver type in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

/**
 * the enum class for different types of nonlinear solver
 */
enum class NonlinearSolverType{
    NEWTON,
    NEWTONLS,
    NEWTONSECANT,
    NEWTONTR,
    NEWTONCG,
    NEWTONGMRES,
    RICHARDSON,
    NASM,
    ASPIN,
    NMS,
    FAS,
    BFGS,
    BROYDEN,
    BADBROYDEN,
    // for user-defined nonlinear solver
    USER1,
    USER2,
    USER3,
    USER4,
    USER5
};