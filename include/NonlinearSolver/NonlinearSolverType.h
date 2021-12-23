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
//+++ Date   : 2020.07.12
//+++ Purpose: define the nonlinear solver type in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

enum class NonlinearSolverType{
    NEWTON,
    NEWTONLS,
    NEWTONSECANT,
    NEWTONTR,
    NEWTONCG,
    NEWTONGMRES,
    BFGS,
    BROYDEN,
    BADBROYDEN
};
