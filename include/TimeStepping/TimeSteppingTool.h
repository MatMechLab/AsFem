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
//+++ Date   : 2022.08.25
//+++ Purpose: Implement aux functions for time stepping
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "SolutionSystem/SolutionSystem.h"
#include "TimeStepping/TimeSteppingType.h"
#include "FEProblem/FEControlInfo.h"

/**
 * compute the different orders of time derivatives according to choosen time stepping method
 * @param fectrlinfo fe control information
 * @param U the intermediate PETSc Vec for current solution during iteration
 * @param soln the solution system
 */
void computeTimeDerivatives(FEControlInfo &fectrlinfo,const Vec &U,SolutionSystem &soln);
/**
 * compute the different orders of time derivatives according to choosen time stepping method
 * @param fectrlinfo fe control information
 * @param soln the solution system
 */
void computeTimeDerivatives(FEControlInfo &fectrlinfo,SolutionSystem &soln);