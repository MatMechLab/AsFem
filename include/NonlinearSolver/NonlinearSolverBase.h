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
//+++ Date   : 2022.08.12
//+++ Purpose: define the nonlinear solver abstrct class in AsFem
//+++          all the nonlinear solver should inherit from here
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "Utils/MessagePrinter.h"

#include "FECell/FECell.h"
#include "DofHandler/DofHandler.h"
#include "FE/FE.h"
#include "ElmtSystem/ElmtSystem.h"
#include "MateSystem/MateSystem.h"
#include "FESystem/FESystem.h"
#include "BCSystem/BCSystem.h"
#include "SolutionSystem/SolutionSystem.h"
#include "EquationSystem/EquationSystem.h"
#include "LinearSolver/LinearSolver.h"
#include "FEProblem/FEControlInfo.h"


/**
 * This class defines the abstract class for the nonlinear solver in AsFem
 */
class NonlinearSolverBase{
public:
    /**
     * solve the nonlinear equation, if success then return true
     * @param t_FECell the fe cell class
     * @param t_DofHandler the dof class
     * @param t_FE the fe class
     * @param t_ElmtSystem the element system class
     * @param t_MateSystem the material system class
     * @param t_FESystem the fe system class
     * @param t_BCSystem the boundary condition system
     * @param t_SolnSystem the solution system class
     * @param t_EqSystem the equation system class
     * @param t_LinearSolver the linear solver system
     * @param t_FECtrlInfo the fe control info
     */
    virtual bool solve(FECell &t_FECell,
                       DofHandler &t_DofHandler,
                       FE &t_FE,
                       ElmtSystem &t_ElmtSystem,
                       MateSystem &t_MateSystem,
                       FESystem &t_FESystem,
                       BCSystem &t_BCSystem,
                       SolutionSystem &t_SolnSystem,
                       EquationSystem &t_EqSystem,
                       LinearSolver &t_LinearSolver,
                       FEControlInfo &t_FECtrlInfo)=0;
};