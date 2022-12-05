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
//+++ Date   : 2022.08.12
//+++ Purpose: define the nonlinear solver abstrct class in AsFem
//+++          all the nonlinear solver should inherit from here
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "Utils/MessagePrinter.h"

#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "FE/FE.h"
#include "ElmtSystem/ElmtSystem.h"
#include "MateSystem/MateSystem.h"
#include "FESystem/FESystem.h"
#include "BCSystem/BCSystem.h"
#include "ProjectionSystem/ProjectionSystem.h"
#include "SolutionSystem/SolutionSystem.h"
#include "EquationSystem/EquationSystem.h"
#include "FEProblem/FEControlInfo.h"


/**
 * This class defines the abstract class for the nonlinear solver in AsFem
 */
class NonlinearSolverBase{
public:
    /**
     * solve the nonlinear equation, if success then return true
     * @param t_mesh the mesh class
     * @param t_dofhandler the dof class
     * @param t_fe the fe class
     * @param t_elmtsyste the element system class
     * @param t_matesystem the material system class
     * @param t_fesystem the fe system class
     * @param t_bcsystem the boundary condition system
     * @param t_solutionsystem the solution system class
     * @param t_equationsystem the equation system class
     * @param t_fectrlinfo the fe control info
     */
    virtual bool solve(Mesh &t_mesh,DofHandler &t_dofhandler,FE &t_fe,
                       ElmtSystem &t_elmtsyste,MateSystem &t_matesystem,
                       FESystem &t_fesystem,
                       BCSystem &t_bcsystem,
                       SolutionSystem &t_solutionsystem,
                       EquationSystem &t_equationsystem,
                       FEControlInfo &t_fectrlinfo)=0;
};