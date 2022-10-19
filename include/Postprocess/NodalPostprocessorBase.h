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
//+++ Date   : 2022.09.28
//+++ Purpose: Base class for the calculation of the specific nodal value
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "DofHandler/DofHandler.h"
#include "ProjectionSystem/ProjectionSystem.h"
#include "SolutionSystem/SolutionSystem.h"

#include "Utils/MessagePrinter.h"

/**
 * This class defines the abstract class for nodal type postprocess
 */
class NodalPostprocessorBase{
protected:
    /**
     * compute the nodal value for nodal pps
     * @param dofid the local dof id, start from 1
     * @param t_parameters the parameters from json
     * @param t_dofhandler the dofHandler class
     * @param t_soln the solution class
     * @param t_projsystem the projection class
     */
    virtual double computeNodalValue(const int &dofid,
                                     const nlohmann::json &t_parameters,
                                     const DofHandler &t_dofhandler,
                                     SolutionSystem &t_soln,
                                     ProjectionSystem &t_projsystem)=0;
};