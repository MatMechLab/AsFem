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
//+++ Date   : 2022.10.10
//+++ Purpose: Base class for the calculation of volume integral type
//+++          postprocess ('volume' mesh integration)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "SolutionSystem/SolutionSystem.h"
#include "ProjectionSystem/ProjectionSystem.h"

#include "MathUtils/Vector.h"

#include "Utils/MessagePrinter.h"

/**
 * This class defines the abstract class for volume integral type postprocess
 */
class VolumeIntegralPostprocessorBase{
protected:
    /**
     * compute the integral value for volume integral pps
     * @param dofid the global dof id, start from 1
     * @param nodeid the global node id, starts from 1
     * @param t_parameters the parameters from json
     * @param t_elmtinfo the local element info structure
     * @param t_shp the local shape function
     * @param t_soln the solution class
     * @param t_projsystem the projection class
     */
    virtual double computeVolumeIntegralValue(const int &dofid,
                                              const int &nodeid,
                                              const nlohmann::json &t_parameters,
                                              const LocalElmtInfo &t_elmtinfo,
                                              const LocalShapeFun &t_shp,
                                              SolutionSystem &t_soln,
                                              ProjectionSystem &t_projsystem)=0;

};