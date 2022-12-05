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
//+++ Date   : 2022.08.10
//+++ Purpose: Defines the abstract class for free energy materials
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "MathUtils/Vector3d.h"
#include "MathUtils/VectorXd.h"
#include "MathUtils/MatrixXd.h"

/**
 * This class defines the necessary functions for free energy materials
*/
class FreeEnergyMaterialBase{
protected:
    /**
     * Calculate the free energy and its first/second order (partial) derivatives.
     * @param t_parameters the json parameters read from input file
     * @param t_args the variables for the free energy expression, it could be concentration or order parameters
     * @param t_F the system free energy
     * @param t_dFdargs the first order derivatives of F with respect to its own args
     * @param t_d2Fdargs2 the second order (partial) derivatives of F with respect to different args, off-diagnoal part for partial derivatives
    */
    virtual void computeFreeEnergyAndDerivatives(const nlohmann::json &t_parameters,
                                                 const VectorXd &t_args,
                                                 VectorXd       &t_F,
                                                 VectorXd       &t_dFdargs,
                                                 MatrixXd       &t_d2Fdargs2)=0;
};