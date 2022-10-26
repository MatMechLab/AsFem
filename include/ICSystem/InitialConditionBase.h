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
//+++ Purpose: Define abstract class for initial conditions
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "petsc.h"

#include "nlohmann/json.hpp"

#include "MathUtils/Vector3d.h"
#include "MathUtils/VectorXd.h"

#include "Utils/MessagePrinter.h"
#include "Utils/JsonUtils.h"

using std::string;
using std::vector;


/**
 * The abstract class for different initial conditions
 */
class InitialConditionBase{
protected:
    /**
     * compute the initial condition value for each dof (of each node)
     * @param t_params parameters read from json input file
     * @param icvalue the initial condition value
     * @param dim the dimensions of current mesh
     * @param dofs the dofs number of current IC
     * @param nodecoords the current node's coordinates
     * @param localU the initial condition solution of current IC
     */
    virtual void computeInitialValue(const nlohmann::json &t_params,
                                     const double &icvalue,
                                     const int &dim,
                                     const int &dofs,
                                     const Vector3d &nodecoords,
                                     VectorXd &localU)=0;

};
