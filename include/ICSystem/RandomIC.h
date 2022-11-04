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
//+++ Purpose: Implement the constant IC, the given dofs will be
//+++          assigned with a constant value
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "ICSystem/InitialConditionBase.h"

class RandomIC:public InitialConditionBase{
// public:
    // ~RandomIC(){
    //     PetscRandomDestroy(&m_rnd);
    // }
protected:
    /**
     * compute the initial condition value for each dof (of each node)
     * @param t_params parameters read from json input file
     * @param icvalue the initial condition value
     * @param dim the dimensions of current mesh
     * @param dofs the dofs number of current IC
     * @param nodecoords the current node's coordinates
     * @param localU the initial condition values of current IC
     */
    virtual void computeInitialValue(const nlohmann::json &t_params,
                                     const double &icvalue,
                                     const int &dim,
                                     const int &dofs,
                                     const Vector3d &nodecoords,
                                     VectorXd &localU) override;

};