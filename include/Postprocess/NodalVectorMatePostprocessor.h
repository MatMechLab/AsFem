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
//+++ Purpose: Get the nodal vector material value via its node id for pps
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "Postprocess/NodalPostprocessorBase.h"

/**
 * This class access the nodal's vector material value via its global node id
 */
class NodalVectorMatePostprocessor:public NodalPostprocessorBase{
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
                                     ProjectionSystem &t_projsystem) override;

private:
    int m_nodeid=0;/**< the global node id */
    int m_component=1;/**< the component of the specific vector */
    string m_vectormatename="";/**< the string name of the scalar material */
    Vector3d m_pps_value=0.0;/**< the postprocess result */

};