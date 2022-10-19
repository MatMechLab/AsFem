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
//+++ Purpose: Get the nodal rank4 material value via its node id for pps
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "Postprocess/NodalPostprocessorBase.h"

/**
 * This class access the nodal's rank4 material value via its global node id
 */
class NodalRank4MatePostprocessor:public NodalPostprocessorBase{
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
    int m_i=1;/**< the i index of the specific rank4 tensor */
    int m_j=1;/**< the j index of the specific rank4 tensor */
    int m_k=1;/**< the k index of the specific rank4 tensor */
    int m_l=1;/**< the l index of the specific rank4 tensor */
    string m_rank4matename="";/**< the string name of the rank-4 material */
    Rank4Tensor m_pps_value=0.0;/**< the postprocess result */

};