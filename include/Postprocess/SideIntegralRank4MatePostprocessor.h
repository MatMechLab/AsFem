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
//+++ Date   : 2022.09.30
//+++ Purpose: This class calculates the integration of the specific
//+++          rank-4 material on the specific side
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "Postprocess/SideIntegralPostprocessorBase.h"

/**
 * This class implements the integration of specific rank-4 material on specific side
 */
class SideIntegralRank4MatePostprocessor:public SideIntegralPostprocessorBase{
protected:
    /**
     * compute the nodal value for nodal pps
     * @param dofid the global dof id, start from 1
     * @param nodeid the global node id, starts from 1
     * @param t_parameters the parameters from json
     * @param t_elmtinfo the local element info
     * @param t_shp the local shape function
     * @param t_soln the solution class
     * @param t_projsystem the projection class
     */
    virtual double computeSideIntegralValue(const int &dofid,
                                            const int &nodeid,
                                            const nlohmann::json &t_parameters,
                                            const LocalElmtInfo &t_elmtinfo,
                                            const LocalShapeFun &t_shp,
                                            SolutionSystem &t_soln,
                                            ProjectionSystem &t_projsystem) override;

private:
    double m_ppsvalue=0.0;/** the postprocess result */  
    int m_i=1;/**< the i index of the specific rank4 tensor */
    int m_j=1;/**< the j index of the specific rank4 tensor */
    int m_k=1;/**< the k index of the specific rank4 tensor */
    int m_l=1;/**< the l index of the specific rank4 tensor */
    string m_rank4matename="";/**< the string name of the rank-4 material */
    Rank4Tensor m_rank4value;            
};