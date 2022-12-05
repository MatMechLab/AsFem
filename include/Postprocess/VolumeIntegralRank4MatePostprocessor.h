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
//+++ Date   : 2022.10.10
//+++ Purpose: This class calculates the volume integral for the
//+++          specific rank-4 material on specific domain
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "Postprocess/VolumeIntegralPostprocessorBase.h"

/**
 * This class implements the volume integral calculation of specific rank-4 material on the
 * specific domain
 */
class VolumeIntegralRank4MatePostprocessor:public VolumeIntegralPostprocessorBase{
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
                                              ProjectionSystem &t_projsystem) override;

private:
    int m_i=1;/**< the i index of the specific rank4 tensor */
    int m_j=1;/**< the j index of the specific rank4 tensor */
    int m_k=1;/**< the k index of the specific rank4 tensor */
    int m_l=1;/**< the l index of the specific rank4 tensor */
    string m_rank4matename="";/**< the string name of the scalar material */
    Rank4Tensor m_rank4value;/**< the rank-4 tensor material value */
    double m_ppsvalue=0.0;/** the postprocess result */

};