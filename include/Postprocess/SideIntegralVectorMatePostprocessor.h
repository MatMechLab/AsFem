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
//+++          vector material on the specific side
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "Postprocess/SideIntegralPostprocessorBase.h"

/**
 * This class implements the integration of specific vector material on specific side
 */
class SideIntegralVectorMatePostprocessor:public SideIntegralPostprocessorBase{
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
    Vector3d m_vectorvalue=0.0;/**< the vector value */    
    int m_component=0;/**< the component of the specific vector material */
    string m_vectormatename="";/**< string name of the vector material */                
};