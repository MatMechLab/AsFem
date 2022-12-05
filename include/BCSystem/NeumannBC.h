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
//+++ Date   : 2021.10.06
//+++ Purpose: implement the Neumann type boundary condition
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "BCSystem/IntegrateBCBase.h"


/**
 * This class implement the flux type bc, the final contribution of this class is:
 * \f\int_{\partial\Omega}\nabla u\cdot\vec{n}N^{I}dS=\int_{\partial\Omega}J N^{I}dS\f$
 */
class NeumannBC:public IntegrateBCBase{
public:
    /**
     * execute the boundary integration calculation for different integrated boundary conditions.
     * \f$\int_{\partial\Omega}f(u,\nabla u,\vec{n})dS\f$ is the final output
     * @param t_calctype the calculation type for either residual or jacbobian
     * @param t_bcvalue the boundary value defined in the input file
     * @param t_json the json content for boundary condition parameters
     * @param t_elmtinfo the basic information for current element
     * @param t_soln the solution of current element
     * @param normal the normal vector of current gauss point in the current element
     * @param t_shp the local shape function values
     * @param ctan the time derivative related coefficient array
     * @param localK the local K matrix
     * @param localR the locak R matrix
     */
    virtual void computeBCValue(const FECalcType &t_calctype,const double &t_bcvalue,
                                const nlohmann::json t_json,
                                const LocalElmtInfo &t_elmtinfo,
                                const LocalElmtSolution &t_soln,
                                const Vector3d &normal,
                                const LocalShapeFun &t_shp,
                                const double (&ctan)[3],
                                MatrixXd &localK,
                                VectorXd &localR) override;
private:
    /**
     * calculate the residual of current boundary element
     * \f$\int_{\partial\Omega}R_{u}dS\f$ is the final output
     * @param t_bcvalue the boundary value defined in the input file
     * @param t_json the boundary condition related parameters (json content)
     * @param t_elmtinfo the basic information for current element
     * @param t_elmtsoln the solution of current element
     * @param normal the normal vector of current gauss point in the current element
     * @param shp the local shape function values
     * @param localR the locak R matrix
     */
    virtual void computeResidual(const double &t_bcvalue,
                                 const nlohmann::json t_json,
                                 const LocalElmtInfo &t_elmtinfo,
                                 const LocalElmtSolution &t_elmtsoln,
                                 const Vector3d &normal,
                                 const LocalShapeFun &t_shp,
                                 VectorXd &localR) override; 
    
    /**
     * calculate the jacbobian contribution of current boundary element 
     * \f$\int_{\partial\Omega}\frac{\partial R{u}}{\partial u}dS\f$ is the final output
     * @param t_bcvalue the boundary value defined in the input file
     * @param t_json the boundary condition related parameters (json content)
     * @param t_elmtinfo the basic information for current element
     * @param t_soln the solution of current element
     * @param normal the normal vector of current gauss point in the current element
     * @param t_shp the local shape function values
     * @param localK the local K matrix
     */
    virtual void computeJacobian(const double &t_bcvalue,
                                 const nlohmann::json &t_json,
                                 const LocalElmtInfo &t_elmtinfo,
                                 const LocalElmtSolution &t_soln,
                                 const Vector3d &normal,
                                 const LocalShapeFun &t_shp,
                                 const double (&ctan)[3],
                                 MatrixXd &localK) override;
};