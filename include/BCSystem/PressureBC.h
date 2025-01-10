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
//+++ Date   : 2022.08.25
//+++ Purpose: implement the pressure type boundary condition for
//+++          solid mechanics problem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "BCSystem/IntegrateBCBase.h"


/**
 * This class implement the pressure type bc, the final contribution of this class is:
 * \f\int_{\partial\Omega}p\vec{n}N^{I}dS=\int_{\partial\Omega}\vec{t} N^{I}dS\f$
 */
class PressureBC:public IntegrateBCBase{
public:
    /**
     * execute the boundary integration calculation for different integrated boundary conditions.
     * \f$\int_{\partial\Omega}f(u,\nabla u,\vec{n})dS\f$ is the final output
     * @param CalcType the calculation type for either residual or jacbobian
     * @param BCValue the boundary value defined in the input file
     * @param Params the json content for boundary condition parameters
     * @param ElmtInfo the basic information for current element
     * @param ElmtSoln the solution of current element
     * @param Normal the normal vector of current gauss point in the current element
     * @param Shp the local shape function values
     * @param Ctan the time derivative related coefficient array
     * @param LocalK the local K matrix
     * @param LocalR the locak R matrix
     */
    virtual void computeBCValue(const FECalcType &CalcType,
                                const double &BCValue,
                                const nlohmann::json Params,
                                const LocalElmtInfo &ElmtInfo,
                                const LocalElmtSolution &ElmtSoln,
                                const Vector3d &Normal,
                                const LocalShapeFun &Shp,
                                const double (&Ctan)[3],
                                MatrixXd &LocalK,
                                VectorXd &LocalR) override;
private:
    /**
     * calculate the residual of current boundary element
     * \f$\int_{\partial\Omega}R_{u}dS\f$ is the final output
     * @param BCValue the boundary value defined in the input file
     * @param Params the boundary condition related parameters (json content)
     * @param ElmtInfo the basic information for current element
     * @param ElmtSoln the solution of current element
     * @param Normal the normal vector of current gauss point in the current element
     * @param Shp the local shape function values
     * @param LocalR the locak R matrix
     */
    virtual void computeResidual(const double &BCValue,
                                 const nlohmann::json Params,
                                 const LocalElmtInfo &ElmtInfo,
                                 const LocalElmtSolution &ElmtSoln,
                                 const Vector3d &Normal,
                                 const LocalShapeFun &Shp,
                                 VectorXd &LocalR) override; 
    
    /**
     * calculate the jacbobian contribution of current boundary element 
     * \f$\int_{\partial\Omega}\frac{\partial R{u}}{\partial u}dS\f$ is the final output
     * @param BCValue the boundary value defined in the input file
     * @param Params the boundary condition related parameters (json content)
     * @param ElmtInfo the basic information for current element
     * @param ElmtSoln the solution of current element
     * @param Normal the normal vector of current gauss point in the current element
     * @param Shp the local shape function values
     * @param LocalK the local K matrix
     */
    virtual void computeJacobian(const double &BCValue,
                                 const nlohmann::json &Params,
                                 const LocalElmtInfo &ElmtInfo,
                                 const LocalElmtSolution &ElmtSoln,
                                 const Vector3d &Normal,
                                 const LocalShapeFun &Shp,
                                 const double (&Ctan)[3],
                                 MatrixXd &LocalK) override;
private:
    Vector3d m_traction;/**< the traction vector based given pressure */
    int m_component;/**< the displacement component index */
    double m_pressure;/**< the given pressure scalar */
};