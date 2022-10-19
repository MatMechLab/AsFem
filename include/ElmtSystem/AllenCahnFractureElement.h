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
//+++ Date   : 2021.12.25
//+++ Purpose: implement the residual and jacobian for allen-cahn 
//+++          type phase-field fracture model
//+++          1) dD/dt=-L*(delta f/delta D)
//+++          2) div(\Sigma)=0
//+++ Ref    : A continuum phase field model for fracture, by Prof. M\"uller
//+++ DOI    : https://doi.org/10.1016/j.engfracmech.2010.08.009
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "ElmtSystem/BulkElmtBase.h"


/**
 * This class implement the calculation of allen-cahn type phase-field fracture equations
 */
class AllenCahnFractureElement:public BulkElmtBase{
public:
    /**
     * This function responsible for the residual, jacobian
     * @param t_calctype the calculation type of FEM analysis, i.e., residual-calc, jacobian-calc, projection-calc
     * @param t_elmtinfo the structure which contains the nodes numer, dimension, dofs num, quadrature point coordinates information
     * @param t_ctan 1x3 vector, where ctan[0] is responsible for the non-time-derivative part in the K matrix, while ctan[1] and ctan[2] represent the coeffecient for the 1st and 2nd order time derivatives in the K matrix
     * @param t_soln the solution structure, which contains the local displacement 'u' and velocity 'v' vector, as well as their derivatives
     * @param t_shp the shape function structure, which stores the test/trial function value and their derivatives
     * @param t_mate_old the materials of previous step
     * @param t_mate the materials of current step
     * @param localK the K matrix of local element, it could be 1x1 or 5x5 depends on your local nDOfs
     * @param localR the Residual vector of local element,it could 1 or 5 depends on your local nDofs
     */
    virtual void computeAll(const FECalcType &t_calctype,const LocalElmtInfo &t_elmtinfo,const double (&ctan)[3],
            const LocalElmtSolution &t_soln,const LocalShapeFun &t_shp,
            const MaterialsContainer &t_mate_old,const MaterialsContainer &t_mate,
            MatrixXd &localK,VectorXd &localR) override;

protected:
    /**
     * This function is responsible for the local residual vector calculation
     * @param elmtinfo the structure which contains the nodes numer, dimension, dofs num, quadrature point coordinates information
     * @param soln the solution structure, which contains the local displacement 'u' and velocity 'v' vector, as well as their derivatives
     * @param shp the shape function structure, which stores the test/trial function value and their derivatives
     * @param t_mate_old the materials of previous step
     * @param t_mate the materials of current step
     * @param localR the Residual vector of local element,it could 1 or 5 depends on your local nDofs 
     */
    virtual void computeResidual(const LocalElmtInfo &t_elmtinfo,
                                 const LocalElmtSolution &t_soln,
                                 const LocalShapeFun &t_shp,
                                 const MaterialsContainer &t_mate_old,const MaterialsContainer &t_mate,
                                 VectorXd &localR) override;

    /**
     * This function responsible for the  jacobian matrix calculation
     * @param t_elmtinfo the structure which contains the nodes numer, dimension, dofs num, quadrature point coordinates information
     * @param ctan 1x3 vector, where ctan[0] is responsible for the non-time-derivative part in the K matrix, while ctan[1] represents the coeffecient for the time derivatives in the K matrix
     * @param t_soln the solution structure, which contains the local displacement 'u' and velocity 'v' vector, as well as their derivatives
     * @param t_shp the shape function structure, which stores the test/trial function value and their derivatives
     * @param t_mate_old the materials of previous step
     * @param t_mate the materials of current step
     * @param localK the K matrix of local element, it could be 1x1 or 5x5 depends on your local nDOfs
     */
    virtual void computeJacobian(const LocalElmtInfo &t_elmtinfo,const double (&ctan)[3],
                                 const LocalElmtSolution &t_soln,
                                 const LocalShapeFun &t_shp,
                                 const MaterialsContainer &t_mate_old,const MaterialsContainer &t_mate,
                                 MatrixXd &localK) override;
private:

    /**
     * The degradation function
     * @param x the damage variable, 0-> undamage,1->full damage
    */
    double g(const double &x){
        return (1.0-x)*(1.0-x);
    }

    /**
     * The 1st order derivative of the degradation function
     * @param x the damage variable, 0-> undamage,1->full damage
    */
    double dg(const double &x){
        return 2.0*(x-1.0);
    }

    /**
     * The 2nd order derivative of the degradation function
     * @param x the damage variable, 0-> undamage,1->full damage
    */
    double d2g(const double &x){
        if(x){}
        return 2.0;
    }
};