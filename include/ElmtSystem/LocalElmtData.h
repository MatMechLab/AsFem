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
//+++ Date   : 2021.08.12
//+++ Purpose: Define the data binder structure for the calculation
//+++          in local element
//+++          This structure can reduce the arguments of our local
//+++          element calculation!
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <vector>

#include "MathUtils/Vector3d.h"

using std::vector;

/**
 * This structure stores the basic information
 * for each local element's calculation.
 * This is based on the quadrature point, not the nodal point !!! 
 */
struct LocalElmtInfo{
    int m_dim; /**< the dimension of current element*/ 
    int m_nodesnum;/**< the nodes number of current element*/
    int m_dofsnum;/**< the total DoFs of current element*/
    double m_t;/**< the current time */
    double m_dt;/**< the current delta t */
    Vector3d m_gpCoords0;/**< initial/reference coordinates of current gauss point */
    Vector3d m_gpCoords;/**< current coordinates of current gauss point */
};

/**
 * This structure stores the local shape function, for instance, the shape function value, and its derivatives
 * Again, they are assigned for each quadrature point!!!
 */
struct LocalShapeFun{
    double m_test;  /**< the test function value of current gauss point*/
    double m_trial; /**< the trial function value of current gauss point*/
    Vector3d m_grad_test;/**< the gradient of the test function for current gauss point*/
    Vector3d m_grad_trial;/**< the gradient of the trial function for current gauss point*/
    Vector3d m_grad_test_current;/**< the gradient of the test function based on the current configuration for current gauss point*/
    Vector3d m_grad_trial_current;/**< the gradient of the trial function based on the current configuration for current gauss point*/

};

/**
 * This structure stores the displacement 'u', the velocity 'v' and their gradient 
 * for the local element
 */
struct LocalElmtSolution{
    vector<double> m_gpU;/**< 'displacement' of current gauss point */
    vector<double> m_gpUold;/**< previous 'displacement' of current gauss point */
    vector<double> m_gpUolder;/**< pre-previous 'displacement' of current gauss point */
    vector<double> m_gpV;/**< 'velocity' of current gauss point */
    vector<double> m_gpA;/**< 'acceleration' of current gauss point */

    vector<Vector3d> m_gpGradU;/**< 'displacement' gradient of current gauss point (reference configuration) */
    vector<Vector3d> m_gpGradUold;/**< old 'displacement' gradient of current gauss point (reference configuration) */
    vector<Vector3d> m_gpGradUolder;/**< older 'displacement' gradient of current gauss point (reference configuration) */
    vector<Vector3d> m_gpGradV;/**< 'velocity' gradient of current gauss point (reference configuration) */

    vector<Vector3d> m_gpgradu;/**< 'displacement' gradient of current gauss point in current configuration */
    vector<Vector3d> m_gpgradv;/**< 'velocity' gradient of current gauss point in current configuration */
};