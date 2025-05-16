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
    int m_Dim; /**< the dimension of current element*/ 
    int m_NodesNum;/**< the nodes number of current element*/
    int m_DofsNum;/**< the total DoFs of current element*/
    double m_T;/**< the current time */
    double m_Dt;/**< the current delta t */
    Vector3d m_QpCoords0;/**< initial/reference coordinates of current gauss point */
    Vector3d m_QpCoords;/**< current coordinates of current gauss point */
    int m_ElmtID;/**< current element id (global one) */
    int m_ElmtsNum;/**< the total elements num */
    int m_QpointsNum;/**< the total qpoints num of current element */
    int m_QpointID;/**< current qpoint id (local one) */
    double m_TempVal;/**< temp variables to store whatever required by users */
};

/**
 * This structure stores the local shape function, for instance, the shape function value, and its derivatives
 * Again, they are assigned for each quadrature point!!!
 */
struct LocalShapeFun{
    double m_Test;  /**< the test function value of current gauss point*/
    double m_Trial; /**< the trial function value of current gauss point*/
    Vector3d m_GradTest;/**< the gradient of the test function for current gauss point*/
    Vector3d m_GradTrial;/**< the gradient of the trial function for current gauss point*/
    Vector3d m_GradTestCurrent;/**< the gradient of the test function based on the current configuration for current gauss point*/
    Vector3d m_GradTrialCurrent;/**< the gradient of the trial function based on the current configuration for current gauss point*/

};

/**
 * This structure stores the displacement 'u', the velocity 'v' and their gradient 
 * for the local element
 */
struct LocalElmtSolution{
    vector<double> m_QpU;/**< 'displacement' of current gauss point */
    vector<double> m_QpUold;/**< previous 'displacement' of current gauss point */
    vector<double> m_QpUolder;/**< pre-previous 'displacement' of current gauss point */
    vector<double> m_QpV;/**< 'velocity' of current gauss point */
    vector<double> m_QpA;/**< 'acceleration' of current gauss point */

    vector<Vector3d> m_QpGradU;/**< 'displacement' gradient of current gauss point (reference configuration) */
    vector<Vector3d> m_QpGradUold;/**< old 'displacement' gradient of current gauss point (reference configuration) */
    vector<Vector3d> m_QpGradUolder;/**< older 'displacement' gradient of current gauss point (reference configuration) */
    vector<Vector3d> m_QpGradV;/**< 'velocity' gradient of current gauss point (reference configuration) */

    vector<Vector3d> m_Qpgradu;/**< 'displacement' gradient of current gauss point in current configuration */
    vector<Vector3d> m_Qpgradv;/**< 'velocity' gradient of current gauss point in current configuration */
};