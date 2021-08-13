//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
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
#include "Utils/Vector3d.h"

using namespace std;

/**
 * This structure store the basic information
 * for each local element's calculation.
 * This is based on the quadrature point, not nodal point !!! 
 */
struct LocalElmtInfo{
    int nDim; /**< the dimension of current element*/ 
    int nNodes;/**< the nodes number of current element*/
    int nDofs;/**< the total DoFs of current element*/
    double t;/**< the current time */
    double dt;/**< the current delta t */
    Vector3d gpCoords;/**< the current gauss point's coordinates (x,y,z)*/
};

/**
 * This structure stores the local shape function, for instance, the shape function value, and its derivatives
 * Again, they are assigned for each quadrature point!!!
 */
struct LocalShapeFun{
    double test;  /**< the test function value of current gauss point*/
    double trial; /**< the trial function value of current gauss point*/
    Vector3d grad_test;/**< the gradient of the test function for current gauss point*/
    Vector3d grad_trial;/**< the gradient of the trial function for current gauss point*/
    Vector3d grad_test_current;/**< the gradient of the test function based on the current configuration for current gauss point*/
    Vector3d grad_trial_current;/**< the gradient of the trial function based on the current configuration for current gauss point*/

};




/**
 * This structure stores the displacement 'u', the velocity 'v' and their gradient 
 * for the local element
 */
struct LocalElmtSolution{
    vector<double> gpU;   /**< the solution vector of local displacement 'u'*/
    vector<double> gpUold;/**< the solution vector of local displacement 'u' in the previous step*/
    vector<double> gpV;   /**< the solution vector of local displacement 'v'*/
    vector<double> gpVold;/**< the solution vector of local displacement 'v' in the previous step*/
    vector<Vector3d> gpGradU;   /**< the gradient vector of local disp 'u'*/
    vector<Vector3d> gpGradUold;/**< the gradient vector of local disp 'u' in the previous step*/
    vector<Vector3d> gpGradV;   /**< the gradient vector of local disp 'u'*/
    vector<Vector3d> gpGradVold;/**< the gradient vector of local disp 'u' in the previous step*/
};
