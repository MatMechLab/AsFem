//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.08.11
//+++ Purpose: implement the calculation for 3D Lagrange shape functions 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>

#include "Utils/Vector3d.h"
#include "Utils/MatrixXd.h"
#include "Utils/MessagePrinter.h"

#include "Mesh/MeshType.h"
#include "Mesh/Nodes.h"

using namespace std;

/**
 * This class implement the calculation of 3D lagrange shape functions
 */
class Lagrange3DShapeFun{
public:
    /**
     * Constructor, initiliaze variables
     */
    Lagrange3DShapeFun();
    /**
     * This function calculate the shape function and its derivatives based on the input nodes
     * @param meshtype the meshtype of 1d mesh, i.e., edge2, edge3, ...
     * @param xi the local coordinate \f$\xi\f$, not the global one (\f$X\f$)!!!
     * @param eta the local coordinate \f$\eta\f$, not the global one (\f$Y\f$)!!!
     * @param zeta the local coordinate \f$\zeta\f$, not the global one (\f$Z\f$)!!!
     * @param nodes the current element's nodes, which contains the \f$(x,y,z)\f$ coordinate of each node
     * @param shape_val the double variable vector, which stores the shape function value of each node
     * @param shape_grad the vector of vector3d variable, which stores the derivatives with respect to the local/global derivatives for each node
     * @param detjac the determinte of the jacobian matrix for current element 
     * @param flag if flag=true(default), then we calculate the derivatives with respecte to global coordinate \f$(x,y,z)\f$. Otherwise (flag=false), we calculte the derivatives by using the local coordinates \f$(\xi)\f$
     */ 
    void Calc3DShapeFun(const MeshType &meshtype,const double &xi,const double &eta,const double &zeta,const Nodes &nodes,vector<double> &shape_val,vector<Vector3d> &shape_grad,double &detjac,bool flag=true);

private:

    MatrixXd _Jac,_XJac;/**< temporary matrix for different calculation*/

    double _dxdxi;/**< the local derivative of x over \f$\xi\f$, dxdxi=\f$\frac{\partial x}{\partial\xi}\f$*/
    double _dxdeta;/**< the local derivative of x over \f$\eta\f$, dxdeta=\f$\frac{\partial x}{\partial\eta}\f$*/
    double _dxdzeta;/**< the local derivative of x over \f$\zeta\f$, dxdzeta=\f$\frac{\partial x}{\partial\zeta}\f$*/

    double _dydxi;/**< the local derivative of y over \f$\xi\f$, dydxi=\f$\frac{\partial y}{\partial\xi}\f$*/
    double _dydeta;/**< the local derivative of y over \f$\eta\f$, dydeta=\f$\frac{\partial y}{\partial\eta}\f$*/
    double _dydzeta;/**< the local derivative of y over \f$\zeta\f$, dydzeta=\f$\frac{\partial y}{\partial\zeta}\f$*/

    double _dzdxi;/**< the local derivative of z over \f$\xi\f$, dzdxi=\f$\frac{\partial z}{\partial\xi}\f$*/
    double _dzdeta;/**< the local derivative of z over \f$\eta\f$, dzdeta=\f$\frac{\partial z}{\partial\eta}\f$*/
    double _dzdzeta;/**< the local derivative of z over \f$\zeta\f$, dzdzeta=\f$\frac{\partial z}{\partial\zeta}\f$*/

    int _nNodes;/**< number of node number for current shape function calculation */

    double _tol=1.0e-15; /**< tolerance for error check*/
    double _val;/**< temporary variable*/


};
