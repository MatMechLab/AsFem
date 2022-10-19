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
//+++ Date   : 2022.05.14
//+++ Purpose: 2d shape function
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "FE/ShapeFun2DTri3.h"
#include "FE/ShapeFun2DTri6.h"
#include "FE/ShapeFun2DQuad4.h"
#include "FE/ShapeFun2DQuad8.h"
#include "FE/ShapeFun2DQuad9.h"


#include "Mesh/Nodes.h"
#include "Mesh/MeshType.h"

/**
 * This is the 2d lagrange shape function class, responsible for the calculation of shape function value and
 * its global/local derivatives.
 * It should be mentioned that, this class only offers the implementation instead of data storage.
 */
class ShapeFun2D:public ShapeFun2DTri3,
                 public ShapeFun2DTri6,
                 public ShapeFun2DQuad4,
                 public ShapeFun2DQuad8,
                 public ShapeFun2DQuad9{
public:
    /**
     * constructor
     */
    ShapeFun2D();
    /**
     * calculate the shape function value and its global/local derivatives
     * @param t_meshtype the mesh type
     * @param xi the local coordinate \f$\xi\f$
     * @param eta the local coordinate \f$\eta\f$
     * @param t_nodes the nodal coordinates of current mesh
     * @param flag boolean, true for global derivatives, false for local derivatives
     * @param t_vals the double vector for shape function values
     * @param t_ders the Vector3 vector for shape function derivatives
     * @param jacdet the determinte of jacobian transformation
     */
    void calc2DShapeFun(const MeshType &t_meshtype,const double &xi,const double &eta,
                        const Nodes &t_nodes,const bool &flag,
                        vector<double> &t_vals,
                        vector<Vector3d> &t_ders,
                        double &jacdet);

private:
    double m_dxdxi;/**< \f$\frac{dx}{d\xi}\f$ */
    double m_dydxi;/**< \f$\frac{dy}{d\xi}\f$ */
    double m_dzdxi;/**< \f$\frac{dz}{d\xi}\f$ */

    double m_dxdeta;/**< \f$\frac{dx}{d\eta}\f$ */
    double m_dydeta;/**< \f$\frac{dy}{d\eta}\f$ */
    double m_dzdeta;/**< \f$\frac{dz}{d\eta}\f$ */

    double m_dxidx;/**< \f$\frac{d\xi}{dx}\f$*/
    double m_dxidy;/**< \f$\frac{d\xi}{dy}\f$*/
    double m_dxidz;/**< \f$\frac{d\xi}{dz}\f$*/

    double m_detadx;/**< \f$\frac{d\eta}{dx}\f$*/
    double m_detady;/**< \f$\frac{d\eta}{dy}\f$*/
    double m_detadz;/**< \f$\frac{d\eta}{dz}\f$*/

    double valx,valy;/**< temp variable */
    double jac11,jac12,jac21,jac22;/**< variables for jacobian matrix */
    double xjac11,xjac12,xjac21,xjac22;/**< variables for inverse jacobian matrix */
    int m_nodes;/**< number of nodes */

};