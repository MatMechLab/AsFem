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
//+++ Purpose: 3d shape function
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "FE/ShapeFun3DTet4.h"
#include "FE/ShapeFun3DTet10.h"
#include "FE/ShapeFun3DHex8.h"
#include "FE/ShapeFun3DHex20.h"
#include "FE/ShapeFun3DHex27.h"


#include "Mesh/Nodes.h"
#include "Mesh/MeshType.h"
#include "MathUtils/MatrixXd.h"

/**
 * This is the 3d lagrange shape function class, responsible for the calculation of shape function value and
 * its global/local derivatives.
 * It should be mentioned that, this class only offers the implementation instead of data storage.
 */
class ShapeFun3D:public ShapeFun3DTet4,
                 public ShapeFun3DTet10,
                 public ShapeFun3DHex8,
                 public ShapeFun3DHex20,
                 public ShapeFun3DHex27{
public:
    /**
     * constructor
     */
    ShapeFun3D();
    /**
     * calculate the shape function value and its global/local derivatives
     * @param t_meshtype the mesh type
     * @param xi the local coordinate \f$\xi\f$
     * @param eta the local coordinate \f$\eta\f$
     * @param zeta the local coordinate \f$\zeta\f$
     * @param t_nodes the nodal coordinates of current mesh
     * @param flag boolean, true for global derivatives, false for local derivatives
     * @param t_vals the double vector for shape function values
     * @param t_ders the Vector3 vector for shape function derivatives
     * @param jacdet the determinte of jacobian transformation
     */
    void calc3DShapeFun(const MeshType &t_meshtype,
                        const double &xi,const double &eta,const double &zeta,
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

    double m_dxdzeta;/**< \f$\frac{dx}{d\zeta}\f$ */
    double m_dydzeta;/**< \f$\frac{dy}{d\zeta}\f$ */
    double m_dzdzeta;/**< \f$\frac{dz}{d\zeta}\f$ */

    MatrixXd m_jac;/**< jacobian transformation matrix */
    MatrixXd m_xjac;/**< the inverse of jacobian matrix */
    double temp1,temp2,temp3;/**< temp variables */
    int m_nodes;/**< number of nodes */

};