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
//+++ Date   : 2022.05.22
//+++ Purpose: defines the user-defined shape functions in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "FE/ShapeFunUser1.h"
#include "FE/ShapeFunUser2.h"
#include "FE/ShapeFunUser3.h"
#include "FE/ShapeFunUser4.h"
#include "FE/ShapeFunUser5.h"


#include "Mesh/Nodes.h"
#include "Mesh/MeshType.h"
#include "FE/ShapeFunType.h"
#include "MathUtils/MatrixXd.h"

/**
 * This class implement the user-defined shape function calculation
 */
class ShapeFunUser:public ShapeFunUser1,
                   public ShapeFunUser2,
                   public ShapeFunUser3,
                   public ShapeFunUser4,
                   public ShapeFunUser5{
public:
    /**
     * constructor
     */
    ShapeFunUser();
    /**
     * calculate the shape function value and its global/local derivatives
     * @param t_shp_type the shape function type
     * @param xi the local coordinate \f$\xi\f$
     * @param eta the local coordinate \f$\eta\f$
     * @param zeta the local coordinate \f$\zeta\f$
     * @param t_nodes the nodal coordinates of current mesh
     * @param flag boolean, true for global derivatives, false for local derivatives
     * @param t_vals the double vector for shape function values
     * @param t_ders the Vector3 vector for shape function derivatives
     * @param jacdet the determinte of jacobian transformation
     */
    void calcUserShapeFun(const ShapeFunType &t_shp_type,
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
    int m_nodes;/**< number of nodes */

};