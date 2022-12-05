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
//+++ Purpose: implement the 1d shape fun and its local derivatives
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "FE/ShapeFun1DEdge2.h"
#include "FE/ShapeFun1DEdge3.h"
#include "FE/ShapeFun1DEdge4.h"

#include "Mesh/Nodes.h"
#include "Mesh/MeshType.h"

/**
 * This is the 1d lagrange shape function class, responsible for the calculation of shape function value and
 * its global/local derivatives.
 * It should be mentioned that, this class only offers the implementation instead of data storage.
 */
class ShapeFun1D:public ShapeFun1DEdge2,
                 public ShapeFun1DEdge3,
                 public ShapeFun1DEdge4{
public:
    /**
     * constructor
     */
    ShapeFun1D();
    /**
     * calculate the shape function value and its global/local derivatives
     * @param t_meshtype the mesh type
     * @param xi the local coordinate \f$\xi\f$
     * @param t_nodes the nodal coordinates of current mesh
     * @param flag boolean, true for global derivatives, false for local derivatives
     * @param t_vals the double vector for shape function values
     * @param t_ders the Vector3 vector for shape function derivatives
     * @param jacdet the determinte of jacobian transformation
     */
    void calc1DShapeFun(const MeshType &t_meshtype,const double &xi,const Nodes &t_nodes,const bool &flag,
                        vector<double> &t_vals,
                        vector<Vector3d> &t_ders,
                        double &jacdet);

private:
    double m_dxdxi;/**< \f$\frac{dx}{d\xi}\f$ */
    double m_dydxi;/**< \f$\frac{dy}{d\xi}\f$ */
    double m_dzdxi;/**< \f$\frac{dz}{d\xi}\f$ */
    double val;/**< temp variable */
    int m_nodes;/**< number of nodes */
};