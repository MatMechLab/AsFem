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
//+++ Purpose: implement the edge4 shape fun and its local derivatives
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "FE/ShapeFun1DBase.h"

/**
 * This class implement the edge4's shape function calculation
 */
class ShapeFun1DEdge4:public ShapeFun1DBase{
public:
    ShapeFun1DEdge4();

protected:
    /**
     * this override function responsible for the calculation of 1d shape function and its local derivatives
     * @param xi \f$\eta\f$ for local coordinate
     * @param t_shpvals the vector of shape function values
     * @param t_shpders the vector of shape function's local derivatives, \f$\frac{dN}{d\xi}\f$
     */
    virtual void calc1DShapeValsAndDerivatives(const double &xi,vector<double> &t_shpvals,vector<Vector3d> &t_shpders) override;

};