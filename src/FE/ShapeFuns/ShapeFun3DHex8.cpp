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
//+++ Purpose: 3d hex8 shape function
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/ShapeFun3DHex8.h"

ShapeFun3DHex8::ShapeFun3DHex8(){

}

void ShapeFun3DHex8::calc3DShapeValsAndDerivatives(const double &xi,const double &eta,const double &zeta,
                                                   vector<double> &t_shpvals,
                                                   vector<Vector3d> &t_shpders){
    if(t_shpvals.size()<8 || t_shpders.size()<8){
        MessagePrinter::printErrorTxt("your shape val or derivs vector size is smaller than 8, error detected in ShapeFun3DHex8.cpp");
        MessagePrinter::exitAsFem();
    }
    /**
     * The nodes should look like:
     * bottom layer:
     * 4----3
     * |    |
     * |    |
     * 1----2
     * top layer:
     * 8----7
     * |    |
     * |    |
     * 5----6
    */
    t_shpvals[0] = (1 - xi) * (1 - eta) * (1 - zeta) / 8.0;
    t_shpders[0](1) = -(1 - eta) * (1 - zeta) / 8.0;
    t_shpders[0](2) = -(1 - xi) * (1 - zeta) / 8.0;
    t_shpders[0](3) = -(1 - xi) * (1 - eta) / 8.0;

    t_shpvals[1] = (1 + xi) * (1 - eta) * (1 - zeta) / 8.0;
    t_shpders[1](1) = (1 - eta) * (1 - zeta) / 8.0;
    t_shpders[1](2) = -(1 + xi) * (1 - zeta) / 8.0;
    t_shpders[1](3) = -(1 + xi) * (1 - eta) / 8.0;

    t_shpvals[2] = (1 + xi) * (1 + eta) * (1 - zeta) / 8.0;
    t_shpders[2](1) = (1 + eta) * (1 - zeta) / 8.0;
    t_shpders[2](2) = (1 + xi) * (1 - zeta) / 8.0;
    t_shpders[2](3) = -(1 + xi) * (1 + eta) / 8.0;

    t_shpvals[3] = (1 - xi) * (1 + eta) * (1 - zeta) / 8.0;
    t_shpders[3](1) = -(1 + eta) * (1 - zeta) / 8.0;
    t_shpders[3](2) = (1 - xi) * (1 - zeta) / 8.0;
    t_shpders[3](3) = -(1 - xi) * (1 + eta) / 8.0;

    t_shpvals[4] = (1 - xi) * (1 - eta) * (1 + zeta) / 8.0;
    t_shpders[4](1) = -(1 - eta) * (1 + zeta) / 8.0;
    t_shpders[4](2) = -(1 - xi) * (1 + zeta) / 8.0;
    t_shpders[4](3) = (1 - xi) * (1 - eta) / 8.0;

    t_shpvals[5] = (1 + xi) * (1 - eta) * (1 + zeta) / 8.0;
    t_shpders[5](1) = (1 - eta) * (1 + zeta) / 8.0;
    t_shpders[5](2) = -(1 + xi) * (1 + zeta) / 8.0;
    t_shpders[5](3) = (1 + xi) * (1 - eta) / 8.0;

    t_shpvals[6] = (1 + xi) * (1 + eta) * (1 + zeta) / 8.0;
    t_shpders[6](1) = (1 + eta) * (1 + zeta) / 8.0;
    t_shpders[6](2) = (1 + xi) * (1 + zeta) / 8.0;
    t_shpders[6](3) = (1 + xi) * (1 + eta) / 8.0;

    t_shpvals[7] = (1 - xi) * (1 + eta) * (1 + zeta) / 8.0;
    t_shpders[7](1) = -(1 + eta) * (1 + zeta) / 8.0;
    t_shpders[7](2) = (1 - xi) * (1 + zeta) / 8.0;
    t_shpders[7](3) = (1 - xi) * (1 + eta) / 8.0;

}