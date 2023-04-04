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
//+++ Purpose: 2d tri6 shape function
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/ShapeFun2DTri6.h"

ShapeFun2DTri6::ShapeFun2DTri6(){}

void ShapeFun2DTri6::calc2DShapeValsAndDerivatives(const double &xi,const double &eta,vector<double> &t_shpvals,vector<Vector3d> &t_shpders){
    if(t_shpvals.size()<6 || t_shpders.size()<6){
        MessagePrinter::printErrorTxt("your shape val or derivs vector size is smaller than 6, error detected in Shape2DTri6.cpp");
        MessagePrinter::exitAsFem();
    }
    // taken from: http://www.sd.ruhr-uni-bochum.de/downloads/Shape_funct.pdf
    /**
     * VTK Cell type: vtkQuadraticTriangle
     * The nodes should look like:
     * 3
     * |\
     * | \
     * 6  5
     * |   \
     * |    \
     * 1--4--2
    */
    t_shpvals[1-1]=(1.0-xi-eta)*(1.0-2*xi-2*eta);
    t_shpders[1-1](1)=-3.0+4.0*eta+4.0*xi;
    t_shpders[1-1](2)=-3.0+4.0*eta+4.0*xi;
    t_shpders[1-1](3)= 0.0;

    t_shpvals[2-1]= xi*(2.0*xi-1.0);
    t_shpders[2-1](1)= 4.0*xi-1.0;
    t_shpders[2-1](2)= 0.0;
    t_shpders[2-1](3)= 0.0;

    t_shpvals[3-1]= eta*(2.0*eta-1.0);
    t_shpders[3-1](1)= 0.0;
    t_shpders[3-1](2)= 4.0*eta-1.0;
    t_shpders[3-1](3)= 0.0;

    t_shpvals[4-1]= 4.0*xi*(1.0-xi-eta);
    t_shpders[4-1](1)= 4.0*(1.0-eta-2.0*xi);
    t_shpders[4-1](2)=-4.0*xi;
    t_shpders[4-1](3)= 0.0;

    t_shpvals[5-1]= 4.0*xi*eta;
    t_shpders[5-1](1)= 4.0*eta;
    t_shpders[5-1](2)= 4.0*xi;
    t_shpders[5-1](3)= 0.0;

    t_shpvals[6-1]= 4.0*eta*(1.0-xi-eta);
    t_shpders[6-1](1)=-4.0*eta;
    t_shpders[6-1](2)= 4.0*(1.0-2.0*eta-xi);
    t_shpders[6-1](3)= 0.0;

}
