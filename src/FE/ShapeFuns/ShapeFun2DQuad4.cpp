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
//+++ Purpose: 2d quad4 shape function
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/ShapeFun2DQuad4.h"

ShapeFun2DQuad4::ShapeFun2DQuad4(){}

void ShapeFun2DQuad4::calc2DShapeValsAndDerivatives(const double &xi,const double &eta,vector<double> &t_shpvals,vector<Vector3d> &t_shpders){
    if(t_shpvals.size()<4 || t_shpders.size()<4){
        MessagePrinter::printErrorTxt("your shape val or derivs vector size is smaller than 4, error detected in ShapeFun2DQuad4.cpp");
        MessagePrinter::exitAsFem();
    }
    /**
     * VTK Cell type: vtkQuad
     * The nodes should look like:
     * 4----3
     * |    |
     * |    |
     * 1----2
    */
    t_shpvals[1-1]=(1.0-xi)*(1.0-eta)/4.0;
    t_shpders[1-1](1)=(eta-1.0)/4.0;
    t_shpders[1-1](2)=(xi -1.0)/4.0;
    t_shpders[1-1](3)= 0.0;

    t_shpvals[2-1]=(1.0+xi)*(1.0-eta)/4.0;
    t_shpders[2-1](1)= (1.0-eta)/4.0;
    t_shpders[2-1](2)=-(1.0+xi )/4.0;
    t_shpders[2-1](3)= 0.0;

    t_shpvals[3-1]=(1.0+xi)*(1.0+eta)/4.0;
    t_shpders[3-1](1)= (1.0+eta)/4.0;
    t_shpders[3-1](2)= (1.0+xi )/4.0;
    t_shpders[3-1](3)= 0.0;

    t_shpvals[4-1]=(1.0-xi)*(1.0+eta)/4.0;
    t_shpders[4-1](1)=-(1.0+eta)/4.0;
    t_shpders[4-1](2)= (1.0-xi )/4.0;
    t_shpders[4-1](3)= 0.0;

}
