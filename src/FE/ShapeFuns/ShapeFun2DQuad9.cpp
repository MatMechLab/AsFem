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
//+++ Purpose: 2d quad9 shape function
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/ShapeFun2DQuad9.h"

ShapeFun2DQuad9::ShapeFun2DQuad9(){}

void ShapeFun2DQuad9::calc2DShapeValsAndDerivatives(const double &xi,const double &eta,vector<double> &t_shpvals,vector<Vector3d> &t_shpders){
    if(t_shpvals.size()<9 || t_shpders.size()<9){
        MessagePrinter::printErrorTxt("your shape val or derivs vector size is smaller than 9, error detected in ShapeFun2DQuad9.cpp");
        MessagePrinter::exitAsFem();
    }
    /**
     * VTK Cell type: vtkQuadraticQuad
     * The nodes should look like:
     * 4---7---3
     * |   |   |
     * 8---9---6
     * |   |   |
     * 1---5---2
    */
    t_shpvals[1-1]=(xi*xi-xi )*(eta*eta-eta)/4.0;
    t_shpders[1-1](1)= (2.0*xi-1.0)*(eta*eta-eta)/4.0;
    t_shpders[1-1](2)= (xi*xi-xi  )*(2.0*eta-1.0)/4.0;
    t_shpders[1-1](3)= 0.0;

    t_shpvals[2-1]=(xi*xi+xi )*(eta*eta-eta)/4.0;
    t_shpders[2-1](1)= (2.0*xi+1.0)*(eta*eta-eta)/4.0;
    t_shpders[2-1](2)= (xi*xi+xi  )*(2.0*eta-1.0)/4.0;
    t_shpders[2-1](3)= 0.0;

    t_shpvals[3-1]=(xi*xi+xi )*(eta*eta+eta)/4.0;
    t_shpders[3-1](1)= (2.0*xi+1.0)*(eta*eta+eta)/4.0;
    t_shpders[3-1](2)= (xi*xi+xi  )*(2.0*eta+1.0)/4.0;
    t_shpders[3-1](3)= 0.0;

    t_shpvals[4-1]=(xi*xi-xi )*(eta*eta+eta)/4.0;
    t_shpders[4-1](1)= (2.0*xi-1.0)*(eta*eta+eta)/4.0;
    t_shpders[4-1](2)= (xi*xi-xi  )*(2.0*eta+1.0)/4.0;
    t_shpders[4-1](3)= 0.0;

    t_shpvals[5-1]=(1.0-xi*xi)*(eta*eta-eta)/2.0;
    t_shpders[5-1](1)=-xi*(eta*eta-eta);
    t_shpders[5-1](2)= (1.0-xi*xi )*(2.0*eta-1.0)/2.0;
    t_shpders[5-1](3)= 0.0;

    t_shpvals[6-1]=(xi*xi+xi )*(1.0-eta*eta)/2.0;
    t_shpders[6-1](1)= (2.0*xi+1.0)*(1.0-eta*eta)/2.0;
    t_shpders[6-1](2)=-(xi*xi+xi )*eta;
    t_shpders[6-1](3)= 0.0;

    t_shpvals[7-1]=(1.0-xi*xi)*(eta*eta+eta)/2.0;
    t_shpders[7-1](1)=-xi*(eta*eta+eta);
    t_shpders[7-1](2)= (1.0-xi*xi )*(2.0*eta+1.0)/2.0;
    t_shpders[7-1](3)= 0.0;

    t_shpvals[8-1]=(xi*xi-xi )*(1.0-eta*eta)/2.0;
    t_shpders[8-1](1)= (2.0*xi-1.0)*(1.0-eta*eta)/2.0;
    t_shpders[8-1](2)=-(xi*xi-xi )*eta;
    t_shpders[8-1](3)= 0.0;

    t_shpvals[9-1]=(1.0-xi*xi)*(1.0-eta*eta);
    t_shpders[9-1](1)=-2.0*xi*(1.0-eta*eta);
    t_shpders[9-1](2)=-2.0*eta*(1.0-xi*xi);
    t_shpders[9-1](3)= 0.0;

}
