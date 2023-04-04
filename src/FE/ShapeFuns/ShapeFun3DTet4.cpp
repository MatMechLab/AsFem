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
//+++ Purpose: 3d tet4 shape function
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/ShapeFun3DTet4.h"

ShapeFun3DTet4::ShapeFun3DTet4(){

}

void ShapeFun3DTet4::calc3DShapeValsAndDerivatives(const double &xi,const double &eta,const double &zeta,vector<double> &t_shpvals,vector<Vector3d> &t_shpders){
    if(t_shpvals.size()<4 || t_shpders.size()<4){
        MessagePrinter::printErrorTxt("your shape val or derivs vector size is smaller than 4, error detected in ShapeFun3DTet4.cpp");
        MessagePrinter::exitAsFem();
    }

    // It should be mentioned that, tet4 mesh has different node ordering, here we use the one defined in Gmsh
    // https://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
    /**
     * VTK Cell type: vtkTetra
     * bottom layer:
     * 3
     * |\
     * | \
     * |  \
     * |   \
     * 1----2
     * 
     * top layer:
     * 4*
    */

    t_shpvals[0]=1.0-xi-eta-zeta;
    t_shpders[0](1)=-1.0;
    t_shpders[0](2)=-1.0;
    t_shpders[0](3)=-1.0;

    t_shpvals[1]=xi;
    t_shpders[1](1)=1.0;
    t_shpders[1](2)=0.0;
    t_shpders[1](3)=0.0;

    t_shpvals[2]=eta;
    t_shpders[2](1)= 0.0;
    t_shpders[2](2)= 1.0;
    t_shpders[2](3)= 0.0;

    t_shpvals[3]=zeta;
    t_shpders[3](1)= 0.0;
    t_shpders[3](2)= 0.0;
    t_shpders[3](3)= 1.0;

}