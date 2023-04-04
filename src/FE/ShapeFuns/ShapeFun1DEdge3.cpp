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
//+++ Purpose: implement the edge3 shape fun and its local derivatives
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/ShapeFun1DEdge3.h"

ShapeFun1DEdge3::ShapeFun1DEdge3(){
}

void ShapeFun1DEdge3::calc1DShapeValsAndDerivatives(const double &xi,vector<double> &t_shpvals,vector<Vector3d> &t_shpders){
    if(t_shpvals.size()<3 || t_shpders.size()<3){
        MessagePrinter::printErrorTxt("your shape val or derivs vector size is smaller than 3, error detected in Shape1DEdge3.cpp");
        MessagePrinter::exitAsFem();
    }
    /**
     * The nodes should look like:
     * 1---2---3
    */
    t_shpvals[1-1]=0.5*xi*(xi-1.0);
    t_shpders[1-1](1)= 0.5*(2.0*xi-1.0);
    t_shpders[1-1](2)= 0.0;
    t_shpders[1-1](3)= 0.0;

    t_shpvals[2-1]=-(xi+1.0)*(xi-1.0);
    t_shpders[2-1](1)=-2.0*xi;
    t_shpders[2-1](2)= 0.0;
    t_shpders[2-1](3)= 0.0;

    t_shpvals[3-1]=0.5*xi*(xi+1.0);
    t_shpders[3-1](1)= 0.5*(2.0*xi+1.0);
    t_shpders[3-1](2)= 0.0;
    t_shpders[3-1](3)= 0.0;
    
}