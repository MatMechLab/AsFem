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
//+++ Date   : 2022.05.13
//+++ Purpose: implement the edge2 shape fun and its local derivatives
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/ShapeFun1DEdge2.h"

ShapeFun1DEdge2::ShapeFun1DEdge2(){

}

void ShapeFun1DEdge2::calc1DShapeValsAndDerivatives(const double &xi,vector<double> &t_shpvals,vector<Vector3d> &t_shpders){
    if(t_shpvals.size()<2 || t_shpders.size()<2){
        MessagePrinter::printErrorTxt("your shape val or derivs vector size is smaller than 2, error detected in Shape1DEdge2.cpp");
        MessagePrinter::exitAsFem();
    }
    /**
     * The nodes should look like:
     * 1---2
    */
    t_shpvals[1-1]=0.5*(1.0-xi);
    t_shpders[1-1](1)=-0.5;
    t_shpders[1-1](2)= 0.0;
    t_shpders[1-1](3)= 0.0;

    t_shpvals[2-1]=0.5*(1.0+xi);
    t_shpders[2-1](1)= 0.5;
    t_shpders[2-1](2)= 0.0;
    t_shpders[2-1](3)= 0.0;
    
}