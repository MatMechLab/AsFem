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
//+++ Purpose: defines the user-defined (1) shape functions in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/ShapeFunUser1.h"

void ShapeFunUser1::calcShapeValsAndDerivatives(const double &xi,const double &eta,const double &zeta,vector<double> &t_shpvals,vector<Vector3d> &t_shpders){
    //****************************************
    //*** get rid of unused warning
    //****************************************
    if(xi||eta||zeta||t_shpvals.size()||t_shpders.size()){}

    //********************************************
    //*** please provide 
    //***   1) shape function value 
    //***   2) their local derivatives (w.r.t xi,eta,zeta not x,y,z !!!)
    //*** this subroutine will be called in ShapeFunUser class
    //********************************************

}