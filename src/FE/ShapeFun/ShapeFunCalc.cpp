//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/walkandthinker/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "FE/ShapeFun.h"

void ShapeFun::Calc(const double &xi,const Nodes &nodes){
    Compute1DLagrangeShapeFun(xi,nodes);
}

//*************************************
void ShapeFun::Calc(const double &xi,const double &eta,const Nodes &nodes){
    Compute2DLagrangeShapeFun(xi,eta,nodes);
}

//**************************************
void ShapeFun::Calc(const double &xi,const double &eta,const double &zeta,const Nodes &nodes){
    Compute3DLagrangeShapeFun(xi,eta,zeta,nodes);
}