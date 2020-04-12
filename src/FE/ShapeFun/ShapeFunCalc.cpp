//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "FE/ShapeFun.h"

void ShapeFun::Calc(const double &xi,const Nodes &nodes,const bool &flag){
    Compute1DLagrangeShapeFun(xi,nodes,flag);
}

//*************************************
void ShapeFun::Calc(const double &xi,const double &eta,const Nodes &nodes,const bool &flag){
    Compute2DLagrangeShapeFun(xi,eta,nodes,flag);
}

//**************************************
void ShapeFun::Calc(const double &xi,const double &eta,const double &zeta,const Nodes &nodes,const bool &flag){
    Compute3DLagrangeShapeFun(xi,eta,zeta,nodes,flag);
}