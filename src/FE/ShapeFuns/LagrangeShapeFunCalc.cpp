//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.11.29
//+++ Purpose: implement the general FE shape functions for FEM calculation
//+++          in this code AsFem offer you:
//+++            1) lagrange shape function in 1d case,i.e. edge2,edge3,edge4  
//+++            2) lagrange shape function for 2d case, i.e. quad4,8,9 and tri3,6 mesh
//+++            3) lagrange shape function for 3d case, i.e. hex8,20,27 and tet4,10 mesh
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/LagrangeShapeFun.h"

void LagrangeShapeFun::Calc(const double &xi,const Nodes &nodes,const bool &flag){
    Compute1DLagrangeShapeFun(xi,nodes,flag);
}

//*************************************
void LagrangeShapeFun::Calc(const double &xi,const double &eta,const Nodes &nodes,const bool &flag){
    Compute2DLagrangeShapeFun(xi,eta,nodes,flag);
}

//**************************************
void LagrangeShapeFun::Calc(const double &xi,const double &eta,const double &zeta,const Nodes &nodes,const bool &flag){
    Compute3DLagrangeShapeFun(xi,eta,zeta,nodes,flag);
}