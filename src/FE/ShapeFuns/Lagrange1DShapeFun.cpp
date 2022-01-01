//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.08.08
//+++ Purpose: implement the calculation for 1D Lagrange shape functions 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/Lagrange1DShapeFun.h"

Lagrange1DShapeFun::Lagrange1DShapeFun(){
    _dxdxi=0.0;_dydxi=0.0;_dzdxi=0.0;
    _nNodes=0;
    _tol=1.0e-15;
}
//******************************************************************
void Lagrange1DShapeFun::Calc1DShapeFun(const MeshType &meshtype, const double &xi, const Nodes &nodes, vector<double> &shape_val, vector<Vector3d> &shape_grad,double &detjac,bool flag){
    switch (meshtype){
        case MeshType::EDGE2:
            shape_val[ 0]=0.5*(1.0-xi);
            shape_grad[0].setZero();
            shape_grad[0](1)=-0.5;

            shape_val[ 1]=0.5*(1.0+xi);
            shape_grad[1].setZero();
            shape_grad[1](1)=0.5;

            _nNodes=2;
            break;
        case MeshType::EDGE3:
            shape_val[ 0]=0.5*xi*(xi-1.0);
            shape_grad[0].setZero();
            shape_grad[0](1)=0.5*(2.0*xi-1.0);

            shape_val[ 1]=-(xi + 1.0)*(xi - 1.0);
            shape_grad[1].setZero();
            shape_grad[1](1)=-2.0*xi;

            shape_val[ 2]=0.5*xi*(xi + 1.0);
            shape_grad[2].setZero();
            shape_grad[2](1)=0.5*(2.0*xi + 1.0);
            
            _nNodes=3;
            break;
        case MeshType::EDGE4:
            shape_val[ 0]=-(3.0*xi+1.0)*(3.0*xi-1.0)*(    xi-1.0)/16.0;
            shape_grad[0].setZero();
            shape_grad[0](1)=-27.0*xi*xi/16.0+9.0*xi/8.0+ 1.0/16.0;

            shape_val[ 1]= (3.0*xi+3.0)*(3.0*xi-1.0)*(3.0*xi-3.0)/16.0;
            shape_grad[1].setZero();
            shape_grad[1](1)=-27.0*xi*xi/16.0+9.0*xi/8.0+ 1.0/16.0;

            shape_val[ 2]=-(3.0*xi+3.0)*(3.0*xi+1.0)*(3.0*xi-3.0)/16.0;
            shape_grad[2].setZero();
            shape_grad[2](1)=-81.0*xi*xi/16.0-9.0*xi/8.0+27.0/16.0;
            
            shape_val[ 3]= (    xi+1.0)*(3.0*xi+1.0)*(3.0*xi-1.0)/16.0;
            shape_grad[3].setZero();
            shape_grad[3](1)= 27.0*xi*xi/16.0+9.0*xi/8.0- 1.0/16.0;

            _nNodes=4;
            break;
        default:
            MessagePrinter::PrintErrorTxt("unsupported mesh type in Lagrange1DShapeFun calculation");
            MessagePrinter::AsFem_Exit();
            break;
    }

    _dxdxi=0.0;_dydxi=0.0;_dzdxi=0.0;
    for(int i=1;i<=_nNodes;i++){
        _dxdxi+=shape_grad[i-1](1)*nodes(i,1);
        _dydxi+=shape_grad[i-1](1)*nodes(i,2);
        _dzdxi+=shape_grad[i-1](1)*nodes(i,3);
    }

    detjac=sqrt(_dxdxi*_dxdxi+_dydxi*_dydxi+_dzdxi*_dzdxi);

    if(abs(detjac)<1.0e-15){
        MessagePrinter::PrintErrorTxt("singular element in 1D case, this error occurs in your 1D shape function calculation");
        MessagePrinter::AsFem_Exit();
    }
    if(flag){
        // for the derivatives in global coordinates (x,y,z)
        for(int i=1;i<=_nNodes;i++){
            _val=shape_grad[i-1](1)/detjac;// dN/dX
            shape_grad[i-1](1)=_val*_dxdxi/detjac; // x-component
            shape_grad[i-1](2)=_val*_dydxi/detjac; // y-component
            shape_grad[i-1](3)=_val*_dzdxi/detjac; // z-component
        }
    }

}
