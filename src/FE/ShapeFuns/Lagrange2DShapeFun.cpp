//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.08.09
//+++ Purpose: implement the calculation for 2D Lagrange shape functions 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/Lagrange2DShapeFun.h"


Lagrange2DShapeFun::Lagrange2DShapeFun(){
    _dxdxi=0.0;_dxdeta=0.0;
    _dydxi=0.0;_dydeta=0.0;
    _dzdxi=0.0;_dzdeta=0.0;

    _nNodes=0;
    _tol=1.0e-15;

    _Jac2.Resize(2,2);_XJac2.Resize(2,2);
    _Vec32.Resize(3,2);_Vec21.Resize(2,1);
    _dN.Resize(3,1);_dN2.Resize(2,1);
}
//**********************************************
void Lagrange2DShapeFun::Calc2DShapeFun(const MeshType &meshtype, const double &xi, const double &eta, const Nodes &nodes, vector<double> &shape_val, vector<Vector3d> &shape_grad, double &detjac,bool flag){
    switch (meshtype) {
        case MeshType::TRI3:
            _nNodes=3;
            shape_val[ 0]=1.0-xi-eta;
            shape_grad[0](1)=-1.0;
            shape_grad[0](2)=-1.0;
            shape_grad[0](3)=0.0;

            shape_val[ 1]=xi;
            shape_grad[1](1)=1.0;
            shape_grad[1](2)=0.0;
            shape_grad[1](3)=0.0;

            shape_val[ 2]=eta;
            shape_grad[2](1)=0.0;
            shape_grad[2](2)=1.0;
            shape_grad[2](3)=0.0;
            break;
        case MeshType::TRI6:
            // taken from: http://www.sd.ruhr-uni-bochum.de/downloads/Shape_funct.pdf
            _nNodes=6;
            shape_val[ 0]= (1.0-xi-eta)*(1.0-2*xi-2*eta);
            shape_grad[0](1)=-3.0+4.0*eta+4.0*xi;
            shape_grad[0](2)=-3.0+4.0*eta+4.0*xi;
            shape_grad[0](3)=0.0;

            shape_val[ 1]= xi*(2.0*xi-1.0);
            shape_grad[1](1)= 4.0*xi-1.0;
            shape_grad[1](2)= 0.0;
            shape_grad[1](3)=0.0;

            shape_val[ 2]= eta*(2.0*eta-1.0);
            shape_grad[2](1)= 0.0;
            shape_grad[2](2)= 4.0*eta-1.0;
            shape_grad[2](3)=0.0;

            shape_val[ 3]= 4.0*xi*(1.0-xi-eta);
            shape_grad[3](1)= 4.0*(1.0-eta-2*xi);
            shape_grad[3](2)=-4.0*xi;
            shape_grad[3](3)=0.0;

            shape_val[ 4]= 4.0*xi*eta;
            shape_grad[4](1)= 4.0*eta;
            shape_grad[4](2)= 4.0*xi;
            shape_grad[4](3)=0.0;

            shape_val[ 5]= 4.0*eta*(1.0-xi-eta);
            shape_grad[5](1)=-4.0*eta;
            shape_grad[5](2)= 4.0*(1-2*eta-xi);
            shape_grad[5](3)=0.0;
            break;
        case MeshType::QUAD4:
            _nNodes=4;
            shape_val[0]=(1.0-xi)*(1.0-eta)/4.0;
            shape_val[1]=(1.0+xi)*(1.0-eta)/4.0;
            shape_val[2]=(1.0+xi)*(1.0+eta)/4.0;
            shape_val[3]=(1.0-xi)*(1.0+eta)/4.0;

            shape_grad[0](1)= (eta-1.0)/4.0;
            shape_grad[0](2)= (xi -1.0)/4.0;
            shape_grad[0](3)=0.0;

            shape_grad[1](1)= (1.0-eta)/4.0;
            shape_grad[1](2)=-(1.0+xi )/4.0;
            shape_grad[1](3)=0.0;

            shape_grad[2](1)= (1.0+eta)/4.0;
            shape_grad[2](2)= (1.0+xi )/4.0;
            shape_grad[2](3)=0.0;

            shape_grad[3](1)=-(1.0+eta)/4.0;
            shape_grad[3](2)= (1.0-xi )/4.0;
            shape_grad[3](3)=0.0;
            break;
        case MeshType::QUAD8:
            // 2D-8Nodes rectangle element
            _nNodes=8;
            shape_val[0]=(1.0-xi)*(1.0-eta)*(-xi-eta-1.0)/4.0;
            shape_val[1]=(1.0+xi)*(1.0-eta)*( xi-eta-1.0)/4.0;
            shape_val[2] =(1.0+xi)*(1.0+eta)*( xi+eta-1.0)/4.0;
            shape_val[3]=(1.0-xi)*(1.0+eta)*(-xi+eta-1.0)/4.0;
            shape_val[4]=(1.0-xi*xi)*(1.0-eta    )/2.0;
            shape_val[5]=(1.0+xi   )*(1.0-eta*eta)/2.0;
            shape_val[6]=(1.0-xi*xi)*(1.0+eta    )/2.0;
            shape_val[7]=(1.0-xi   )*(1.0-eta*eta)/2.0;

            // derivatives over xi and eta
            shape_grad[0](1)=(1.0-eta)*(2.0*xi+eta)/4.0;
            shape_grad[0](2)=(1.0-xi )*(xi+2.0*eta)/4.0;
            shape_grad[0](3)=0.0;

            shape_grad[1](1)=(1.0-eta)*(2.0*xi-eta)/4.0;
            shape_grad[1](2)=(1.0+xi )*(2.0*eta-xi)/4.0;
            shape_grad[1](3)=0.0;

            shape_grad[2](1)=(1.0+eta)*(2.0*xi+eta)/4.0;
            shape_grad[2](2)=(1.0+xi )*(xi+2.0*eta)/4.0;
            shape_grad[2](3)=0.0;

            shape_grad[3](1)=(1.0+eta)*(2.0*xi-eta)/4.0;
            shape_grad[3](2)=(1.0-xi )*(2.0*eta-xi)/4.0;
            shape_grad[3](3)=0.0;

            shape_grad[4](1)=xi*(eta-1.0);
            shape_grad[4](2)=(xi*xi-1.0)/2.0;
            shape_grad[4](3)=0.0;

            shape_grad[5](1)=(1.0-eta*eta)/2.0;
            shape_grad[5](2)=-(1.0+xi)*eta;
            shape_grad[5](3)=0.0;

            shape_grad[6](1)=-xi*(1.0+eta);
            shape_grad[6](2)=(1.0-xi*xi)/2.0;
            shape_grad[6](3)=0.0;

            shape_grad[7](1)=(eta*eta-1.0)/2.0;
            shape_grad[7](2)=(xi-1.0)*eta;
            shape_grad[7](3)=0.0;

            break;
        case MeshType::QUAD9:
            // 2D-9Nodes rectangle element
            _nNodes=9;
            shape_val[0]=(xi*xi-xi )*(eta*eta-eta)/4.0;
            shape_val[1]=(xi*xi+xi )*(eta*eta-eta)/4.0;
            shape_val[2]=(xi*xi+xi )*(eta*eta+eta)/4.0;
            shape_val[3]=(xi*xi-xi )*(eta*eta+eta)/4.0;
            shape_val[4]=(1.0-xi*xi)*(eta*eta-eta)/2.0;
            shape_val[5]=(xi*xi+xi )*(1.0-eta*eta)/2.0;
            shape_val[6]=(1.0-xi*xi)*(eta*eta+eta)/2.0;
            shape_val[7]=(xi*xi-xi )*(1.0-eta*eta)/2.0;
            shape_val[8]=(1.0-xi*xi)*(1.0-eta*eta);

            shape_grad[0](1)=(2.0*xi-1.0)*(eta*eta-eta)/4.0;
            shape_grad[0](2)=(xi*xi-xi  )*(2.0*eta-1.0)/4.0;
            shape_grad[0](3)=0.0;

            shape_grad[1](1)=(2.0*xi+1.0)*(eta*eta-eta)/4.0;
            shape_grad[1](2)=(xi*xi+xi  )*(2.0*eta-1.0)/4.0;
            shape_grad[1](3)=0.0;

            shape_grad[2](1)=(2.0*xi+1.0)*(eta*eta+eta)/4.0;
            shape_grad[2](2)=(xi*xi+xi  )*(2.0*eta+1.0)/4.0;
            shape_grad[2](3)=0.0;

            shape_grad[3](1)=(2.0*xi-1.0)*(eta*eta+eta)/4.0;
            shape_grad[3](2)=(xi*xi-xi  )*(2.0*eta+1.0)/4.0;
            shape_grad[3](3)=0.0;

            shape_grad[4](1)=-xi*(eta*eta-eta);
            shape_grad[4](2)=(1.0-xi*xi )*(2.0*eta-1.0)/2.0;
            shape_grad[4](3)=0.0;

            shape_grad[5](1)=(2.0*xi+1.0)*(1.0-eta*eta)/2.0;
            shape_grad[5](2)=-(xi*xi+xi )*eta;
            shape_grad[5](3)=0.0;

            shape_grad[6](1)=-xi*(eta*eta+eta);
            shape_grad[6](2)=(1.0-xi*xi )*(2.0*eta+1.0)/2.0;
            shape_grad[6](3)=0.0;

            shape_grad[7](1)=(2.0*xi-1.0)*(1.0-eta*eta)/2.0;
            shape_grad[7](2)=-(xi*xi-xi )*eta;
            shape_grad[7](3)=0.0;

            shape_grad[8](1)=-2.0*xi*(1.0-eta*eta);
            shape_grad[8](2)=-2.0*eta*(1.0-xi*xi);
            shape_grad[8](3)=0.0;
            break;
        default:
            MessagePrinter::PrintErrorTxt("unsupported mesh type in Lagrange2DShapeFun.cpp");
            MessagePrinter::AsFem_Exit();
            break;
    }

    _dxdxi=0.0;_dydxi=0.0;_dzdxi=0.0;
    _dxdeta=0.0;_dydeta=0.0;_dzdeta=0.0;
    for(int i=1;i<=_nNodes;i++){
        _dxdxi+=shape_grad[i-1](1)*nodes(i,1);
        _dydxi+=shape_grad[i-1](1)*nodes(i,2);
        _dzdxi+=shape_grad[i-1](1)*nodes(i,3);

        _dxdeta+=shape_grad[i-1](2)*nodes(i,1);
        _dydeta+=shape_grad[i-1](2)*nodes(i,2);
        _dzdeta+=shape_grad[i-1](2)*nodes(i,3);
    }
    _Vec32(1,1)=_dxdxi;_Vec32(1,2)=_dxdeta;
    _Vec32(2,1)=_dydxi;_Vec32(2,2)=_dydeta;
    _Vec32(3,1)=_dzdxi;_Vec32(3,2)=_dzdeta;

    // here the jacobian must contains the contribution from z-axis, even though
    // you only have a "2d" elemet!!!
    detjac=(_dydxi*_dzdeta-_dydeta*_dzdxi)*(_dydxi*_dzdeta-_dydeta*_dzdxi)
           +(_dzdxi*_dxdeta-_dzdeta*_dxdxi)*(_dzdxi*_dxdeta-_dzdeta*_dxdxi)
           +(_dxdxi*_dydeta-_dxdeta*_dydxi)*(_dxdxi*_dydeta-_dxdeta*_dydxi);

    detjac=sqrt(detjac);

    if(detjac<_tol){
        MessagePrinter::PrintErrorTxt("singular element in 2D case, this error occurs in your 2D shape function calculation");
        MessagePrinter::AsFem_Exit();
    }

    if(flag){

        double vx,vy,vz;
        vx=_dydxi*_dzdeta-_dydeta*_dzdxi/detjac;
        vy=_dzdxi*_dxdeta-_dzdeta*_dxdxi/detjac;
        vz=_dxdxi*_dydeta-_dxdeta*_dydxi/detjac;

        if(abs(vx)<_tol && abs(vy)<_tol && abs(vz)>_tol){
            // for x-y plane
            _Jac2(1,1)=_dxdxi ;_Jac2(1,2)=_dydxi;
            _Jac2(2,1)=_dxdeta;_Jac2(2,2)=_dydeta;

            _XJac2=_Jac2.Inverse();
            for(int i=1;i<=_nNodes;i++){
                _dN2(1,1)=shape_grad[i-1](1)*_XJac2(1,1)+shape_grad[i-1](2)*_XJac2(1,2);
                _dN2(2,1)=shape_grad[i-1](1)*_XJac2(2,1)+shape_grad[i-1](2)*_XJac2(2,2);

                shape_grad[i-1](1)=_dN2(1,1);
                shape_grad[i-1](2)=_dN2(2,1);
                shape_grad[i-1](3)=0.0;
            }
        }
        else if(abs(vx)<_tol && abs(vz)<_tol && abs(vy)>_tol){
            // for z-x plane
            _Jac2(1,1)=_dzdxi ;_Jac2(1,2)=_dxdxi;
            _Jac2(2,1)=_dzdeta;_Jac2(2,2)=_dxdeta;

            _XJac2=_Jac2.Inverse();
            for(int i=1;i<=_nNodes;i++){
                _dN2(1,1)=shape_grad[i-1](1)*_XJac2(1,1)+shape_grad[i-1](2)*_XJac2(1,2);
                _dN2(2,1)=shape_grad[i-1](1)*_XJac2(2,1)+shape_grad[i-1](2)*_XJac2(2,2);

                shape_grad[i-1](1)=_dN2(2,1);
                shape_grad[i-1](2)=0.0;
                shape_grad[i-1](3)=_dN2(1,1);
            }
        }
        else if(vy<_tol && vz<_tol && vx>_tol){
            // for y-z plane
            _Jac2(1,1)=_dydxi ;_Jac2(1,2)=_dzdxi;
            _Jac2(2,1)=_dydeta;_Jac2(2,2)=_dzdeta;

            _XJac2=_Jac2.Inverse();
            for(int i=1;i<=_nNodes;i++){
                _dN2(1,1)=shape_grad[i-1](1)*_XJac2(1,1)+shape_grad[i-1](2)*_XJac2(1,2);
                _dN2(2,1)=shape_grad[i-1](1)*_XJac2(2,1)+shape_grad[i-1](2)*_XJac2(2,2);

                shape_grad[i-1](1)=0.0;
                shape_grad[i-1](2)=_dN2(1,1);
                shape_grad[i-1](3)=_dN2(2,1);
            }
        }
        else{
            _Jac2=_Vec32.Transpose()*_Vec32;

            _XJac2=_Jac2.Inverse();
        
            for(int i=1;i<=_nNodes;i++){
                _dN2(1,1)=shape_grad[i-1](1)*_XJac2(1,1)+shape_grad[i-1](2)*_XJac2(1,2);
                _dN2(2,1)=shape_grad[i-1](1)*_XJac2(2,1)+shape_grad[i-1](2)*_XJac2(2,2);

                _dN=_Vec32*_dN2;

                shape_grad[i-1](1)=_dN(1,1);
                shape_grad[i-1](2)=_dN(2,1);
                shape_grad[i-1](3)=_dN(3,1);
            }
            
        }
        
    }
}
