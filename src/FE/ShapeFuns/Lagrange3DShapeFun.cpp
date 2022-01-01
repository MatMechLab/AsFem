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
//+++ Date   : 2021.08.11
//+++ Purpose: implement the calculation for 3D Lagrange shape functions 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/Lagrange3DShapeFun.h"

Lagrange3DShapeFun::Lagrange3DShapeFun(){
    _dxdxi=0.0;_dxdeta=0.0;_dxdzeta=0.0;
    _dydxi=0.0;_dydeta=0.0;_dydzeta=0.0;
    _dzdxi=0.0;_dzdeta=0.0;_dzdzeta=0.0;

    _Jac.Resize(3,3);_XJac.Resize(3,3);
}
//************************************************
void Lagrange3DShapeFun::Calc3DShapeFun(const MeshType &meshtype, const double &xi, const double &eta, const double &zeta, const Nodes &nodes, vector<double> &shape_val, vector<Vector3d> &shape_grad, double &detjac,bool flag){
    switch (meshtype) {
        case MeshType::TET4:{
            _nNodes=4;
            const double sqrt2=sqrt(2.0);
            const double sqrt3=sqrt(3.0);

            shape_val[ 0]=(3.0+8.0*xi-2.0*sqrt2*zeta)/12.0;
            shape_grad[0](1)= 8.0/12.0;
            shape_grad[0](2)= 0.0;
            shape_grad[0](3)=-2.0*sqrt2/12.0;

            shape_val[ 1 ]=(3.0-4.0*xi-4.0*sqrt3*eta-2.0*sqrt2*zeta)/12.0;
            shape_grad[1](1)=-4.0/12.0;
            shape_grad[1](2)=-4.0*sqrt3/12.0;
            shape_grad[1](3)=-2.0*sqrt2/12.0;

            shape_val[ 2]=(3.0-4.0*xi+4.0*sqrt3*eta-2.0*sqrt2*zeta)/12.0;
            shape_grad[2](1)=-4.0/12.0;
            shape_grad[2](2)= 4.0*sqrt3/12.0;
            shape_grad[2](3)=-2.0*sqrt2/12.0;

            shape_val[ 3]=(1.0+2.0*sqrt2*zeta)/4.0;
            shape_grad[3](1)=0.0;
            shape_grad[3](2)=0.0;
            shape_grad[3](3)=2.0*sqrt2/4.0;
            break;
        }
        case MeshType::HEX8:{
            _nNodes=8;
            shape_val[ 0] = (1 - xi) * (1 - eta) * (1 - zeta) / 8.0;
            shape_grad[0](1) = -(1 - eta) * (1 - zeta) / 8.0;
            shape_grad[0](2) = -(1 - xi) * (1 - zeta) / 8.0;
            shape_grad[0](3) = -(1 - xi) * (1 - eta) / 8.0;

            shape_val[ 1] = (1 + xi) * (1 - eta) * (1 - zeta) / 8.0;
            shape_grad[1](1) = (1 - eta) * (1 - zeta) / 8.0;
            shape_grad[1](2) = -(1 + xi) * (1 - zeta) / 8.0;
            shape_grad[1](3) = -(1 + xi) * (1 - eta) / 8.0;

            shape_val[ 2] = (1 + xi) * (1 + eta) * (1 - zeta) / 8.0;
            shape_grad[2](1) = (1 + eta) * (1 - zeta) / 8.0;
            shape_grad[2](2) = (1 + xi) * (1 - zeta) / 8.0;
            shape_grad[2](3) = -(1 + xi) * (1 + eta) / 8.0;

            shape_val[ 3] = (1 - xi) * (1 + eta) * (1 - zeta) / 8.0;
            shape_grad[3](1) = -(1 + eta) * (1 - zeta) / 8.0;
            shape_grad[3](2) = (1 - xi) * (1 - zeta) / 8.0;
            shape_grad[3](3) = -(1 - xi) * (1 + eta) / 8.0;

            shape_val[ 4] = (1 - xi) * (1 - eta) * (1 + zeta) / 8.0;
            shape_grad[4](1) = -(1 - eta) * (1 + zeta) / 8.0;
            shape_grad[4](2) = -(1 - xi) * (1 + zeta) / 8.0;
            shape_grad[4](3) = (1 - xi) * (1 - eta) / 8.0;

            shape_val[ 5] = (1 + xi) * (1 - eta) * (1 + zeta) / 8.0;
            shape_grad[5](1) = (1 - eta) * (1 + zeta) / 8.0;
            shape_grad[5](2) = -(1 + xi) * (1 + zeta) / 8.0;
            shape_grad[5](3) = (1 + xi) * (1 - eta) / 8.0;

            shape_val[ 6] = (1 + xi) * (1 + eta) * (1 + zeta) / 8.0;
            shape_grad[6](1) = (1 + eta) * (1 + zeta) / 8.0;
            shape_grad[6](2) = (1 + xi) * (1 + zeta) / 8.0;
            shape_grad[6](3) = (1 + xi) * (1 + eta) / 8.0;

            shape_val[ 7] = (1 - xi) * (1 + eta) * (1 + zeta) / 8.0;
            shape_grad[7](1) = -(1 + eta) * (1 + zeta) / 8.0;
            shape_grad[7](2) = (1 - xi) * (1 + zeta) / 8.0;
            shape_grad[7](3) = (1 - xi) * (1 + eta) / 8.0;
            break;
        }
        case MeshType::HEX20:{
            _nNodes=20;
            const double XI[]={0.0,
                             -1.0,1.0,1.0,-1.0,//1-4
                             -1.0,1.0,1.0,-1.0,//5-8
                              0.0,1.0,0.0,-1.0,//9-12
                              0.0,1.0,0.0,-1.0,//13-16
                             -1.0,1.0,1.0,-1.0//17-20
                             };
            const double ETA[]={0.0,
                               -1.0,-1.0,1.0,1.0,//1-4
                               -1.0,-1.0,1.0,1.0,//5-8
                               -1.0, 0.0,1.0,0.0,//9-12
                               -1.0, 0.0,1.0,0.0,//13-16
                               -1.0,-1.0,1.0,1.0//17-20
                               };
            const double ZETA[]={0.0,
                                -1.0,-1.0,-1.0,-1.0,//1-4
                                 1.0, 1.0, 1.0, 1.0,//5-8
                                -1.0,-1.0,-1.0,-1.0,//9-12
                                 1.0, 1.0, 1.0, 1.0,//13-16
                                 0.0, 0.0, 0.0, 0.0//17-20
                                };

            int i;
            // for corner nodes
            for(i=1;i<=8;++i){
                shape_val[i-1]=(1.0+XI[i]*xi)*(1.0+ETA[i]*eta)*(1.0+ZETA[i]*zeta)*(XI[i]*xi+ETA[i]*eta+ZETA[i]*zeta-2.0)/8.0;

                shape_grad[i-1](1)=(XI[i])*(1.0+ETA[i]*eta)*(1.0+ZETA[i]*zeta)*(XI[i]*xi+ETA[i]*eta+ZETA[i]*zeta-2.0)/8.0
                            +(1.0+XI[i]*xi)*(1.0+ETA[i]*eta)*(1.0+ZETA[i]*zeta)*(XI[i])/8.0;

                shape_grad[i-1](2)=(1.0+XI[i]*xi)*(ETA[i])*(1.0+ZETA[i]*zeta)*(XI[i]*xi+ETA[i]*eta+ZETA[i]*zeta-2.0)/8.0
                            +(1.0+XI[i]*xi)*(1.0+ETA[i]*eta)*(1.0+ZETA[i]*zeta)*(ETA[i])/8.0;

                shape_grad[i-1](3)=(1.0+XI[i]*xi)*(1.0+ETA[i]*eta)*(ZETA[i])*(XI[i]*xi+ETA[i]*eta+ZETA[i]*zeta-2.0)/8.0
                            +(1.0+XI[i]*xi)*(1.0+ETA[i]*eta)*(1.0+ZETA[i]*zeta)*(ZETA[i])/8.0;
            }

            // for midside nodes
            for(i=1;i<=4;++i){
                // for 9,11,13,15
                shape_val[ 8+2*i-1-1]=(1.0-xi*xi)*(1.0+ETA[8+2*i-1]*eta)*(1.0+ZETA[8+2*i-1]*zeta)/4.0;
                shape_grad[8+2*i-1-1](1)=(-2.0*xi)*(1.0+ETA[8+2*i-1]*eta)*(1.0+ZETA[8+2*i-1]*zeta)/4.0;
                shape_grad[8+2*i-1-1](2)=(1.0-xi*xi)*(ETA[8+2*i-1])*(1.0+ZETA[8+2*i-1]*zeta)/4.0;
                shape_grad[8+2*i-1-1](3)=(1.0-xi*xi)*(1.0+ETA[8+2*i-1]*eta)*(ZETA[8+2*i-1])/4.0;
 
                // for 10,12,14,16
                shape_val[ 8+2*i-1]=(1.0-eta*eta)*(1.0+XI[8+2*i]*xi)*(1.0+ZETA[8+2*i]*zeta)/4.0;
                shape_grad[8+2*i-1](1)=(1.0-eta*eta)*(XI[8+2*i])*(1.0+ZETA[8+2*i]*zeta)/4.0;
                shape_grad[8+2*i-1](2)=(-2.0*eta)*(1.0+XI[8+2*i]*xi)*(1.0+ZETA[8+2*i]*zeta)/4.0;
                shape_grad[8+2*i-1](3)=(1.0-eta*eta)*(1.0+XI[8+2*i]*xi)*(ZETA[8+2*i])/4.0;

                // for 17,18,19,20
                shape_val[ 16+i-1]=(1.0-zeta*zeta)*(1.0+XI[16+i]*xi)*(1.0+ETA[16+i]*eta)/4.0;
                shape_grad[16+i-1](1)=(1.0-zeta*zeta)*(XI[16+i])*(1.0+ETA[16+i]*eta)/4.0;
                shape_grad[16+i-1](2)=(1.0-zeta*zeta)*(1.0+XI[16+i]*xi)*(ETA[16+i])/4.0;
                shape_grad[16+i-1](3)=(-2.0*zeta)*(1.0+XI[16+i]*xi)*(1.0+ETA[16+i]*eta)/4.0;
            }
            break;
        }
        case MeshType::HEX27:{
            _nNodes=27;
            shape_val[ 1-1] = xi * (xi - 1) * eta * (eta - 1) * zeta * (zeta - 1) / 8.0;
            shape_grad[1-1](1) = (2 * xi - 1) * eta * (eta - 1) * zeta * (zeta - 1) / 8.0;
            shape_grad[1-1](2) = xi * (xi - 1) * (2 * eta - 1) * zeta * (zeta - 1) / 8.0;
            shape_grad[1-1](3) = xi * (xi - 1) * eta * (eta - 1) * (2 * zeta - 1) / 8.0;

            shape_val[ 2-1] = xi * (xi + 1) * eta * (eta - 1) * zeta * (zeta - 1) / 8.0;
            shape_grad[2-1](1) = (2 * xi + 1) * eta * (eta - 1) * zeta * (zeta - 1) / 8.0;
            shape_grad[2-1](2) = xi * (xi + 1) * (2 * eta - 1) * zeta * (zeta - 1) / 8.0;
            shape_grad[2-1](3) = xi * (xi + 1) * eta * (eta - 1) * (2 * zeta - 1) / 8.0;

            shape_val[ 3-1] = xi * (xi + 1) * eta * (eta + 1) * zeta * (zeta - 1) / 8.0;
            shape_grad[3-1](1) = (2 * xi + 1) * eta * (eta + 1) * zeta * (zeta - 1) / 8.0;
            shape_grad[3-1](2) = xi * (xi + 1) * (2 * eta + 1) * zeta * (zeta - 1) / 8.0;
            shape_grad[3-1](3) = xi * (xi + 1) * eta * (eta + 1) * (2 * zeta - 1) / 8.0;

            shape_val[ 4-1] = xi * (xi - 1) * eta * (eta + 1) * zeta * (zeta - 1) / 8.0;
            shape_grad[4-1](1) = (2 * xi - 1) * eta * (eta + 1) * zeta * (zeta - 1) / 8.0;
            shape_grad[4-1](2) = xi * (xi - 1) * (2 * eta + 1) * zeta * (zeta - 1) / 8.0;
            shape_grad[4-1](3) = xi * (xi - 1) * eta * (eta + 1) * (2 * zeta - 1) / 8.0;

            shape_val[ 5-1] = xi * (xi - 1) * eta * (eta - 1) * zeta * (zeta + 1) / 8.0;
            shape_grad[5-1](1) = (2 * xi - 1) * eta * (eta - 1) * zeta * (zeta + 1) / 8.0;
            shape_grad[5-1](2) = xi * (xi - 1) * (2 * eta - 1) * zeta * (zeta + 1) / 8.0;
            shape_grad[5-1](3) = xi * (xi - 1) * eta * (eta - 1) * (2 * zeta + 1) / 8.0;

            shape_val[ 6-1] = xi * (xi + 1) * eta * (eta - 1) * zeta * (zeta + 1) / 8.0;
            shape_grad[6-1](1) = (2 * xi + 1) * eta * (eta - 1) * zeta * (zeta + 1) / 8.0;
            shape_grad[6-1](2) = xi * (xi + 1) * (2 * eta - 1) * zeta * (zeta + 1) / 8.0;
            shape_grad[6-1](3) = xi * (xi + 1) * eta * (eta - 1) * (2 * zeta + 1) / 8.0;

            shape_val[ 7-1] = xi * (xi + 1) * eta * (eta + 1) * zeta * (zeta + 1) / 8.0;
            shape_grad[7-1](1) = (2 * xi + 1) * eta * (eta + 1) * zeta * (zeta + 1) / 8.0;
            shape_grad[7-1](2) = xi * (xi + 1) * (2 * eta + 1) * zeta * (zeta + 1) / 8.0;
            shape_grad[7-1](3) = xi * (xi + 1) * eta * (eta + 1) * (2 * zeta + 1) / 8.0;

            shape_val[ 8-1] = xi * (xi - 1) * eta * (eta + 1) * zeta * (zeta + 1) / 8.0;
            shape_grad[8-1](1) = (2 * xi - 1) * eta * (eta + 1) * zeta * (zeta + 1) / 8.0;
            shape_grad[8-1](2) = xi * (xi - 1) * (2 * eta + 1) * zeta * (zeta + 1) / 8.0;
            shape_grad[8-1](3) = xi * (xi - 1) * eta * (eta + 1) * (2 * zeta + 1) / 8.0;

            shape_val[ 9-1] = (1 - xi * xi) * eta * (eta - 1) * zeta * (zeta - 1) / 4.0;
            shape_grad[9-1](1) = -xi * eta * (eta - 1) * zeta * (zeta - 1) / 2.0;
            shape_grad[9-1](2) = (1 - xi * xi) * (2 * eta - 1) * zeta * (zeta - 1) / 4.0;
            shape_grad[9-1](3) = (1 - xi * xi) * eta * (eta - 1) * (2 * zeta - 1) / 4.0;

            shape_val[ 10-1] = xi * (xi + 1) * (1 - eta * eta) * zeta * (zeta - 1) / 4.0;
            shape_grad[10-1](1) = (2 * xi + 1) * (1 - eta * eta) * zeta * (zeta - 1) / 4.0;
            shape_grad[10-1](2) = -xi * (xi + 1) * eta * zeta * (zeta - 1) / 2.0;
            shape_grad[10-1](3) = xi * (xi + 1) * (1 - eta * eta) * (2 * zeta - 1) / 4.0;

            shape_val[ 11-1] = (1 - xi * xi) * eta * (eta + 1) * zeta * (zeta - 1) / 4.0;
            shape_grad[11-1](1) = -xi * eta * (eta + 1) * zeta * (zeta - 1) / 2.0;
            shape_grad[11-1](2) = (1 - xi * xi) * (2 * eta + 1) * zeta * (zeta - 1) / 4.0;
            shape_grad[11-1](3) = (1 - xi * xi) * eta * (eta + 1) * (2 * zeta - 1) / 4.0;

            shape_val[ 12-1] = xi * (xi - 1) * (1 - eta * eta) * zeta * (zeta - 1) / 4.0;
            shape_grad[12-1](1) = (2 * xi - 1) * (1 - eta * eta) * zeta * (zeta - 1) / 4.0;
            shape_grad[12-1](2) = -xi * (xi - 1) * eta * zeta * (zeta - 1) / 2.0;
            shape_grad[12-1](3) = xi * (xi - 1) * (1 - eta * eta) * (2 * zeta - 1) / 4.0;

            shape_val[ 13-1] = (1 - xi * xi) * eta * (eta - 1) * zeta * (zeta + 1) / 4.0;
            shape_grad[13-1](1) = -xi * eta * (eta - 1) * zeta * (zeta + 1) / 2.0;
            shape_grad[13-1](2) = (1 - xi * xi) * (2 * eta - 1) * zeta * (zeta + 1) / 4.0;
            shape_grad[13-1](3) = (1 - xi * xi) * eta * (eta - 1) * (2 * zeta + 1) / 4.0;

            shape_val[ 14-1] = xi * (xi + 1) * (1 - eta * eta) * zeta * (zeta + 1) / 4.0;
            shape_grad[14-1](1) = (2 * xi + 1) * (1 - eta * eta) * zeta * (zeta + 1) / 4.0;
            shape_grad[14-1](2) =-xi * (xi + 1) * eta * zeta * (zeta + 1) / 2.0;
            shape_grad[14-1](3) = xi * (xi + 1) * (1 - eta * eta) * (2 * zeta + 1) / 4.0;

            shape_val[ 15-1] = (1 - xi * xi) * eta * (eta + 1) * zeta * (zeta + 1) / 4.0;
            shape_grad[15-1](1) = -xi * eta * (eta + 1) * zeta * (zeta + 1) / 2.0;
            shape_grad[15-1](2) = (1 - xi * xi) * (2 * eta + 1) * zeta * (zeta + 1) / 4.0;
            shape_grad[15-1](3) = (1 - xi * xi) * eta * (eta + 1) * (2 * zeta + 1) / 4.0;

            shape_val[ 16-1] = xi * (xi - 1) * (1 - eta * eta) * zeta * (zeta + 1) / 4.0;
            shape_grad[16-1](1) = (2 * xi - 1) * (1 - eta * eta) * zeta * (zeta + 1) / 4.0;
            shape_grad[16-1](2) = -xi * (xi - 1) * eta * zeta * (zeta + 1) / 2.0;
            shape_grad[16-1](3) = xi * (xi - 1) * (1 - eta * eta) * (2 * zeta + 1) / 4.0;

            shape_val[ 17-1] = xi * (xi - 1) * eta * (eta - 1) * (1 - zeta * zeta) / 4.0;
            shape_grad[17-1](1) = (2 * xi - 1) * eta * (eta - 1) * (1 - zeta * zeta) / 4.0;
            shape_grad[17-1](2) = xi * (xi - 1) * (2 * eta - 1) * (1 - zeta * zeta) / 4.0;
            shape_grad[17-1](3) = -xi * (xi - 1) * eta * (eta - 1) * zeta / 2.0;

            shape_val[ 18-1] = xi * (xi + 1) * eta * (eta - 1) * (1 - zeta * zeta) / 4.0;
            shape_grad[18-1](1) = (2 * xi + 1) * eta * (eta - 1) * (1 - zeta * zeta) / 4.0;
            shape_grad[18-1](2) = xi * (xi + 1) * (2 * eta - 1) * (1 - zeta * zeta) / 4.0;
            shape_grad[18-1](3) = -xi * (xi + 1) * eta * (eta - 1) * zeta / 2.0;

            shape_val[ 19-1] = xi * (xi + 1) * eta * (eta + 1) * (1 - zeta * zeta) / 4.0;
            shape_grad[19-1](1) = (2 * xi + 1) * eta * (eta + 1) * (1 - zeta * zeta) / 4.0;
            shape_grad[19-1](2) = xi * (xi + 1) * (2 * eta + 1) * (1 - zeta * zeta) / 4.0;
            shape_grad[19-1](3) = -xi * (xi + 1) * eta * (eta + 1) * zeta / 2.0;

            shape_val[ 20-1] = xi * (xi - 1) * eta * (eta + 1) * (1 - zeta * zeta) / 4.0;
            shape_grad[20-1](1) = (2 * xi - 1) * eta * (eta + 1) * (1 - zeta * zeta) / 4.0;
            shape_grad[20-1](2) = xi * (xi - 1) * (2 * eta + 1) * (1 - zeta * zeta) / 4.0;
            shape_grad[20-1](3) = -xi * (xi - 1) * eta * (eta + 1) * zeta / 2.0;

            shape_val[ 21-1] = xi * (xi - 1) * (1 - eta * eta) * (1 - zeta * zeta) / 2.0;
            shape_grad[21-1](1) = (2 * xi - 1) * (1 - eta * eta) * (1 - zeta * zeta) / 2.0;
            shape_grad[21-1](2) = -xi * (xi - 1) * eta * (1 - zeta * zeta);
            shape_grad[21-1](3) = -xi * (xi - 1) * (1 - eta * eta) * zeta;

            shape_val[ 22-1] = xi * (xi + 1) * (1 - eta * eta) * (1 - zeta * zeta) / 2.0;
            shape_grad[22-1](1) = (2 * xi + 1) * (1 - eta * eta) * (1 - zeta * zeta) / 2.0;
            shape_grad[22-1](2) = -xi * (xi + 1) * eta * (1 - zeta * zeta);
            shape_grad[22-1](3) = -xi * (xi + 1) * (1 - eta * eta) * zeta;

            shape_val[ 23-1] = (1 - xi * xi) * eta * (eta - 1) * (1 - zeta * zeta) / 2.0;
            shape_grad[23-1](1) = -xi * eta * (eta - 1) * (1 - zeta * zeta);
            shape_grad[23-1](2) = (1 - xi * xi) * (2 * eta - 1) * (1 - zeta * zeta) / 2.0;
            shape_grad[23-1](3) = -(1 - xi * xi) * eta * (eta - 1) * zeta;

            shape_val[ 24-1] = (1 - xi * xi) * eta * (eta + 1) * (1 - zeta * zeta) / 2.0;
            shape_grad[24-1](1) = -xi * eta * (eta + 1) * (1 - zeta * zeta);
            shape_grad[24-1](2) = (1 - xi * xi) * (2 * eta + 1) * (1 - zeta * zeta) / 2.0;
            shape_grad[24-1](3) = -(1 - xi * xi) * eta * (eta + 1) * zeta;

            shape_val[ 25-1] = (1 - xi * xi) * (1 - eta * eta) * zeta * (zeta - 1) / 2.0;
            shape_grad[25-1](1) = -xi * (1 - eta * eta) * zeta * (zeta - 1);
            shape_grad[25-1](2) = -(1 - xi * xi) * eta * zeta * (zeta - 1);
            shape_grad[25-1](3) = (1 - xi * xi) * (1 - eta * eta) * (2 * zeta - 1) / 2.0;

            shape_val[ 26-1] = (1 - xi * xi) * (1 - eta * eta) * zeta * (zeta + 1) / 2.0;
            shape_grad[26-1](1) = -xi * (1 - eta * eta) * zeta * (zeta + 1);
            shape_grad[26-1](2) = -(1 - xi * xi) * eta * zeta * (zeta + 1);
            shape_grad[26-1](3) = (1 - xi * xi) * (1 - eta * eta) * (2 * zeta + 1) / 2.0;

            shape_val[ 27-1] = (1 - xi * xi) * (1 - eta * eta) * (1 - zeta * zeta);
            shape_grad[27-1](1) = -2 * xi * (1 - eta * eta) * (1 - zeta * zeta);
            shape_grad[27-1](2) = -2 * (1 - xi * xi) * eta * (1 - zeta * zeta);
            shape_grad[27-1](3) = -2 * (1 - xi * xi) * (1 - eta * eta) * zeta;

            break;
        }
        default:
            MessagePrinter::PrintErrorTxt("unsupported mesh type in Lagrange3DShapeFun calculation");
            MessagePrinter::AsFem_Exit();
            break;
    }

    _dxdxi  =0.0;_dydxi  =0.0;_dzdxi  =0.0;
    _dxdeta =0.0;_dydeta =0.0;_dzdeta =0.0;
    _dxdzeta=0.0;_dydzeta=0.0;_dzdzeta=0.0;
    for(int i=1;i<=_nNodes;i++){
        _dxdxi+=shape_grad[i-1](1)*nodes(i,1);
        _dydxi+=shape_grad[i-1](1)*nodes(i,2);
        _dzdxi+=shape_grad[i-1](1)*nodes(i,3);

        _dxdeta+=shape_grad[i-1](2)*nodes(i,1);
        _dydeta+=shape_grad[i-1](2)*nodes(i,2);
        _dzdeta+=shape_grad[i-1](2)*nodes(i,3);

        _dxdzeta+=shape_grad[i-1](3)*nodes(i,1);
        _dydzeta+=shape_grad[i-1](3)*nodes(i,2);
        _dzdzeta+=shape_grad[i-1](3)*nodes(i,3);
    }

    _Jac(1,1)=  _dxdxi;_Jac(1,2)=  _dydxi;_Jac(1,3)=  _dzdxi;
    _Jac(2,1)= _dxdeta;_Jac(2,2)= _dydeta;_Jac(2,3)= _dzdeta;
    _Jac(3,1)=_dxdzeta;_Jac(3,2)=_dydzeta;_Jac(3,3)=_dzdzeta;
    
    detjac=_Jac.Det();
    
    if(abs(detjac)<1.0e-15){
        MessagePrinter::PrintErrorTxt("singular element in 3D case, this error occurs in your 3D shape function calculation");
        MessagePrinter::AsFem_Exit();
    }

    if(flag){
        _XJac=_Jac.Inverse();
        double temp1,temp2,temp3;
        for(int i=1;i<=_nNodes;i++){
            temp1 =shape_grad[i-1](1)*_XJac(1,1)
              +shape_grad[i-1](2)*_XJac(1,2)
              +shape_grad[i-1](3)*_XJac(1,3);
            temp2 =shape_grad[i-1](1)*_XJac(2,1)
              +shape_grad[i-1](2)*_XJac(2,2)
              +shape_grad[i-1](3)*_XJac(2,3);

            temp3 =shape_grad[i-1](1)*_XJac(3,1)
              +shape_grad[i-1](2)*_XJac(3,2)
              +shape_grad[i-1](3)*_XJac(3,3);
            
            shape_grad[i-1](1) = temp1;
            shape_grad[i-1](2) = temp2;
            shape_grad[i-1](3) = temp3;
        }
    }

}
