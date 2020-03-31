//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "FE/ShapeFun.h"

void ShapeFun::Compute3DLagrangeShapeFun(const double &xi,
                                         const double &eta,
                                         const double &zeta,
                                         const Nodes &nodes,
                                         bool flag){
    switch (GetMeshType()){
        case MeshType::HEX8:
        {
            (*this)(1,0) = (1 - xi) * (1 - eta) * (1 - zeta) / 8.0;
            (*this)(1,1) = -(1 - eta) * (1 - zeta) / 8.0;
            (*this)(1,2) = -(1 - xi) * (1 - zeta) / 8.0;
            (*this)(1,3) = -(1 - xi) * (1 - eta) / 8.0;

            (*this)(2,0) = (1 + xi) * (1 - eta) * (1 - zeta) / 8.0;
            (*this)(2,1) = (1 - eta) * (1 - zeta) / 8.0;
            (*this)(2,2) = -(1 + xi) * (1 - zeta) / 8.0;
            (*this)(2,3) = -(1 + xi) * (1 - eta) / 8.0;

            (*this)(3,0) = (1 + xi) * (1 + eta) * (1 - zeta) / 8.0;
            (*this)(3,1) = (1 + eta) * (1 - zeta) / 8.0;
            (*this)(3,2) = (1 + xi) * (1 - zeta) / 8.0;
            (*this)(3,3) = -(1 + xi) * (1 + eta) / 8.0;

            (*this)(4,0) = (1 - xi) * (1 + eta) * (1 - zeta) / 8.0;
            (*this)(4,1) = -(1 + eta) * (1 - zeta) / 8.0;
            (*this)(4,2) = (1 - xi) * (1 - zeta) / 8.0;
            (*this)(4,3) = -(1 - xi) * (1 + eta) / 8.0;

            (*this)(5,0) = (1 - xi) * (1 - eta) * (1 + zeta) / 8.0;
            (*this)(5,1) = -(1 - eta) * (1 + zeta) / 8.0;
            (*this)(5,2) = -(1 - xi) * (1 + zeta) / 8.0;
            (*this)(5,3) = (1 - xi) * (1 - eta) / 8.0;

            (*this)(6,0) = (1 + xi) * (1 - eta) * (1 + zeta) / 8.0;
            (*this)(6,1) = (1 - eta) * (1 + zeta) / 8.0;
            (*this)(6,2) = -(1 + xi) * (1 + zeta) / 8.0;
            (*this)(6,3) = (1 + xi) * (1 - eta) / 8.0;

            (*this)(7,0) = (1 + xi) * (1 + eta) * (1 + zeta) / 8.0;
            (*this)(7,1) = (1 + eta) * (1 + zeta) / 8.0;
            (*this)(7,2) = (1 + xi) * (1 + zeta) / 8.0;
            (*this)(7,3) = (1 + xi) * (1 + eta) / 8.0;

            (*this)(8,0) = (1 - xi) * (1 + eta) * (1 + zeta) / 8.0;
            (*this)(8,1) = -(1 + eta) * (1 + zeta) / 8.0;
            (*this)(8,2) = (1 - xi) * (1 + zeta) / 8.0;
            (*this)(8,3) = (1 - xi) * (1 + eta) / 8.0;
            break;
        }
        case MeshType::HEX20:
        {
            /* this seems dosen't work
            const double XI[] = {0.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0};
            const double ETA[] = {0.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0};
            const double ZETA[] = {0.0,-1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,1.0};
            const double M1[] = {0.0,-1.0,1.0,-1.0,1.0};
            const double M2[] = {0.0,-1.0,-1.0,1.0,1.0};
            double xi0,eta0,zeta0,xi1,eta1,zeta1;
            int i;
            // for 8-corner node
            for (i = 1; i <= 8; ++i){
                xi0 = XI[i] * xi;
                eta0 = ETA[i] * eta;
                zeta0 = ZETA[i] * zeta;
                xi1 = 0.5 + 0.5 * xi0;
                eta1 = 0.5 + 0.5 * eta0;
                zeta1 = 0.5 + 0.5 * zeta0;
                (*this)(i,0) = xi1*eta1*zeta1*(xi0+eta0+zeta0-2.0);
                (*this)(i,1) = eta1*zeta1*(xi+0.5*XI[i]*(eta0+zeta0-1.0));
                (*this)(i,2) = xi1*zeta1*(eta+0.5*ETA[i]*(zeta0+xi0-1.0));
                (*this)(i,3) = xi1*eta1*(zeta+0.5*ZETA[i]*(xi0+eta0-1.0));
            }
            // for middle point
            const int ix[] = {0, 9, 11, 13, 15};
            const int iy[] = {0, 12, 10, 16, 14};
            const int iz[] = {0, 17, 18, 20, 19};
            for (i=1;i<=4;i++){
                xi1   = (1.0-xi*xi)*0.25;
                eta1  = 1.0+M1[i]*eta;
                zeta1 = 1.0+M2[i]*zeta;
                (*this)(ix[i],0)= xi1*eta1*zeta1;
                (*this)(ix[i],1)=-0.5*xi*eta1*zeta1;
                (*this)(ix[i],2)= xi1*M1[i]*zeta1;
                (*this)(ix[i],3)= xi1*eta1*M2[i];

                xi1  = 1.0+M1[i]*xi;
                eta1 = (1.0-eta*eta)*0.25;
                (*this)(iy[i],0)= xi1*eta1*zeta1;
                (*this)(iy[i],1)= M1[i]*eta1*zeta1;
                (*this)(iy[i],2)=-0.5*xi1*eta*zeta1;
                (*this)(iy[i],3)= xi1*eta1*M2[i];

                eta1  = 1.0+M2[i]*eta;
                zeta1 = (1.0-zeta*zeta)*0.25;
                (*this)(iz[i],0)= xi1*eta1*zeta1;
                (*this)(iz[i],1)= M1[i]*eta1*zeta1;
                (*this)(iz[i],2)= xi1*M2[i]*zeta1;
                (*this)(iz[i],3)=-0.5*xi1*eta1*zeta;
            }
            */
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
                (*this)(i,0)=(1.0+XI[i]*xi)*(1.0+ETA[i]*eta)*(1.0+ZETA[i]*zeta)*(XI[i]*xi+ETA[i]*eta+ZETA[i]*zeta-2.0)/8.0;
                
                (*this)(i,1)=(XI[i])*(1.0+ETA[i]*eta)*(1.0+ZETA[i]*zeta)*(XI[i]*xi+ETA[i]*eta+ZETA[i]*zeta-2.0)/8.0
                            +(1.0+XI[i]*xi)*(1.0+ETA[i]*eta)*(1.0+ZETA[i]*zeta)*(XI[i])/8.0;
                
                (*this)(i,2)=(1.0+XI[i]*xi)*(ETA[i])*(1.0+ZETA[i]*zeta)*(XI[i]*xi+ETA[i]*eta+ZETA[i]*zeta-2.0)/8.0
                            +(1.0+XI[i]*xi)*(1.0+ETA[i]*eta)*(1.0+ZETA[i]*zeta)*(ETA[i])/8.0;
                
                (*this)(i,3)=(1.0+XI[i]*xi)*(1.0+ETA[i]*eta)*(ZETA[i])*(XI[i]*xi+ETA[i]*eta+ZETA[i]*zeta-2.0)/8.0
                            +(1.0+XI[i]*xi)*(1.0+ETA[i]*eta)*(1.0+ZETA[i]*zeta)*(ZETA[i])/8.0;
            }

            // for midside nodes
            for(i=1;i<=4;++i){
                // for 9,11,13,15
                (*this)(8+2*i-1,0)=(1.0-xi*xi)*(1.0+ETA[8+2*i-1]*eta)*(1.0+ZETA[8+2*i-1]*zeta)/4.0;
                (*this)(8+2*i-1,1)=(-2.0*xi)*(1.0+ETA[8+2*i-1]*eta)*(1.0+ZETA[8+2*i-1]*zeta)/4.0;
                (*this)(8+2*i-1,2)=(1.0-xi*xi)*(ETA[8+2*i-1])*(1.0+ZETA[8+2*i-1]*zeta)/4.0;
                (*this)(8+2*i-1,3)=(1.0-xi*xi)*(1.0+ETA[8+2*i-1]*eta)*(ZETA[8+2*i-1])/4.0;

                // for 10,12,14,16
                (*this)(8+2*i,0)=(1.0-eta*eta)*(1.0+XI[8+2*i]*xi)*(1.0+ZETA[8+2*i]*zeta)/4.0;
                (*this)(8+2*i,1)=(1.0-eta*eta)*(XI[8+2*i])*(1.0+ZETA[8+2*i]*zeta)/4.0;
                (*this)(8+2*i,2)=(-2.0*eta)*(1.0+XI[8+2*i]*xi)*(1.0+ZETA[8+2*i]*zeta)/4.0;
                (*this)(8+2*i,3)=(1.0-eta*eta)*(1.0+XI[8+2*i]*xi)*(ZETA[8+2*i])/4.0;

                // for 17,18,19,20
                (*this)(16+i,0)=(1.0-zeta*zeta)*(1.0+XI[16+i]*xi)*(1.0+ETA[16+i]*eta)/4.0;
                (*this)(16+i,1)=(1.0-zeta*zeta)*(XI[16+i])*(1.0+ETA[16+i]*eta)/4.0;
                (*this)(16+i,2)=(1.0-zeta*zeta)*(1.0+XI[16+i]*xi)*(ETA[16+i])/4.0;
                (*this)(16+i,3)=(-2.0*zeta)*(1.0+XI[16+i]*xi)*(1.0+ETA[16+i]*eta)/4.0;
            }
            break;
        }
        case MeshType::HEX27:
        {
            (*this)(1,0) = xi * (xi - 1) * eta * (eta - 1) * zeta * (zeta - 1) / 8.0;
            (*this)(1,1) = (2 * xi - 1) * eta * (eta - 1) * zeta * (zeta - 1) / 8.0;
            (*this)(1,2) = xi * (xi - 1) * (2 * eta - 1) * zeta * (zeta - 1) / 8.0;
            (*this)(1,3) = xi * (xi - 1) * eta * (eta - 1) * (2 * zeta - 1) / 8.0;

            (*this)(2,0) = xi * (xi + 1) * eta * (eta - 1) * zeta * (zeta - 1) / 8.0;
            (*this)(2,1) = (2 * xi + 1) * eta * (eta - 1) * zeta * (zeta - 1) / 8.0;
            (*this)(2,2) = xi * (xi + 1) * (2 * eta - 1) * zeta * (zeta - 1) / 8.0;
            (*this)(2,3) = xi * (xi + 1) * eta * (eta - 1) * (2 * zeta - 1) / 8.0;

            (*this)(3,0) = xi * (xi + 1) * eta * (eta + 1) * zeta * (zeta - 1) / 8.0;
            (*this)(3,1) = (2 * xi + 1) * eta * (eta + 1) * zeta * (zeta - 1) / 8.0;
            (*this)(3,2) = xi * (xi + 1) * (2 * eta + 1) * zeta * (zeta - 1) / 8.0;
            (*this)(3,3) = xi * (xi + 1) * eta * (eta + 1) * (2 * zeta - 1) / 8.0;

            (*this)(4,0) = xi * (xi - 1) * eta * (eta + 1) * zeta * (zeta - 1) / 8.0;
            (*this)(4,1) = (2 * xi - 1) * eta * (eta + 1) * zeta * (zeta - 1) / 8.0;
            (*this)(4,2) = xi * (xi - 1) * (2 * eta + 1) * zeta * (zeta - 1) / 8.0;
            (*this)(4,3) = xi * (xi - 1) * eta * (eta + 1) * (2 * zeta - 1) / 8.0;

            (*this)(5,0) = xi * (xi - 1) * eta * (eta - 1) * zeta * (zeta + 1) / 8.0;
            (*this)(5,1) = (2 * xi - 1) * eta * (eta - 1) * zeta * (zeta + 1) / 8.0;
            (*this)(5,2) = xi * (xi - 1) * (2 * eta - 1) * zeta * (zeta + 1) / 8.0;
            (*this)(5,3) = xi * (xi - 1) * eta * (eta - 1) * (2 * zeta + 1) / 8.0;

            (*this)(6,0) = xi * (xi + 1) * eta * (eta - 1) * zeta * (zeta + 1) / 8.0;
            (*this)(6,1) = (2 * xi + 1) * eta * (eta - 1) * zeta * (zeta + 1) / 8.0;
            (*this)(6,2) = xi * (xi + 1) * (2 * eta - 1) * zeta * (zeta + 1) / 8.0;
            (*this)(6,3) = xi * (xi + 1) * eta * (eta - 1) * (2 * zeta + 1) / 8.0;

            (*this)(7,0) = xi * (xi + 1) * eta * (eta + 1) * zeta * (zeta + 1) / 8.0;
            (*this)(7,1) = (2 * xi + 1) * eta * (eta + 1) * zeta * (zeta + 1) / 8.0;
            (*this)(7,2) = xi * (xi + 1) * (2 * eta + 1) * zeta * (zeta + 1) / 8.0;
            (*this)(7,3) = xi * (xi + 1) * eta * (eta + 1) * (2 * zeta + 1) / 8.0;

            (*this)(8,0) = xi * (xi - 1) * eta * (eta + 1) * zeta * (zeta + 1) / 8.0;
            (*this)(8,1) = (2 * xi - 1) * eta * (eta + 1) * zeta * (zeta + 1) / 8.0;
            (*this)(8,2) = xi * (xi - 1) * (2 * eta + 1) * zeta * (zeta + 1) / 8.0;
            (*this)(8,3) = xi * (xi - 1) * eta * (eta + 1) * (2 * zeta + 1) / 8.0;

            (*this)(9,0) = (1 - xi * xi) * eta * (eta - 1) * zeta * (zeta - 1) / 4.0;
            (*this)(9,1) = -xi * eta * (eta - 1) * zeta * (zeta - 1) / 2.0;
            (*this)(9,2) = (1 - xi * xi) * (2 * eta - 1) * zeta * (zeta - 1) / 4.0;
            (*this)(9,3) = (1 - xi * xi) * eta * (eta - 1) * (2 * zeta - 1) / 4.0;

            (*this)(10,0) = xi * (xi + 1) * (1 - eta * eta) * zeta * (zeta - 1) / 4.0;
            (*this)(10,1) = (2 * xi + 1) * (1 - eta * eta) * zeta * (zeta - 1) / 4.0;
            (*this)(10,2) = -xi * (xi + 1) * eta * zeta * (zeta - 1) / 2.0;
            (*this)(10,3) = xi * (xi + 1) * (1 - eta * eta) * (2 * zeta - 1) / 4.0;

            (*this)(11,0) = (1 - xi * xi) * eta * (eta + 1) * zeta * (zeta - 1) / 4.0;
            (*this)(11,1) = -xi * eta * (eta + 1) * zeta * (zeta - 1) / 2.0;
            (*this)(11,2) = (1 - xi * xi) * (2 * eta + 1) * zeta * (zeta - 1) / 4.0;
            (*this)(11,3) = (1 - xi * xi) * eta * (eta + 1) * (2 * zeta - 1) / 4.0;

            (*this)(12,0) = xi * (xi - 1) * (1 - eta * eta) * zeta * (zeta - 1) / 4.0;
            (*this)(12,1) = (2 * xi - 1) * (1 - eta * eta) * zeta * (zeta - 1) / 4.0;
            (*this)(12,2) = -xi * (xi - 1) * eta * zeta * (zeta - 1) / 2.0;
            (*this)(12,3) = xi * (xi - 1) * (1 - eta * eta) * (2 * zeta - 1) / 4.0;

            (*this)(13,0) = (1 - xi * xi) * eta * (eta - 1) * zeta * (zeta + 1) / 4.0;
            (*this)(13,1) = -xi * eta * (eta - 1) * zeta * (zeta + 1) / 2.0;
            (*this)(13,2) = (1 - xi * xi) * (2 * eta - 1) * zeta * (zeta + 1) / 4.0;
            (*this)(13,3) = (1 - xi * xi) * eta * (eta - 1) * (2 * zeta + 1) / 4.0;

            (*this)(14,0) = xi * (xi + 1) * (1 - eta * eta) * zeta * (zeta + 1) / 4.0;
            (*this)(14,1) = (2 * xi + 1) * (1 - eta * eta) * zeta * (zeta + 1) / 4.0;
            (*this)(14,2) =-xi * (xi + 1) * eta * zeta * (zeta + 1) / 2.0;
            (*this)(14,3) = xi * (xi + 1) * (1 - eta * eta) * (2 * zeta + 1) / 4.0;

            (*this)(15,0) = (1 - xi * xi) * eta * (eta + 1) * zeta * (zeta + 1) / 4.0;
            (*this)(15,1) = -xi * eta * (eta + 1) * zeta * (zeta + 1) / 2.0;
            (*this)(15,2) = (1 - xi * xi) * (2 * eta + 1) * zeta * (zeta + 1) / 4.0;
            (*this)(15,3) = (1 - xi * xi) * eta * (eta + 1) * (2 * zeta + 1) / 4.0;

            (*this)(16,0) = xi * (xi - 1) * (1 - eta * eta) * zeta * (zeta + 1) / 4.0;
            (*this)(16,1) = (2 * xi - 1) * (1 - eta * eta) * zeta * (zeta + 1) / 4.0;
            (*this)(16,2) = -xi * (xi - 1) * eta * zeta * (zeta + 1) / 2.0;
            (*this)(16,3) = xi * (xi - 1) * (1 - eta * eta) * (2 * zeta + 1) / 4.0;

            (*this)(17,0) = xi * (xi - 1) * eta * (eta - 1) * (1 - zeta * zeta) / 4.0;
            (*this)(17,1) = (2 * xi - 1) * eta * (eta - 1) * (1 - zeta * zeta) / 4.0;
            (*this)(17,2) = xi * (xi - 1) * (2 * eta - 1) * (1 - zeta * zeta) / 4.0;
            (*this)(17,3) = -xi * (xi - 1) * eta * (eta - 1) * zeta / 2.0;

            (*this)(18,0) = xi * (xi + 1) * eta * (eta - 1) * (1 - zeta * zeta) / 4.0;
            (*this)(18,1) = (2 * xi + 1) * eta * (eta - 1) * (1 - zeta * zeta) / 4.0;
            (*this)(18,2) = xi * (xi + 1) * (2 * eta - 1) * (1 - zeta * zeta) / 4.0;
            (*this)(18,3) = -xi * (xi + 1) * eta * (eta - 1) * zeta / 2.0;

            (*this)(19,0) = xi * (xi + 1) * eta * (eta + 1) * (1 - zeta * zeta) / 4.0;
            (*this)(19,1) = (2 * xi + 1) * eta * (eta + 1) * (1 - zeta * zeta) / 4.0;
            (*this)(19,2) = xi * (xi + 1) * (2 * eta + 1) * (1 - zeta * zeta) / 4.0;
            (*this)(19,3) = -xi * (xi + 1) * eta * (eta + 1) * zeta / 2.0;

            (*this)(20,0) = xi * (xi - 1) * eta * (eta + 1) * (1 - zeta * zeta) / 4.0;
            (*this)(20,1) = (2 * xi - 1) * eta * (eta + 1) * (1 - zeta * zeta) / 4.0;
            (*this)(20,2) = xi * (xi - 1) * (2 * eta + 1) * (1 - zeta * zeta) / 4.0;
            (*this)(20,3) = -xi * (xi - 1) * eta * (eta + 1) * zeta / 2.0;

            (*this)(21,0) = xi * (xi - 1) * (1 - eta * eta) * (1 - zeta * zeta) / 2.0;
            (*this)(21,1) = (2 * xi - 1) * (1 - eta * eta) * (1 - zeta * zeta) / 2.0;
            (*this)(21,2) = -xi * (xi - 1) * eta * (1 - zeta * zeta);
            (*this)(21,3) = -xi * (xi - 1) * (1 - eta * eta) * zeta;

            (*this)(22,0) = xi * (xi + 1) * (1 - eta * eta) * (1 - zeta * zeta) / 2.0;
            (*this)(22,1) = (2 * xi + 1) * (1 - eta * eta) * (1 - zeta * zeta) / 2.0;
            (*this)(22,2) = -xi * (xi + 1) * eta * (1 - zeta * zeta);
            (*this)(22,3) = -xi * (xi + 1) * (1 - eta * eta) * zeta;

            (*this)(23,0) = (1 - xi * xi) * eta * (eta - 1) * (1 - zeta * zeta) / 2.0;
            (*this)(23,1) = -xi * eta * (eta - 1) * (1 - zeta * zeta);
            (*this)(23,2) = (1 - xi * xi) * (2 * eta - 1) * (1 - zeta * zeta) / 2.0;
            (*this)(23,3) = -(1 - xi * xi) * eta * (eta - 1) * zeta;

            (*this)(24,0) = (1 - xi * xi) * eta * (eta + 1) * (1 - zeta * zeta) / 2.0;
            (*this)(24,1) = -xi * eta * (eta + 1) * (1 - zeta * zeta);
            (*this)(24,2) = (1 - xi * xi) * (2 * eta + 1) * (1 - zeta * zeta) / 2.0;
            (*this)(24,3) = -(1 - xi * xi) * eta * (eta + 1) * zeta;

            (*this)(25,0) = (1 - xi * xi) * (1 - eta * eta) * zeta * (zeta - 1) / 2.0;
            (*this)(25,1) = -xi * (1 - eta * eta) * zeta * (zeta - 1);
            (*this)(25,2) = -(1 - xi * xi) * eta * zeta * (zeta - 1);
            (*this)(25,3) = (1 - xi * xi) * (1 - eta * eta) * (2 * zeta - 1) / 2.0;

            (*this)(26,0) = (1 - xi * xi) * (1 - eta * eta) * zeta * (zeta + 1) / 2.0;
            (*this)(26,1) = -xi * (1 - eta * eta) * zeta * (zeta + 1);
            (*this)(26,2) = -(1 - xi * xi) * eta * zeta * (zeta + 1);
            (*this)(26,3) = (1 - xi * xi) * (1 - eta * eta) * (2 * zeta + 1) / 2.0;

            (*this)(27,0) = (1 - xi * xi) * (1 - eta * eta) * (1 - zeta * zeta);
            (*this)(27,1) = -2 * xi * (1 - eta * eta) * (1 - zeta * zeta);
            (*this)(27,2) = -2 * (1 - xi * xi) * eta * (1 - zeta * zeta);
            (*this)(27,3) = -2 * (1 - xi * xi) * (1 - eta * eta) * zeta;
            

            break;
        }
        case MeshType::TET4:
        {
            const double sqrt2=sqrt(2.0);
            const double sqrt3=sqrt(3.0);

            (*this)(1,0)=(3.0+8.0*xi-2.0*sqrt2*zeta)/12.0;
            (*this)(1,1)=8.0/12.0;
            (*this)(1,2)=0.0;
            (*this)(1,3)=-2.0*sqrt2/12.0;

            (*this)(2,0)=(3.0-4.0*xi-4.0*sqrt3*eta-2.0*sqrt2*zeta)/12.0;
            (*this)(2,1)=-4.0/12.0;
            (*this)(2,2)=-4.0*sqrt3/12.0;
            (*this)(2,3)=-2.0*sqrt2/12.0;

            (*this)(3,0)=(3.0-4.0*xi+4.0*sqrt3*eta-2.0*sqrt2*zeta)/12.0;
            (*this)(3,1)=-4.0/12.0;
            (*this)(3,2)=4.0*sqrt3/12.0;
            (*this)(3,3)=-2.0*sqrt2/12.0;

            (*this)(4,0)=(1.0+2.0*sqrt2*zeta)/4.0;
            (*this)(4,1)=0.0;
            (*this)(4,2)=0.0;
            (*this)(4,3)=2.0*sqrt2/4.0;

            break;
        }
        default:
            break;
    }

    _dxdxi=0.0;_dydxi=0.0;_dzdxi=0.0;
    _dxdeta=0.0;_dydeta=0.0;_dzdeta=0.0;
    _dxdzeta=0.0;_dydzeta=0.0;_dzdzeta=0.0;
    for(int i=1;i<=GetShapeFunNums();i++){
        _dxdxi+=(*this)(i,1)*nodes(i,1);
        _dydxi+=(*this)(i,1)*nodes(i,2);
        _dzdxi+=(*this)(i,1)*nodes(i,3);

        _dxdeta+=(*this)(i,2)*nodes(i,1);
        _dydeta+=(*this)(i,2)*nodes(i,2);
        _dzdeta+=(*this)(i,2)*nodes(i,3);

        _dxdzeta+=(*this)(i,3)*nodes(i,1);
        _dydzeta+=(*this)(i,3)*nodes(i,2);
        _dzdzeta+=(*this)(i,3)*nodes(i,3);
    }
    
    _Jac[0][0]=  _dxdxi;_Jac[0][1]=  _dydxi;_Jac[0][2]=  _dzdxi;
    _Jac[1][0]= _dxdeta;_Jac[1][1]= _dydeta;_Jac[1][2]= _dzdeta;
    _Jac[2][0]=_dxdzeta;_Jac[2][1]=_dydzeta;_Jac[2][2]=_dzdzeta;

    // taken from https://en.wikipedia.org/wiki/Rule_of_Sarrus
    _DetJac=_Jac[0][0]*_Jac[1][1]*_Jac[2][2]
           +_Jac[0][1]*_Jac[1][2]*_Jac[2][0]
           +_Jac[0][2]*_Jac[1][0]*_Jac[2][1]
           -_Jac[2][0]*_Jac[1][1]*_Jac[0][2]
           -_Jac[2][1]*_Jac[1][2]*_Jac[0][0]
           -_Jac[2][2]*_Jac[1][0]*_Jac[0][1];

    if(abs(_DetJac)<1.0e-15){
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: singular element in 3D case                          !!!   ***\n");
        Msg_AsFem_Exit();
    }

    

    // taken from: http://mathworld.wolfram.com/MatrixInverse.html
    _XJac[0][0]=(_Jac[1][1]*_Jac[2][2]-_Jac[1][2]*_Jac[2][1])/_DetJac;
    _XJac[0][1]=(_Jac[0][2]*_Jac[2][1]-_Jac[0][1]*_Jac[2][2])/_DetJac;
    _XJac[0][2]=(_Jac[0][1]*_Jac[1][2]-_Jac[0][2]*_Jac[1][1])/_DetJac;

    _XJac[1][0]=(_Jac[1][2]*_Jac[2][0]-_Jac[1][0]*_Jac[2][2])/_DetJac;
    _XJac[1][1]=(_Jac[0][0]*_Jac[2][2]-_Jac[0][2]*_Jac[2][0])/_DetJac;
    _XJac[1][2]=(_Jac[0][2]*_Jac[1][0]-_Jac[0][0]*_Jac[1][2])/_DetJac;

    _XJac[2][0]=(_Jac[1][0]*_Jac[2][1]-_Jac[1][1]*_Jac[2][0])/_DetJac;
    _XJac[2][1]=(_Jac[0][1]*_Jac[2][0]-_Jac[0][0]*_Jac[2][1])/_DetJac;
    _XJac[2][2]=(_Jac[0][0]*_Jac[1][1]-_Jac[0][1]*_Jac[1][0])/_DetJac;

    double temp1,temp2,temp3;
    for(int i=1;i<=GetShapeFunNums();i++){

        if(flag){
            temp1 =(*this)(i,1)*_XJac[0][0]
              +(*this)(i,2)*_XJac[0][1]
              +(*this)(i,3)*_XJac[0][2];
            temp2 =(*this)(i,1)*_XJac[1][0]
              +(*this)(i,2)*_XJac[1][1]
              +(*this)(i,3)*_XJac[1][2];
              
            temp3 =(*this)(i,1)*_XJac[2][0]
              +(*this)(i,2)*_XJac[2][1]
              +(*this)(i,3)*_XJac[2][2];

            (*this)(i,1) = temp1;
            (*this)(i,2) = temp2;
            (*this)(i,3) = temp3;
        }

        _shape_value[i-1]=(*this)(i,0);

        _shape_grad[i-1](1)=(*this)(i,1);
        _shape_grad[i-1](2)=(*this)(i,2);
        _shape_grad[i-1](3)=(*this)(i,3);
    }
}