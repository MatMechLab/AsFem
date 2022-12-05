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
//+++ Purpose: 3d hex20 shape function
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/ShapeFun3DHex20.h"

ShapeFun3DHex20::ShapeFun3DHex20(){

}

void ShapeFun3DHex20::calc3DShapeValsAndDerivatives(const double &xi,const double &eta,const double &zeta,
                                                    vector<double> &t_shpvals,
                                                    vector<Vector3d> &t_shpders){
    if(t_shpvals.size()<20 || t_shpders.size()<20){
        MessagePrinter::printErrorTxt("your shape val or derivs vector size is smaller than 20, error detected in ShapeFun3DHex20.cpp");
        MessagePrinter::exitAsFem();
    }

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
        t_shpvals[i-1]=(1.0+XI[i]*xi)*(1.0+ETA[i]*eta)*(1.0+ZETA[i]*zeta)*(XI[i]*xi+ETA[i]*eta+ZETA[i]*zeta-2.0)/8.0;

        t_shpders[i-1](1)=(XI[i])*(1.0+ETA[i]*eta)*(1.0+ZETA[i]*zeta)*(XI[i]*xi+ETA[i]*eta+ZETA[i]*zeta-2.0)/8.0
                        +(1.0+XI[i]*xi)*(1.0+ETA[i]*eta)*(1.0+ZETA[i]*zeta)*(XI[i])/8.0;

        t_shpders[i-1](2)=(1.0+XI[i]*xi)*(ETA[i])*(1.0+ZETA[i]*zeta)*(XI[i]*xi+ETA[i]*eta+ZETA[i]*zeta-2.0)/8.0
                        +(1.0+XI[i]*xi)*(1.0+ETA[i]*eta)*(1.0+ZETA[i]*zeta)*(ETA[i])/8.0;

        t_shpders[i-1](3)=(1.0+XI[i]*xi)*(1.0+ETA[i]*eta)*(ZETA[i])*(XI[i]*xi+ETA[i]*eta+ZETA[i]*zeta-2.0)/8.0
                        +(1.0+XI[i]*xi)*(1.0+ETA[i]*eta)*(1.0+ZETA[i]*zeta)*(ZETA[i])/8.0;
    }

    // for midside nodes
    for(i=1;i<=4;++i){
        // for 9,11,13,15
        t_shpvals[8+2*i-1-1]=(1.0-xi*xi)*(1.0+ETA[8+2*i-1]*eta)*(1.0+ZETA[8+2*i-1]*zeta)/4.0;
        t_shpders[8+2*i-1-1](1)=(-2.0*xi)*(1.0+ETA[8+2*i-1]*eta)*(1.0+ZETA[8+2*i-1]*zeta)/4.0;
        t_shpders[8+2*i-1-1](2)=(1.0-xi*xi)*(ETA[8+2*i-1])*(1.0+ZETA[8+2*i-1]*zeta)/4.0;
        t_shpders[8+2*i-1-1](3)=(1.0-xi*xi)*(1.0+ETA[8+2*i-1]*eta)*(ZETA[8+2*i-1])/4.0;
 
        // for 10,12,14,16
        t_shpvals[8+2*i-1]=(1.0-eta*eta)*(1.0+XI[8+2*i]*xi)*(1.0+ZETA[8+2*i]*zeta)/4.0;
        t_shpders[8+2*i-1](1)=(1.0-eta*eta)*(XI[8+2*i])*(1.0+ZETA[8+2*i]*zeta)/4.0;
        t_shpders[8+2*i-1](2)=(-2.0*eta)*(1.0+XI[8+2*i]*xi)*(1.0+ZETA[8+2*i]*zeta)/4.0;
        t_shpders[8+2*i-1](3)=(1.0-eta*eta)*(1.0+XI[8+2*i]*xi)*(ZETA[8+2*i])/4.0;

        // for 17,18,19,20
        t_shpvals[16+i-1]=(1.0-zeta*zeta)*(1.0+XI[16+i]*xi)*(1.0+ETA[16+i]*eta)/4.0;
        t_shpders[16+i-1](1)=(1.0-zeta*zeta)*(XI[16+i])*(1.0+ETA[16+i]*eta)/4.0;
        t_shpders[16+i-1](2)=(1.0-zeta*zeta)*(1.0+XI[16+i]*xi)*(ETA[16+i])/4.0;
        t_shpders[16+i-1](3)=(-2.0*zeta)*(1.0+XI[16+i]*xi)*(1.0+ETA[16+i]*eta)/4.0;
    }

}