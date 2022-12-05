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
//+++ Date   : 2022.06.04
//+++ Purpose: implement the 1d gauss point generator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/QPoint1DGenerator.h"

void QPoint1DGenerator::generateQPoints(const int &t_order,const MeshType &t_meshtype,int &t_ngp,vector<double> &t_qpoints){
    if(t_meshtype==MeshType::EDGE2){}// for 1d case, we don't need to consider the mesh type
    // here the numbers are taken from:
    //    https://pomax.github.io/bezierinfo/legendre-gauss.html
    // here (orders+1)/2=nodes
    t_qpoints.clear();
    switch (t_order)
    {
    case 0:
    case 1:{
        t_ngp=1;
        t_qpoints.resize(t_ngp*2,0.0);
        t_qpoints[(1-1)*2+0]= 2.0;
        t_qpoints[(1-1)*2+1]= 0.0;
        break;
        }
    case 2:
    case 3:{
        t_ngp=2;
        t_qpoints.resize(t_ngp*2,0.0);
        t_qpoints[(1-1)*2+0]= 1.0;
        t_qpoints[(1-1)*2+1]=-sqrt(1.0/3.0);

        t_qpoints[(2-1)*2+0]= 1.0;
        t_qpoints[(2-1)*2+1]= sqrt(1.0/3.0);
        break;
        }
    case 4:
    case 5:{
        t_ngp=3;
        t_qpoints.resize(t_ngp*2,0.0);
        t_qpoints[(1-1)*2+0]= 5.0/9.0;
        t_qpoints[(1-1)*2+1]=-sqrt(0.6);

        t_qpoints[(2-1)*2+0]= 8.0/(9.0);
        t_qpoints[(2-1)*2+1]= 0.0;

        t_qpoints[(3-1)*2+0]= 5.0/(9.0);
        t_qpoints[(3-1)*2+1]= sqrt(0.6);
        break;
        }
    case 6:
    case 7:{
        t_ngp=4;
        t_qpoints.resize(t_ngp*2,0.0);
        const double t=sqrt(4.8);
        t_qpoints[(1-1)*2+1]=-sqrt((3.0+t)/7.0);
        t_qpoints[(2-1)*2+1]=-sqrt((3.0-t)/7.0);
        t_qpoints[(3-1)*2+1]= sqrt((3.0-t)/7.0);
        t_qpoints[(4-1)*2+1]= sqrt((3.0+t)/7.0);
        const double w=1.0/3.0/t;
        t_qpoints[(1-1)*2+0]=0.5-w;
        t_qpoints[(2-1)*2+0]=0.5+w;
        t_qpoints[(3-1)*2+0]=0.5+w;
        t_qpoints[(4-1)*2+0]=0.5-w;

        break;
        }
    case 8:
    case 9:{
        t_ngp=5;
        t_qpoints.resize(t_ngp*2,0.0);
        const double t = sqrt(1120.0);

        t_qpoints[(1-1)*2+1] =-sqrt((70.0+t)/126.0);
        t_qpoints[(2-1)*2+1] =-sqrt((70.0-t)/126.0);
        t_qpoints[(3-1)*2+1] = 0.000000000000000000;
        t_qpoints[(4-1)*2+1] = sqrt((70.0-t)/126.0);
        t_qpoints[(5-1)*2+1] = sqrt((70.0+t)/126.0);

        t_qpoints[(1-1)*2+0] = (21.0*t+117.60)/(t*(70.0+t));
        t_qpoints[(2-1)*2+0] = (21.0*t-117.60)/(t*(70.0-t));
        t_qpoints[(3-1)*2+0] = 2.0 *(1.0-t_qpoints[(1-1)*2+0]-t_qpoints[(2-1)*2+0]);
        t_qpoints[(4-1)*2+0] = t_qpoints[(2-1)*2+0];
        t_qpoints[(5-1)*2+0] = t_qpoints[(1-1)*2+0];
        
        break;
        }
    case 10:
    case 11:{
        t_ngp=6;
        t_qpoints.resize(t_ngp*2,0.0);

        t_qpoints[(1-1)*2+1] =-0.9324695142031521;
        t_qpoints[(2-1)*2+1] =-0.6612093864662645;
        t_qpoints[(3-1)*2+1] =-0.2386191860831969;
        t_qpoints[(4-1)*2+1] = 0.2386191860831969;
        t_qpoints[(5-1)*2+1] = 0.6612093864662645;
        t_qpoints[(6-1)*2+1] = 0.9324695142031521;

        t_qpoints[(1-1)*2+0] = 0.1713244923791704;
        t_qpoints[(2-1)*2+0] = 0.3607615730481386;
        t_qpoints[(3-1)*2+0] = 0.4679139345726910;
        t_qpoints[(4-1)*2+0] = 0.4679139345726910;
        t_qpoints[(5-1)*2+0] = 0.3607615730481386;
        t_qpoints[(6-1)*2+0] = 0.1713244923791704;
        
        break;
        }
    case 12:
    case 13:{
        t_ngp=7;
        t_qpoints.resize(t_ngp*2,0.0);

        t_qpoints[(1-1)*2+1] =-0.9491079123427585;
        t_qpoints[(2-1)*2+1] =-0.7415311855993945;
        t_qpoints[(3-1)*2+1] =-0.4058451513773972;
        t_qpoints[(4-1)*2+1] = 0.0000000000000000;
        t_qpoints[(5-1)*2+1] = 0.4058451513773972;
        t_qpoints[(6-1)*2+1] = 0.7415311855993945;
        t_qpoints[(7-1)*2+1] = 0.9491079123427585;

        t_qpoints[(1-1)*2+0] = 0.1294849661688697;
        t_qpoints[(2-1)*2+0] = 0.2797053914892766;
        t_qpoints[(3-1)*2+0] = 0.3818300505051189;
        t_qpoints[(4-1)*2+0] = 0.4179591836734694;
        t_qpoints[(5-1)*2+0] = 0.3818300505051189;
        t_qpoints[(6-1)*2+0] = 0.2797053914892766;
        t_qpoints[(7-1)*2+0] = 0.1294849661688697;
        
        break;
        }
    case 14:
    case 15:{
        t_ngp=8;
        t_qpoints.resize(t_ngp*2,0.0);

        t_qpoints[(1-1)*2+1] =-0.9602898564975363;
        t_qpoints[(2-1)*2+1] =-0.7966664774136267;
        t_qpoints[(3-1)*2+1] =-0.5255324099163290;
        t_qpoints[(4-1)*2+1] =-0.1834346424956498;
        t_qpoints[(5-1)*2+1] = 0.1834346424956498;
        t_qpoints[(6-1)*2+1] = 0.5255324099163290;
        t_qpoints[(7-1)*2+1] = 0.7966664774136267;
        t_qpoints[(8-1)*2+1] = 0.9602898564975363;

        t_qpoints[(1-1)*2+0] = 0.1012285362903763;
        t_qpoints[(2-1)*2+0] = 0.2223810344533745;
        t_qpoints[(3-1)*2+0] = 0.3137066458778873;
        t_qpoints[(4-1)*2+0] = 0.3626837833783620;
        t_qpoints[(5-1)*2+0] = 0.3626837833783620;
        t_qpoints[(6-1)*2+0] = 0.3137066458778873;
        t_qpoints[(7-1)*2+0] = 0.2223810344533745;
        t_qpoints[(8-1)*2+0] = 0.1012285362903763;
        
        break;
        }
    case 16:
    case 17:{
        t_ngp=9;
        t_qpoints.resize(t_ngp*2,0.0);

        t_qpoints[(1-1)*2+1] =-0.9681602395076261;
        t_qpoints[(2-1)*2+1] =-0.8360311073266358;
        t_qpoints[(3-1)*2+1] =-0.6133714327005904;
        t_qpoints[(4-1)*2+1] =-0.3242534234038089;
        t_qpoints[(5-1)*2+1] = 0.0000000000000000;
        t_qpoints[(6-1)*2+1] = 0.3242534234038089;
        t_qpoints[(7-1)*2+1] = 0.6133714327005904;
        t_qpoints[(8-1)*2+1] = 0.8360311073266358;
        t_qpoints[(9-1)*2+1] = 0.9681602395076261;

        t_qpoints[(1-1)*2+0] = 0.0812743883615744;
        t_qpoints[(2-1)*2+0] = 0.1806481606948574;
        t_qpoints[(3-1)*2+0] = 0.2606106964029354;
        t_qpoints[(4-1)*2+0] = 0.3123470770400029;
        t_qpoints[(5-1)*2+0] = 0.3302393550012598;
        t_qpoints[(6-1)*2+0] = 0.3123470770400029;
        t_qpoints[(7-1)*2+0] = 0.2606106964029354;
        t_qpoints[(8-1)*2+0] = 0.1806481606948574;
        t_qpoints[(9-1)*2+0] = 0.0812743883615744;
        
        break;
        }
    default:{
        MessagePrinter::printErrorTxt("order="+to_string(t_order)+" is not supported in QPoint1DGenerator");
        MessagePrinter::exitAsFem();
        break;
        }
    }
}