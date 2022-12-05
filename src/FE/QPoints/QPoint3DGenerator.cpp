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
//+++ Purpose: implement the 3d gauss point generator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/QPoint3DGenerator.h"

void QPoint3DGenerator::generateQPoints(const int &t_order,const MeshType &t_meshtype,int &t_ngp,vector<double> &t_qpoints){ 
    switch (t_meshtype)
    {
    case MeshType::HEX8:
    case MeshType::HEX20:
    case MeshType::HEX27:
    {
        // for regular mesh, we can use the simple tensor product
        int ngp,l;
        vector<double> qp1d;
        QPoint1DGenerator q1dgen;
        q1dgen.generateQPoints(t_order,t_meshtype,ngp,qp1d);
        t_ngp=ngp*ngp*ngp;
        t_qpoints.resize(t_ngp*(1+3),0.0);
        l=0;
        for(int k=0;k<ngp;k++){
            for(int j=0;j<ngp;j++){
                for(int i=0;i<ngp;i++){
                    t_qpoints[l*4+1]=qp1d[i*2+1];
                    t_qpoints[l*4+2]=qp1d[j*2+1];
                    t_qpoints[l*4+3]=qp1d[k*2+1];
                    t_qpoints[l*4+0]=qp1d[i*2+0]*qp1d[j*2+0]*qp1d[k*2+0];
                    l+=1;
                }
            }
        }
        qp1d.clear();
        if(l!=t_ngp){
            MessagePrinter::printErrorTxt("generated qpoints number is less than expected, error detected in QPoint3DGenerator.cpp");
            MessagePrinter::exitAsFem();
        }
        return;break;
    }
    case MeshType::TET4:
    case MeshType::TET10:
    {
        // taken from "The Finite Element Method" by Zienkiewicz & Taylor
        switch (t_order)
        {
        case 0:
        case 1:
        {
            t_ngp=1;
            t_qpoints.resize(t_ngp*4,0.0);

            t_qpoints[(1-1)*4+1] = 0.25;
            t_qpoints[(1-1)*4+2] = 0.25;
            t_qpoints[(1-1)*4+3] = 0.25;
            t_qpoints[(1-1)*4+0] = 1.0/6.0;

            break;
        }
        case 2:
        {
            t_ngp=4;
            t_qpoints.resize(t_ngp*4,0.0);
            const double a = (5.0+3.0*sqrt(5.0))/20.0;
            const double b = (5.0-sqrt(5.0))/20.0;

            t_qpoints[(1-1)*4+1] = a;
            t_qpoints[(1-1)*4+2] = b;
            t_qpoints[(1-1)*4+3] = b;
            t_qpoints[(1-1)*4+0] = 1.0/24.0;

            t_qpoints[(2-1)*4+1] = b;
            t_qpoints[(2-1)*4+2] = a;
            t_qpoints[(2-1)*4+3] = b;
            t_qpoints[(2-1)*4+0] = 1.0/24.0;

            t_qpoints[(3-1)*4+1] = b;
            t_qpoints[(3-1)*4+2] = b;
            t_qpoints[(3-1)*4+3] = a;
            t_qpoints[(3-1)*4+0] = 1.0/24.0;

            t_qpoints[(4-1)*4+1] = b;
            t_qpoints[(4-1)*4+2] = b;
            t_qpoints[(4-1)*4+3] = b;
            t_qpoints[(4-1)*4+0] = 1.0/24.0;

            break;
        }
        case 3:
        {
            t_ngp=5;
            t_qpoints.resize(t_ngp*4,0.0);

            t_qpoints[(1-1)*4+1] = 0.25;
            t_qpoints[(1-1)*4+2] = 0.25;
            t_qpoints[(1-1)*4+3] = 0.25;
            t_qpoints[(1-1)*4+0] = -2.0/15.0;

            t_qpoints[(2-1)*4+1] = 1.0/6.0;
            t_qpoints[(2-1)*4+2] = 1.0/6.0;
            t_qpoints[(2-1)*4+3] = 1.0/6.0;
            t_qpoints[(2-1)*4+0] = 3.0/40.0;

            t_qpoints[(3-1)*4+1] = 1.0/6.0;
            t_qpoints[(3-1)*4+2] = 1.0/6.0;
            t_qpoints[(3-1)*4+3] = 0.5;
            t_qpoints[(3-1)*4+0] = 3.0/40.0;

            t_qpoints[(4-1)*4+1] = 1.0/6.0;
            t_qpoints[(4-1)*4+2] = 0.5;
            t_qpoints[(4-1)*4+3] = 1.0/6.0;
            t_qpoints[(4-1)*4+0] = 3.0/40.0;

            t_qpoints[(5-1)*4+1] = 0.5;
            t_qpoints[(5-1)*4+2] = 1.0/6.0;
            t_qpoints[(5-1)*4+3] = 1.0/6.0;
            t_qpoints[(5-1)*4+0] = 3.0/40.0;

            break;
        }
        case 4:
        case 5:
        {
            t_ngp=15;
            t_qpoints.resize(t_ngp*4,0.0);

            const double a=0.25;
            const double d=(5.0-sqrt(15.0))/20.0;
            const double e=(5.0+sqrt(15.0))/20.0;
            const double b1=(7.0+sqrt(15.0))/34.0;
            const double b2=(7.0-sqrt(15.0))/34.0;
            const double c1=(13.0-3.0*sqrt(15.0))/34.0;
            const double c2=(13.0+3.0*sqrt(15.0))/34.0;
            const double w1=8.0/405.0;
            const double w2=(2665.0-14.0*sqrt(15.0))/226800.0;
            const double w3=(2665.0+14.0*sqrt(15.0))/226800.0;
            const double w4=5.0/567.0;

            t_qpoints[(1-1)*4+1] = a;
            t_qpoints[(1-1)*4+2] = a;
            t_qpoints[(1-1)*4+3] = a;
            t_qpoints[(1-1)*4+0] = w1;
            //********************************
            t_qpoints[(2-1)*4+1] = b1;
            t_qpoints[(2-1)*4+2] = b1;
            t_qpoints[(2-1)*4+3] = b1;
            t_qpoints[(2-1)*4+0] = w2;

            t_qpoints[(3-1)*4+1] = b1;
            t_qpoints[(3-1)*4+2] = b1;
            t_qpoints[(3-1)*4+3] = c1;
            t_qpoints[(3-1)*4+0] = w2;

            t_qpoints[(4-1)*4+1] = b1;
            t_qpoints[(4-1)*4+2] = c1;
            t_qpoints[(4-1)*4+3] = b1;
            t_qpoints[(4-1)*4+0] = w2;

            t_qpoints[(5-1)*4+1] = c1;
            t_qpoints[(5-1)*4+2] = b1;
            t_qpoints[(5-1)*4+3] = b1;
            t_qpoints[(5-1)*4+0] = w2;
            //******************************
            t_qpoints[(6-1)*4+1] = b2;
            t_qpoints[(6-1)*4+2] = b2;
            t_qpoints[(6-1)*4+3] = b2;
            t_qpoints[(6-1)*4+0] = w3;

            t_qpoints[(7-1)*4+1] = b2;
            t_qpoints[(7-1)*4+2] = b2;
            t_qpoints[(7-1)*4+3] = c2;
            t_qpoints[(7-1)*4+0] = w3;

            t_qpoints[(8-1)*4+1] = b2;
            t_qpoints[(8-1)*4+2] = c2;
            t_qpoints[(8-1)*4+3] = b2;
            t_qpoints[(8-1)*4+0] = w3;

            t_qpoints[(9-1)*4+1] = c2;
            t_qpoints[(9-1)*4+2] = b2;
            t_qpoints[(9-1)*4+3] = b2;
            t_qpoints[(9-1)*4+0] = w3;
            //************************************
            t_qpoints[(10-1)*4+1] = d;
            t_qpoints[(10-1)*4+2] = d;
            t_qpoints[(10-1)*4+3] = e;
            t_qpoints[(10-1)*4+0] = w4;

            t_qpoints[(11-1)*4+1] = d;
            t_qpoints[(11-1)*4+2] = e;
            t_qpoints[(11-1)*4+3] = d;
            t_qpoints[(11-1)*4+0] = w4;

            t_qpoints[(12-1)*4+1] = e;
            t_qpoints[(12-1)*4+2] = d;
            t_qpoints[(12-1)*4+3] = d;
            t_qpoints[(12-1)*4+0] = w4;

            t_qpoints[(13-1)*4+1] = d;
            t_qpoints[(13-1)*4+2] = e;
            t_qpoints[(13-1)*4+3] = e;
            t_qpoints[(13-1)*4+0] = w4;

            t_qpoints[(14-1)*4+1] = e;
            t_qpoints[(14-1)*4+2] = d;
            t_qpoints[(14-1)*4+3] = e;
            t_qpoints[(14-1)*4+0] = w4;

            t_qpoints[(15-1)*4+1] = e;
            t_qpoints[(15-1)*4+2] = e;
            t_qpoints[(15-1)*4+3] = d;
            t_qpoints[(15-1)*4+0] = w4;

            break;
        }
        default:
        {
            MessagePrinter::printErrorTxt("order="+to_string(t_order)+" is unsupported in QPoint3DGenerator for tetrahedral mesh");
            MessagePrinter::exitAsFem();
            break;
        }
        }// end-of-order-switch
        return;
    }
    default:
    {
        MessagePrinter::printErrorTxt("unsupported mesh type in QPoint3DGenerator");
        MessagePrinter::exitAsFem();
        break;
    }
    return;
    }
}