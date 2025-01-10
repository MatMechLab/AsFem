//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2022.06.04
//+++ Purpose: implement the 2d gauss point generator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/QPoint2DGenerator.h"

void QPoint2DGenerator::generateQPoints(const int &t_order,const MeshType &t_meshtype,int &t_ngp,vector<double> &t_qpoints){
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    switch (t_meshtype)
    {
    case MeshType::QUAD4:
    case MeshType::QUAD8:
    case MeshType::QUAD9:
    {
        // for regular mesh, we can use the simple tensor product
        int ngp,k;
        vector<double> qp1d;
        QPoint1DGenerator q1dgen;
        q1dgen.generateQPoints(t_order,t_meshtype,ngp,qp1d);
        t_ngp=ngp*ngp;
        t_qpoints.resize(t_ngp*(1+2),0.0);
        k=0;
        for(int j=0;j<ngp;j++){
            for(int i=0;i<ngp;i++){
                t_qpoints[k*3+1]=qp1d[i*2+1];
                t_qpoints[k*3+2]=qp1d[j*2+1];
                t_qpoints[k*3+0]=qp1d[i*2+0]*qp1d[j*2+0];
                k+=1;
            }
        }
        qp1d.clear();
        if(k!=t_ngp){
            MessagePrinter::printErrorTxt("generated qpoints number is less than expected, error detected in QPoint2DGenerator.cpp");
            MessagePrinter::exitAsFem();
        }
        return;
    }
    case MeshType::TRI3:
    case MeshType::TRI6:
    {
        // for the details, one is referred to:
        // http://www.ce.memphis.edu/7111/notes/class_notes/chapter_03d_slides.pdf
        // Thanks for Xiaoyuan's comments, now, the weights are already divided by 2, the correct one is 0.5.
        switch (t_order)
        {
        case 0:
        case 1:
        {
            t_ngp=1;
            t_qpoints.resize(t_ngp*3,0.0);

            t_qpoints[(1-1)*3+1]= 1.0/3.0;
            t_qpoints[(1-1)*3+2]= 1.0/3.0;
            t_qpoints[(1-1)*3+0]= 0.5;
            break;
        }
        case 2:
        {
            t_ngp=3;
            t_qpoints.resize(t_ngp*3,0.0);

            t_qpoints[(1-1)*3+1]= 1.0/6.0;
            t_qpoints[(1-1)*3+2]= 1.0/6.0;
            t_qpoints[(1-1)*3+0]= 1.0/6.0;

            t_qpoints[(2-1)*3+1]= 2.0/3.0;
            t_qpoints[(2-1)*3+2]= 1.0/6.0;
            t_qpoints[(2-1)*3+0]= 1.0/6.0;

            t_qpoints[(3-1)*3+1]= 1.0/6.0;
            t_qpoints[(3-1)*3+2]= 2.0/3.0;
            t_qpoints[(3-1)*3+0]= 1.0/6.0;

            break;
        }
        case 3:
        {
            t_ngp=4;
            t_qpoints.resize(t_ngp*3,0.0);

            t_qpoints[(1-1)*3+1] = 1.550510257216821e-01;
            t_qpoints[(1-1)*3+2] = 1.785587282636164e-01;
            t_qpoints[(1-1)*3+0] = 1.590206908719885e-01;

            t_qpoints[(2-1)*3+1] = 6.449489742783178e-01;
            t_qpoints[(2-1)*3+2] = 7.503111022260811e-02;
            t_qpoints[(2-1)*3+0] = 9.097930912801141e-02;

            t_qpoints[(3-1)*3+1] = 1.550510257216821e-01;
            t_qpoints[(3-1)*3+2] = 6.663902460147013e-01;
            t_qpoints[(3-1)*3+0] = 1.590206908719885e-01;

            t_qpoints[(4-1)*3+1] = 6.449489742783178e-01;
            t_qpoints[(4-1)*3+2] = 2.800199154990740e-01;
            t_qpoints[(4-1)*3+0] = 9.097930912801141e-02;

            break;
        }
        case 4:
        {
            t_ngp=6;
            t_qpoints.resize(t_ngp*3,0.0);

            t_qpoints[(1-1)*3+1] = 0.0915762135;
            t_qpoints[(1-1)*3+2] = 0.8168475730;
            t_qpoints[(1-1)*3+0] = 0.1099517437*0.5;

            t_qpoints[(2-1)*3+1] = 0.0915762135;
            t_qpoints[(2-1)*3+2] = 0.0915762135;
            t_qpoints[(2-1)*3+0] = 0.1099517437*0.5;

            t_qpoints[(3-1)*3+1] = 0.8168475730;
            t_qpoints[(3-1)*3+2] = 0.0915762135;
            t_qpoints[(3-1)*3+0] = 0.1099517437*0.5;

            t_qpoints[(4-1)*3+1] = 0.4459484909;
            t_qpoints[(4-1)*3+2] = 0.1081030182;
            t_qpoints[(4-1)*3+0] = 0.2233815897*0.5;

            t_qpoints[(5-1)*3+1] = 0.4459484909;
            t_qpoints[(5-1)*3+2] = 0.4459484909;
            t_qpoints[(5-1)*3+0] = 0.2233815897*0.5;

            t_qpoints[(6-1)*3+1] = 0.1081030182;
            t_qpoints[(6-1)*3+2] = 0.4459484909;
            t_qpoints[(6-1)*3+0] = 0.2233815897*0.5;

            break;
        }
        case 5:
        {
            t_ngp=7;
            t_qpoints.resize(t_ngp*3,0.0);

            t_qpoints[(1-1)*3+1] = 1.0 / 3.0;
            t_qpoints[(1-1)*3+2] = 1.0 / 3.0;
            t_qpoints[(1-1)*3+0] = 0.2250000000*0.5;

            t_qpoints[(2-1)*3+1] = 0.1012865073;
            t_qpoints[(2-1)*3+2] = 0.7974269854;
            t_qpoints[(2-1)*3+0] = 0.1259391805*0.5;

            t_qpoints[(3-1)*3+1] = 0.1012865073;
            t_qpoints[(3-1)*3+2] = 0.1012865073;
            t_qpoints[(3-1)*3+0] = 0.1259391805*0.5;

            t_qpoints[(4-1)*3+1] = 0.7974269854;
            t_qpoints[(4-1)*3+2] = 0.1012865073;
            t_qpoints[(4-1)*3+0] = 0.1259391805*0.5;

            t_qpoints[(5-1)*3+1] = 0.0597158718;
            t_qpoints[(5-1)*3+2] = 0.4701420641;
            t_qpoints[(5-1)*3+0] = 0.1323941528*0.5;

            t_qpoints[(6-1)*3+1] = 0.4701420641;
            t_qpoints[(6-1)*3+2] = 0.4701420641;
            t_qpoints[(6-1)*3+0] = 0.1323941528*0.5;

            t_qpoints[(7-1)*3+1] = 0.4701420641;
            t_qpoints[(7-1)*3+2] = 0.0597158718;
            t_qpoints[(7-1)*3+0] = 0.1323941528*0.5;

            break;
        }
        default:{
            MessagePrinter::printErrorTxt("order="+to_string(t_order)+" is unsupported in QPoint2DGenerator for triangle mesh");
            MessagePrinter::exitAsFem();
            break;
        }
        }// end-of-order-switch
        return;
    }
    default:
    {
        MessagePrinter::printErrorTxt("unsupported mesh type in QPoint2DGenerator");
        MessagePrinter::exitAsFem();
        break;
    }
    return;
    }
}