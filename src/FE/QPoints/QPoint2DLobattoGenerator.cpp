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
//+++ Date   : 2022.06.05
//+++ Purpose: implement the 2d lobatto gauss point generator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/QPoint2DLobattoGenerator.h"

void QPoint2DLobattoGenerator::generateQPoints(const int &t_order,const MeshType &t_meshtype,int &t_ngp,vector<double> &t_qpoints){ 
    switch (t_meshtype)
    {
    case MeshType::QUAD4:
    case MeshType::QUAD8:
    case MeshType::QUAD9:
    {
        // for regular mesh, we can use the simple tensor product
        int ngp,k;
        vector<double> qp1d;
        QPoint1DLobattoGenerator q1dgen;
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
            MessagePrinter::printErrorTxt("generated qpoints number is less than expected, error detected in QPoint2DLobattoGenerator.cpp");
            MessagePrinter::exitAsFem();
        }
        return;
    }
    default:
    {
        MessagePrinter::printErrorTxt("unsupported mesh type in QPoint2DLobattoGenerator");
        MessagePrinter::exitAsFem();
        break;
    }
    return;
    }
}