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
//+++ Purpose: implement the 3d lobatto gauss point generator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/QPoint3DLobattoGenerator.h"

void QPoint3DLobattoGenerator::generateQPoints(const int &t_order,const MeshType &t_meshtype,int &t_ngp,vector<double> &t_qpoints){ 
    switch (t_meshtype)
    {
    case MeshType::HEX8:
    case MeshType::HEX20:
    case MeshType::HEX27:
    {
        // for regular mesh, we can use the simple tensor product
        int ngp,l;
        vector<double> qp1d;
        QPoint1DLobattoGenerator q1dgen;
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
            MessagePrinter::printErrorTxt("generated qpoints number is less than expected, error detected in QPoint3DLobattoGenerator.cpp");
            MessagePrinter::exitAsFem();
        }
        return;break;
    }
    default:
    {
        MessagePrinter::printErrorTxt("unsupported mesh type in QPoint3DLobattoGenerator");
        MessagePrinter::exitAsFem();
        break;
    }
    return;
    }
}