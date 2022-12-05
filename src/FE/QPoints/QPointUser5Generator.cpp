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
//+++ Purpose: implement the user5 gauss point generator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/QPointUser5Generator.h"

void QPointUser5Generator::generateQPoints(const int &t_dim,const int &t_order,const MeshType &t_meshtype,int &t_ngp,vector<double> &t_qpoints){
    //************************************************************
    //*** here you should offer your own method for
    //*** the integration point generation
    //*** besides, you should consider about the cases
    //*** in different dimension, i.e., 2d and 3d case.
    //************************************************************

    //********************************
    //*** get rid of unused warning
    //********************************
    t_ngp=0;
    if(t_dim||t_order||t_meshtype==MeshType::NULLTYPE||t_qpoints.size()){}

    // switch (t_order)
    // {
    // case 0:
    // case 1:{
    //     t_ngp=1;
    //     t_qpoints.resize(t_ngp*2,0.0);
    //     t_qpoints[(1-1)*2+0]= 2.0;
    //     t_qpoints[(1-1)*2+1]= 0.0;
    //     break;
    //     }
    // default:{
    //     MessagePrinter::printErrorTxt("order="+to_string(t_order)+" is not supported in QPoint1DGenerator");
    //     MessagePrinter::exitAsFem();
    //     break;
    //     }
    // }

}