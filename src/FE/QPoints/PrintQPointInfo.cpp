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
//+++ Date   : 2020.07.12
//+++ Purpose: implement the printer for some basic information of
//+++          the QPoint class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/QPoint.h"

void QPoint::PrintQPointInfo()const{
    string msg;
    MessagePrinter::PrintNormalTxt("Summary information of qpoints");
    if(_CurrentQPType==QPointType::GAUSSLEGENDRE){
        msg="  qpoint type=Gauss-Legendre";
        MessagePrinter::PrintNormalTxt(msg);
    }
    else if(_CurrentQPType==QPointType::GAUSSLOBATTO){
        msg="  qpoint type=Gauss-Lobatto";
        MessagePrinter::PrintNormalTxt(msg);
    }
    msg="  dim="+to_string(GetDim())+", order="+to_string(GetQpOrder())+", npoints="+to_string(GetQpPointsNum());
    MessagePrinter::PrintNormalTxt(msg);
    MessagePrinter::PrintDashLine();
}
//***************************************************
void QPoint::PrintQPointDetailInfo()const{
    string msg;
    MessagePrinter::PrintNormalTxt("Summary information of qpoints");
    if(_CurrentQPType==QPointType::GAUSSLEGENDRE){
        msg="  qpoint type=Gauss-Legendre";
        MessagePrinter::PrintNormalTxt(msg);
    }
    else if(_CurrentQPType==QPointType::GAUSSLOBATTO){
        msg="  qpoint type=Gauss-Lobatto";
        MessagePrinter::PrintNormalTxt(msg);
    }
    msg="  dim="+to_string(GetDim())+", order="+to_string(GetQpOrder())+", npoints="+to_string(GetQpPointsNum());
    MessagePrinter::PrintNormalTxt(msg);
    for(int i=1;i<=GetQpPointsNum();i++){
        msg.clear();
        msg="  "+to_string(i)+"-th qpoint:";
        msg+="x="+to_string(GetIthQpPointJthCoord(i,1));
        if(GetDim()>=2){
            msg+=", y="+to_string(GetIthQpPointJthCoord(i,2));
        }
        if(GetDim()==3){
            msg+=", z="+to_string(GetIthQpPointJthCoord(i,3));
        }
        msg+=", w="+to_string(GetIthQpPointJthCoord(i,0));
        MessagePrinter::PrintNormalTxt(msg);
    }
    MessagePrinter::PrintDashLine();
}