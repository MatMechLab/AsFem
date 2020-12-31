//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.07.12
//+++ Purpose: implement the general gauss integration class for
//+++          AsFem, here one can use:
//+++            1) standard Gauss-Legendre integration scheme
//+++            2) Gauss-Lobatto integration scheme
//+++            3) if you want to implement your own integration
//+++               scheme, all you need to do is create a new class
//+++               which should inherit from QPointBase, then override
//+++               the CreateQPoints() function
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/QPoint.h"

QPoint::QPoint()
:QPointGaussLegendre(),QPointGaussLobatto(){
    _CurrentQPType=QPointType::GAUSSLEGENDRE;
}
QPoint::QPoint(int dim,int order)
:QPointGaussLegendre(),QPointGaussLobatto(){
    QPointGaussLegendre::SetDim(dim);
    QPointGaussLegendre::SetQPointOrder(order);

    QPointGaussLobatto::SetDim(dim);
    QPointGaussLobatto::SetQPointOrder(order);
}
QPoint::QPoint(int dim,int order,QPointType qptype)
:QPointGaussLegendre(),QPointGaussLobatto(){
    QPointGaussLegendre::SetDim(dim);
    QPointGaussLegendre::SetQPointOrder(order);

    _CurrentQPType=qptype;

    QPointGaussLobatto::SetDim(dim);
    QPointGaussLobatto::SetQPointOrder(order);
}
//*******************************************************
void QPoint::CreateQpoints(MeshType meshtype){
    if(_CurrentQPType==QPointType::GAUSSLEGENDRE){
        QPointGaussLegendre::CreateQpoints(meshtype);
    }
    else if(_CurrentQPType==QPointType::GAUSSLOBATTO){
        QPointGaussLobatto::CreateQpoints(meshtype);
    }
}
void QPoint::Init(){
    if(_CurrentQPType==QPointType::GAUSSLEGENDRE){
        QPointGaussLegendre::Init();
    }
    else if(_CurrentQPType==QPointType::GAUSSLOBATTO){
        QPointGaussLobatto::Init();
    }
    else{
        MessagePrinter::PrintErrorTxt("unsupported gauss point type in QPoint initializing, please check your input file or your code");
        MessagePrinter::AsFem_Exit();
    }
}