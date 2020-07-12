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
//+++ Purpose: implement the Gauss–Lobatto rule for the integration
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/QPointGaussLobatto.h"

QPointGaussLobatto::QPointGaussLobatto()
:QPointBase(){
    _nQpOrder=1;_nQpPoints=1;_nDim=1;
    _QpCoords.clear();
    _QpType=QPointType::GAUSSLOBATTO;

    _HasSettings=true;
    _HasDim=false;
    _HasOrder=false;
}

QPointGaussLobatto::QPointGaussLobatto(int dim,int order)
:QPointBase(){
    _nQpOrder=order;_nDim=dim;
    _QpCoords.clear();
    _QpType=QPointType::GAUSSLOBATTO;

    _HasSettings=true;
    _HasDim=true;
    _HasOrder=true;
}