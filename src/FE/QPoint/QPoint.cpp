//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "FE/QPoint.h"

QPoint::QPoint(){
    _qp_coords.clear();
    _nQpPoints=0;
    _nOrder=0;
    _nDim=0;
    _QPointType="gauss";
    _HasSettings=false;
    _HasDim=false;
    _HasOrder=false;
}
//**********************+
QPoint::QPoint(int dim,int order){
    _qp_coords.clear();
    _nQpPoints=0;
    _QPointType="gauss";
    _HasSettings=false;
    _HasDim=false;
    _HasOrder=false;

    _nDim=dim;_nOrder=order;
}
