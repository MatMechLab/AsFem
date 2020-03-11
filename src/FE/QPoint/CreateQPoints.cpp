//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "FE/QPoint.h"

void QPoint::CreateQPoints(MeshType meshtype){
    if(GetDim()==1){
        if(_QPointType=="gauss"){
            Create1DGaussPoint();
        }
        else if(_QPointType=="gausslobatto"){
            Create1DGaussLobattoPoint();
        }
    }
    else if(GetDim()==2){
        if(_QPointType=="gauss"){
            Create2DGaussPoint(meshtype);
        }
        else if(_QPointType=="gausslobatto"){
            Create2DGaussLobattoPoint(meshtype);
        }
    }
    else if(GetDim()==3){
        if(_QPointType=="gauss"){
            Create3DGaussPoint(meshtype);
        }
        else if(_QPointType=="gausslobatto"){
            Create3DGaussLobattoPoint(meshtype);
        }
    }
}