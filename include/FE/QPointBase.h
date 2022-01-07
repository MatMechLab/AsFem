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
//+++ Purpose: Define the base class for qpoint
//+++          the GaussÃÂ¢ÃÂÃÂLegendre rule, GaussÃÂ¢ÃÂÃÂLobatto rules as well
//+++          as other rules should inherit from this calss
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <string>
#include <cmath>


//**************************************
//*** for AsFem's own header file
//**************************************
#include "Utils/MessagePrinter.h"
#include "Mesh/MeshType.h"
#include "FE/QPointType.h"

class QPointBase{
public:
    QPointBase(){
        _nQpOrder=1;_nQpPoints=1;_nDim=1;
        _QpCoords.clear();
        _QpType=QPointType::GAUSSLEGENDRE;
    }
    virtual void CreateQpoints(MeshType meshtype)=0;


protected:
    QPointType _QpType;
    int _nQpOrder,_nQpPoints,_nDim;
    vector<double> _QpCoords;
};