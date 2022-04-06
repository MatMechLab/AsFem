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

/**
 * this class defines the abstract class for gauss quadrature point generation.
 * All the other gauss point generator should inherit this one!
 */
class QPointBase{
public:
    /**
     * constructor
     */
    QPointBase(){
        _nQpOrder=1;_nQpPoints=1;_nDim=1;
        _QpCoords.clear();
        _QpType=QPointType::GAUSSLEGENDRE;
    }

    /**
     * qpoint generation method, which should be implement by the child class
     * @param meshtype the type of the specific mesh
     */
    virtual void CreateQpoints(MeshType meshtype)=0;


protected:
    QPointType _QpType;
    int _nQpOrder,_nQpPoints,_nDim;
    vector<double> _QpCoords;
};