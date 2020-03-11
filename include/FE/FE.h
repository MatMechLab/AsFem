//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_FE_H
#define ASFEM_FE_H

#include <iostream>
#include <iomanip>
#include <string>

#include "petsc.h"

//**********************************
//*** AsFem's own header file    ***
//**********************************
#include "QPoint.h"
#include "ShapeFun.h"
#include "Mesh/Nodes.h"
#include "Mesh/Mesh.h"

class Mesh;

class FE
{
public:
    FE();
    void SetDim(int dim) {_nDim=dim;}
    void SetDimMin(int dim) {_nDimMin=dim;}
    void SetOrder(int order) {_nOrder=order;}
    void SetBCOrder(int order) {_nBCOrder=order;}
    void SetQPointType(string qtype) {_QPointType=qtype;}
    void InitFE(Mesh &mesh);// init the gauss point system and shape fun system
                            // before real calculcation
                            // 

    inline int GetDim() const {return _nDim;}
    inline int GetDimMin() const {return _nDimMin;}
    inline int GetOrder() const {return _nOrder;}
    inline int GetBCOrder() const {return _nBCOrder;}


    void PrintQPointInfo() const;


public:
    string _QPointType="gauss";
    bool _IsInit=false;
    int _nDim,_nDimMin,_nOrder,_nBCOrder;
    QPoint _qp_bulk,_qp_surface,_qp_line;
    ShapeFun _shp_bulk,_shp_surface,_shp_line;
    Nodes _nodes,_surface_nodes,_line_nodes;
    
};


#endif // ASFEM_FE_H