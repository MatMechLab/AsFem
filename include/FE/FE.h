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
//+++ Purpose: implement the general FE space for FEM calculation
//+++          in AsFem, here one can use:
//+++            1) gauss integration 
//+++            2) shape functions for different mesh
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once


#include "Mesh/Mesh.h"
#include "FE/QPoint.h"

class Mesh;

class FE{
public:
    FE();

    void SetDim(int dim){_nDim=dim;_HasDimSet=true;}
    //***********************************************
    //*** for QPoint
    //***********************************************
    void SetQPointType(QPointType qptype);
    void SetBulkQpOrder(int order);
    void SetBCQpOrder(int order);
    void CreateQPoints(Mesh &mesh);


    //***********************************************
    //*** for get functions
    //***********************************************
    inline int GetDim()const{return _nDim;}

    QPoint& GetBulkQPointPtr(){return _BulkQPoint;}
    QPoint& GetLineQPointPtr(){return _LineQPoint;}
    QPoint& GetSurfaceQPointPtr(){return _SurfaceQPoint;}


    void PrintFEInfo()const;

private:
    int _nDim;
    bool _HasDimSet=false;
    QPoint _BulkQPoint,_LineQPoint,_SurfaceQPoint;
    int _nBulkQpOrder,_nBCQpOrder;
    
};