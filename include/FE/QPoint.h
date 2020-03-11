//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_QPOINT_H
#define ASFEM_QPOINT_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>

#include "petsc.h"

//******************************************
#include "MessagePrinter/MessagePrinter.h"
#include "Mesh/MeshType.h"

using namespace std;

class QPoint
{
public:
    QPoint();
    QPoint(int dim,int order);


    void SetDim(int dim) {_nDim=dim;}
    void SetQPointOrder(int order) {_nOrder=order;}
    void SetQPointType(string type="gauss") {_QPointType=type;}
    void CreateQPoints(MeshType meshtype);

    // get function
    inline int GetDim() const {return _nDim;}
    inline int GetQpOrder() const {return _nOrder;}
    inline int GetQpPointsNum() const {return _nQpPoints;}
    inline string GetQpType()const{return _QPointType;}
    
    inline double operator()(int i,int j) const {return _qp_coords[(i-1)*(_nDim+1)+j];}
    inline double& operator()(int i,int j) {return _qp_coords[(i-1)*(_nDim+1)+j];}

    inline double GetIthQpPointJthCoord(int i,int j) const {return _qp_coords[(i-1)*(_nDim+1)+j];}

private:
    void Create1DGaussPoint();
    void Create1DGaussLobattoPoint();

    void Create2DGaussPoint(MeshType meshtype);
    void Create2DGaussLobattoPoint(MeshType meshtype);

    void Create3DGaussPoint(MeshType meshtype);
    void Create3DGaussLobattoPoint(MeshType meshtype);

private:
    vector<double> _qp_coords;
    int _nQpPoints,_nOrder;
    int _nDim;
    string _QPointType;
    bool _HasSettings=false;
    bool _HasDim,_HasOrder;

};


#endif // ASFEM_QPOINT_H