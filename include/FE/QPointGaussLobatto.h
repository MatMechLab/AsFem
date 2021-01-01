//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.07.12
//+++ Purpose: implement the GaussâLobatto rule for the integration
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "FE/QPointBase.h"

class QPointGaussLobatto:public QPointBase{
public:
    QPointGaussLobatto();
    QPointGaussLobatto(int dim,int order);
    //************************************************
    //*** for basic settings
    //************************************************
    void SetDim(int dim) {_nDim=dim;_HasDim=true;}
    void SetQPointOrder(int order) {_nQpOrder=order;_HasOrder=true;}
    void SetQPointType(QPointType qptype) {_QpType=qptype;_HasSettings=true;}
    void Init();
    //************************************************
    //*** for basic gettings
    //************************************************
    inline int GetDim() const {return _nDim;}
    inline int GetQpOrder() const {return _nQpOrder;}
    inline int GetQpPointsNum() const {return _nQpPoints;}
    inline QPointType GetQpPointType()const{return _QpType;}

    //************************************************
    //*** for operator overload
    //************************************************
    inline double operator()(int i,int j)const{return _QpCoords[(i-1)*(_nDim+1)+j];}
    inline double& operator()(int i,int j){return _QpCoords[(i-1)*(_nDim+1)+j];}

    inline double& GetIthQpPointJthCoord(int i,int j){return _QpCoords[(i-1)*(_nDim+1)+j];}
    inline double  GetIthQpPointJthCoord(int i,int j)const{return _QpCoords[(i-1)*(_nDim+1)+j];}

    virtual void CreateQpoints(MeshType meshtype) override;

private:
    void Create1DGaussPoint();
    void Create2DGaussPoint(MeshType meshtype);
    void Create3DGaussPoint(MeshType meshtype);


    //***********************************************
    //*** for some basic variables
    //***********************************************
    bool _HasSettings=false;
    bool _HasDim=false;
    bool _HasOrder=false;
};