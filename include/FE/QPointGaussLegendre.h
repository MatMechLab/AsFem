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
//+++ Purpose: implement the GaussÃÂ¢ÃÂÃÂLegendre rule for the integration
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "FE/QPointBase.h"

/**
 * this class implement the gauss-legendre point generation
 */
class QPointGaussLegendre:public QPointBase{
public:
    /**
     * constructor for different purpose
     */
    QPointGaussLegendre();
    QPointGaussLegendre(int dim,int order);

    //************************************************
    //*** for basic settings
    //************************************************
    /**
     * set up the dim
     * @param dim dimension numer range from 1~3
     */
    void SetDim(int dim) {
        _nDim=dim;_HasDim=true;
    }
    /**
     * set up the order of gauss integration
     * @param order accuracy order of the gauss integration
     */
    void SetQPointOrder(int order) {
        _nQpOrder=order;_HasOrder=true;
    }
    /**
     * set up the gauss legendre generation method type
     * @param qptype gauss point type
     */
    void SetQPointType(QPointType qptype) {
        _QpType=qptype;_HasSettings=true;
    }
    /**
     * initialize the system
     */
    void Init();
    //************************************************
    //*** for basic gettings
    //************************************************
    /**
     * get the dim of current gauss point
     */
    inline int GetDim() const {return _nDim;}
    /**
     * get the order of current gauss point
     */
    inline int GetQpOrder() const {return _nQpOrder;}
    /**
     * get the gauss point number of current gauss point
     */
    inline int GetQpPointsNum() const {return _nQpPoints;}
    /**
     * get the gauss point type of current gauss point
     */
    inline QPointType GetQpPointType()const{return _QpType;}

    //************************************************
    //*** for operator overload
    //************************************************
    /**
     * get the gauss point cooridnate
     * @param i integer for the i-th qpoint
     * @param j j-th coordinate, 0->weight,1->xi,2->eta,3->zeta
     */
    inline double  operator()(int i,int j)const{return _QpCoords[(i-1)*(_nDim+1)+j];}
    /**
     * get the gauss point cooridnate
     * @param i integer for the i-th qpoint
     * @param j j-th coordinate, 0->weight,1->xi,2->eta,3->zeta
     */
    inline double& operator()(int i,int j){return _QpCoords[(i-1)*(_nDim+1)+j];}

    /**
     * get the gauss point cooridnate
     * @param i integer for the i-th qpoint
     * @param j j-th coordinate, 0->weight,1->xi,2->eta,3->zeta
     */
    inline double& GetIthQpPointJthCoord(int i,int j){return _QpCoords[(i-1)*(_nDim+1)+j];}
    /**
     * get the gauss point cooridnate
     * @param i integer for the i-th qpoint
     * @param j j-th coordinate, 0->weight,1->xi,2->eta,3->zeta
     */
    inline double  GetIthQpPointJthCoord(int i,int j)const{return _QpCoords[(i-1)*(_nDim+1)+j];}

    /**
     * generate the gauss points
     */
    virtual void CreateQpoints(MeshType meshtype) override;

private:
    /**
     * for 1d gauss point generation
     */
    void Create1DGaussPoint();
    /**
     * for 2d gauss point generation
     */
    void Create2DGaussPoint(MeshType meshtype);
    /**
     * for 3d gauss point generation
     */
    void Create3DGaussPoint(MeshType meshtype);


    //***********************************************
    //*** for some basic variables
    //***********************************************
    bool _HasSettings=false;
    bool _HasDim=false;
    bool _HasOrder=false;

};
