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
//+++ Purpose: implement the general gauss integration class for
//+++          AsFem, here one can use:
//+++            1) standard Gauss-Legendre integration scheme
//+++            2) Gauss-Lobatto integration scheme
//+++            3) if you want to implement your own integration
//+++               scheme, all you need to do is create a new class
//+++               which should inherit from QPointBase, then override
//+++               the CreateQPoints() function
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "FE/QPointGaussLegendre.h"
#include "FE/QPointGaussLobatto.h"


class QPoint:public QPointGaussLegendre,public QPointGaussLobatto{
public:
    QPoint();
    QPoint(int dim,int order);
    QPoint(int dim,int order,QPointType qptype);

    //************************************************
    //*** for basic settings
    //************************************************
    void Init();
    void SetDim(int dim){QPointGaussLegendre::SetDim(dim);QPointGaussLobatto::SetDim(dim);}
    void SetQPointOrder(int order){QPointGaussLegendre::SetQPointOrder(order);QPointGaussLobatto::SetQPointOrder(order);}
    void SetQPointType(QPointType qptype){
        QPointGaussLegendre::SetQPointType(qptype);QPointGaussLobatto::SetQPointType(qptype);
        _CurrentQPType=qptype;
    }
    //************************************************
    //*** for basic gettings
    //************************************************
    inline int GetDim() const {
        if(_CurrentQPType==QPointType::GAUSSLEGENDRE){
            return QPointGaussLegendre::GetDim();
        }
        else if(_CurrentQPType==QPointType::GAUSSLOBATTO){
            return QPointGaussLobatto::GetDim();
        }
        return -1;
    }
    inline int GetQpOrder() const {
        if(_CurrentQPType==QPointType::GAUSSLEGENDRE){
            return QPointGaussLegendre::GetQpOrder();
        }
        else if(_CurrentQPType==QPointType::GAUSSLOBATTO){
            return QPointGaussLobatto::GetQpOrder();
        }
        return -1;
    }
    inline int GetQpPointsNum() const {
        if(_CurrentQPType==QPointType::GAUSSLEGENDRE){
            return QPointGaussLegendre::GetQpPointsNum();
        }
        else if(_CurrentQPType==QPointType::GAUSSLOBATTO){
            return QPointGaussLobatto::GetQpPointsNum();
        }
        return -1;
    }
    inline QPointType GetQpPointType()const{return _CurrentQPType;}

    //************************************************
    //*** for operator overload
    //************************************************
    inline double operator()(int i,int j)const{
        if(_CurrentQPType==QPointType::GAUSSLEGENDRE){
            return QPointGaussLegendre::GetIthQpPointJthCoord(i,j);
        }
        else if(_CurrentQPType==QPointType::GAUSSLOBATTO){
            return QPointGaussLobatto::GetIthQpPointJthCoord(i,j);
        }
    }
    inline double& operator()(int i,int j){
        if(_CurrentQPType==QPointType::GAUSSLEGENDRE){
            return QPointGaussLegendre::GetIthQpPointJthCoord(i,j);
        }
        else if(_CurrentQPType==QPointType::GAUSSLOBATTO){
            return QPointGaussLobatto::GetIthQpPointJthCoord(i,j);
        }
        return QPointGaussLegendre::GetIthQpPointJthCoord(i,j);
    }

    inline double GetIthQpPointJthCoord(int i,int j)const{
        if(_CurrentQPType==QPointType::GAUSSLEGENDRE){
            return QPointGaussLegendre::GetIthQpPointJthCoord(i,j);
        }
        else if(_CurrentQPType==QPointType::GAUSSLOBATTO){
            return QPointGaussLobatto::GetIthQpPointJthCoord(i,j);
        }
        return 0.0;
    }

    virtual void CreateQpoints(MeshType meshtype) final;

    void PrintQPointInfo()const;
    void PrintQPointDetailInfo()const;

private:
    QPointType _CurrentQPType;

};