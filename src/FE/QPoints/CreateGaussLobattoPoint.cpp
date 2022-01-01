//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.07.12
//+++ Purpose: generate the GaussÃ¢ÂÂLobatto points for the integration
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/QPointGaussLobatto.h"

void QPointGaussLobatto::CreateQpoints(MeshType meshtype){
    if (!_HasDim && !_HasOrder && !_HasSettings){
        MessagePrinter::PrintErrorTxt("can\'t generate gauss-lobatto points for the integration, neither the dim, the order or the type is given");
        MessagePrinter::AsFem_Exit();
    }
    if (GetDim() == 1){
        Create1DGaussPoint();
    }
    else if (GetDim() == 2){
        Create2DGaussPoint(meshtype);
    }
    else if (GetDim() == 3){
        Create3DGaussPoint(meshtype);
    }
}
//****************************************************************
//*** For 1D Gauss-Legendre point generation
//****************************************************************
void QPointGaussLobatto::Create1DGaussPoint(){
    // For Gauss-Lobatto type integration points generatioin
    // taken from http://mathworld.wolfram.com/LobattoQuadrature.html
    switch (GetQpOrder()){
        case 0:
        case 1:{
            _nQpPoints=2;
            _QpCoords.resize(_nQpPoints*(GetDim()+1),0.0);

            (*this)(1,1)=-1.0;
            (*this)(2,1)= 1.0;

            (*this)(1,0)= 1.0;
            (*this)(2,0)= 1.0;
            return;
        }
        case 2:
        case 3:{
            _nQpPoints=3;
            _QpCoords.resize(_nQpPoints*(GetDim()+1),0.0);

            (*this)(1,1)=-1.0;
            (*this)(2,1)= 0.0;
            (*this)(3,1)= 1.0;

            (*this)(1,0)=1.0/3.0;
            (*this)(2,0)=4.0/3.0;
            (*this)(3,0)=1.0/3.0;
            return;
        }
        case 4:
        case 5:{
            _nQpPoints=4;
            _QpCoords.resize(_nQpPoints*(GetDim()+1),0.0);

            (*this)(1,1)=-1.0;
            (*this)(2,1)=-sqrt(0.2);
            (*this)(3,1)= sqrt(0.2);
            (*this)(4,1)= 1.0;

            (*this)(1,0)=1.0/6.0;
            (*this)(2,0)=5.0/6.0;
            (*this)(3,0)=5.0/6.0;
            (*this)(4,0)=1.0/6.0;
            return;
        }
        case 6:
        case 7:{
            _nQpPoints=5;
            _QpCoords.resize(_nQpPoints*(GetDim()+1),0.0);

            (*this)(1,1)=-1.0;
            (*this)(2,1)=-sqrt(3.0/7.0);
            (*this)(3,1)= 0.0;
            (*this)(4,1)= sqrt(3.0/7.0);
            (*this)(5,1)= 1.0;

            (*this)(1,0)= 0.1;
            (*this)(2,0)=49.0/90.0;
            (*this)(3,0)=32.0/45.0;
            (*this)(4,0)=49.0/90.0;
            (*this)(5,0)=0.1;
            return;
        }
        case 8:
        case 9:{
            _nQpPoints=6;
            _QpCoords.resize(_nQpPoints*(GetDim()+1),0.0);

            (*this)(1,1)=double(-1.0000000000000000000000000000000e+00L);
            (*this)(2,1)=double(-7.6505532392946469285100297395934e-01L);
            (*this)(3,1)=double(-2.8523151648064509631415099404088e-01L);
            (*this)(4,1)=double( 2.8523151648064509631415099404088e-01L);
            (*this)(5,1)=double( 7.6505532392946469285100297395934e-01L);
            (*this)(6,1)=double( 1.0000000000000000000000000000000e+00L);

            (*this)(1,0)=double(6.6666666666666666666666666666667e-02L);
            (*this)(2,0)=double(3.7847495629784698031661280821202e-01L);
            (*this)(3,0)=double(5.5485837703548635301672052512131e-01L);
            (*this)(4,0)=double(5.5485837703548635301672052512131e-01L);
            (*this)(5,0)=double(3.7847495629784698031661280821202e-01L);
            (*this)(6,0)=double(6.6666666666666666666666666666667e-02L);
            return;
        }
        default:
            string msg = "unsupported gauss-lobatto integration order(=" + to_string(GetQpOrder()) + ") for 1D case";
            MessagePrinter::PrintErrorTxt(msg);
            MessagePrinter::AsFem_Exit();
            break;
    }
}
//****************************************************************
//*** For 2D Gauss-Legendre point generation
//****************************************************************
void QPointGaussLobatto::Create2DGaussPoint(MeshType meshtype){
    switch (meshtype){
        case MeshType::QUAD4:
        case MeshType::QUAD8:
        case MeshType::QUAD9:{
            QPointGaussLobatto qp1d(1,GetQpOrder());
            qp1d.SetQPointType(QPointType::GAUSSLOBATTO);
            qp1d.Create1DGaussPoint();
            _nQpPoints=qp1d.GetQpPointsNum()*qp1d.GetQpPointsNum();
            _QpCoords.resize(_nQpPoints*(GetDim()+1),0.0);
            int k=0;
            // for this regular type mesh , it is just simple tensor product
            for(int i=0;i<qp1d.GetQpPointsNum();i++){
                for(int j=0;j<qp1d.GetQpPointsNum();j++){
                    k+=1;
                    (*this)(k,1)=qp1d(i+1,1);
                    (*this)(k,2)=qp1d(j+1,1);
                    (*this)(k,0)=qp1d(i+1,0)*qp1d(j+1,0);
                }
            }
            return;
        }
        default:
            MessagePrinter::PrintErrorTxt("unsupported mesh type for 2D gauss-lobatto integration");
            MessagePrinter::AsFem_Exit();
            break;
    }
}
//****************************************************************
//*** For 3D Gauss-Legendre point generation
//****************************************************************
void QPointGaussLobatto::Create3DGaussPoint(MeshType meshtype){
    switch (meshtype){
        case MeshType::HEX8:
        case MeshType::HEX20:
        case MeshType::HEX27:{
            QPointGaussLobatto qp1d(1,GetQpOrder());
            qp1d.SetQPointType(QPointType::GAUSSLOBATTO);
            qp1d.Create1DGaussPoint();
            _nQpPoints=qp1d.GetQpPointsNum()*qp1d.GetQpPointsNum()*qp1d.GetQpPointsNum();
            _QpCoords.resize(_nQpPoints*(GetDim()+1),0.0);
            int l=0;
            // for this regular type mesh , it is just simple tensor product
            for(int i=0;i<qp1d.GetQpPointsNum();i++){
                for(int j=0;j<qp1d.GetQpPointsNum();j++){
                    for(int k=0;k<qp1d.GetQpPointsNum();k++){
                        l+=1;
                        (*this)(l,1)=qp1d(i+1,1);
                        (*this)(l,2)=qp1d(j+1,1);
                        (*this)(l,3)=qp1d(k+1,1);
                        (*this)(l,0)=qp1d(i+1,0)*qp1d(j+1,0)*qp1d(k+1,0);
                    }
                }
            }
            return;
        }
        default:
            MessagePrinter::PrintErrorTxt("unsupported mesh type for 3D gauss-lobatto integration");
            MessagePrinter::AsFem_Exit();
            break;
    }
}