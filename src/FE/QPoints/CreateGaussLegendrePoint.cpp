//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author  : Yang Bai
//+++ Date    : 2020.07.12
//+++ Reviewer: Xiaoyuan @ 2021.11.21
//+++ Purpose : generate the Gauss Legendre points for the integration
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/QPointGaussLegendre.h"

void QPointGaussLegendre::CreateQpoints(MeshType meshtype){
    if (!_HasDim || !_HasOrder || !_HasSettings){
        MessagePrinter::PrintErrorTxt("can\'t generate gauss-legendre points for the integration, neither the dim, the order or the type is given");
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
void QPointGaussLegendre::Create1DGaussPoint(){
    // taken from https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature
    switch (GetQpOrder())
    {
    case 0:
    case 1:
    {
        _nQpPoints = 1;
        _QpCoords.resize(1*(1+1), 0.0);
        (*this)(1, 0) = 2.0;
        (*this)(1, 1) = 0.0;
        return;
    }
    case 2:
    case 3:
    {
        _nQpPoints = 2;
        _QpCoords.resize(2*(1+1), 0.0);
        (*this)(1, 0) = 1.0;
        (*this)(1, 1) =-sqrt(1.0 / 3.0);

        (*this)(2, 0) = 1.0;
        (*this)(2, 1) = sqrt(1.0 / 3.0);
        return;
    }
    case 4:
    case 5:
    {
        _nQpPoints = 3;
        _QpCoords.resize(3*(1+1), 0.0);
        (*this)(1, 0) = 5.0/9.0;
        (*this)(1, 1) =-sqrt(0.6);

        (*this)(2, 0) = 8.0/9.0;
        (*this)(2, 1) = 0.0;

        (*this)(3, 0) = 5.0/9.0;
        (*this)(3, 1) = sqrt(0.6);
        return;
    }
    case 6:
    case 7:
    {
        _nQpPoints = 4;
        _QpCoords.resize(4*(1+1),0.0);
        const double t1=sqrt(3.0/7.0-(2.0/7.0)*sqrt(6.0/5.0));
        const double t2=sqrt(3.0/7.0+(2.0/7.0)*sqrt(6.0/5.0));
        const double w1=(18.0+sqrt(30.0))/36.0;
        const double w2=(18.0-sqrt(30.0))/36.0;

        (*this)(1,1)=-t2;
        (*this)(2,1)=-t1;
        (*this)(3,1)= t1;
        (*this)(4,1)= t2;

        (*this)(1,0)=w2;
        (*this)(2,0)=w1;
        (*this)(3,0)=w1;
        (*this)(4,0)=w2;
        return;
    }
    case 8:
    case 9:
    {
        _nQpPoints = 5;
        _QpCoords.resize(5*(1+1),0.0);
        const double t1=(1.0/3.0)*sqrt(5.0-2.0*sqrt(10.0/7.0));
        const double t2=(1.0/3.0)*sqrt(5.0+2.0*sqrt(10.0/7.0));
        const double w1=(322.0+13.0*sqrt(70.0))/900.0;
        const double w2=(322.0-13.0*sqrt(70.0))/900.0;

        (*this)(1,1)=-t2;
        (*this)(2,1)=-t1;
        (*this)(3,1)= 0.0;
        (*this)(4,1)= t1;
        (*this)(5,1)= t2;

        (*this)(1,0)=w2;
        (*this)(2,0)=w1;
        (*this)(3,0)=128.0/225.0;
        (*this)(4,0)=w1;
        (*this)(5,0)=w2;
        return;
    }
    case 10:
    case 11:
    {
        // taken from https://pomax.github.io/bezierinfo/legendre-gauss.html
        _nQpPoints = 6;
        _QpCoords.resize(6*(GetDim()+1), 0.0);
        (*this)(1,1)=-0.9324695142031521;
        (*this)(2,1)=-0.6612093864662645;
        (*this)(3,1)=-0.2386191860831969;
        (*this)(4,1)= 0.2386191860831969;
        (*this)(5,1)= 0.6612093864662645;
        (*this)(6,1)= 0.9324695142031521;

        (*this)(1,0)=0.1713244923791704;
        (*this)(2,0)=0.3607615730481386;
        (*this)(3,0)=0.4679139345726910;
        (*this)(4,0)=0.4679139345726910;
        (*this)(5,0)=0.3607615730481386;
        (*this)(6,0)=0.1713244923791704;
        return;
    }
    case 12:
    case 13:
    {
        // taken from https://pomax.github.io/bezierinfo/legendre-gauss.html
        _nQpPoints = 7;
        _QpCoords.resize(7*(GetDim()+1), 0.0);
        (*this)(1,1)=-0.9491079123427585;
        (*this)(2,1)=-0.7415311855993945;
        (*this)(3,1)=-0.4058451513773972;
        (*this)(4,1)= 0.0;
        (*this)(5,1)= 0.4058451513773972;
        (*this)(6,1)= 0.7415311855993945;
        (*this)(7,1)= 0.9491079123427585;

        (*this)(1,0)=0.1294849661688697;
        (*this)(2,0)=0.2797053914892766;
        (*this)(3,0)=0.3818300505051189;
        (*this)(4,0)=0.4179591836734694;
        (*this)(5,0)=0.3818300505051189;
        (*this)(6,0)=0.2797053914892766;
        (*this)(7,0)=0.1294849661688697;
        return;
    }
    case 14:
    case 15:
    {
        // taken from https://pomax.github.io/bezierinfo/legendre-gauss.html
        _nQpPoints = 8;
        _QpCoords.resize(8*(GetDim()+1),0.0);
        (*this)(1,1)=-0.9602898564975363;
        (*this)(2,1)=-0.7966664774136267;
        (*this)(3,1)=-0.5255324099163290;
        (*this)(4,1)=-0.1834346424956498;
        (*this)(5,1)= 0.1834346424956498;
        (*this)(6,1)= 0.5255324099163290;
        (*this)(7,1)= 0.7966664774136267;
        (*this)(8,1)= 0.9602898564975363;

        (*this)(1,0)=0.1012285362903763;
        (*this)(2,0)=0.2223810344533745;
        (*this)(3,0)=0.3137066458778873;
        (*this)(4,0)=0.3626837833783620;
        (*this)(5,0)=0.3626837833783620;
        (*this)(6,0)=0.3137066458778873;
        (*this)(7,0)=0.2223810344533745;
        (*this)(8,0)=0.1012285362903763;
        return;
    }
    default:
        string msg = "unsupported gauss integration order(=" + to_string(GetQpOrder()) + ") for 1D case";
        MessagePrinter::PrintErrorTxt(msg);
        MessagePrinter::AsFem_Exit();
        break;
    }
}
//****************************************************************
//*** For 2D Gauss-Legendre point generation
//****************************************************************
void QPointGaussLegendre::Create2DGaussPoint(MeshType meshtype)
{
    switch (meshtype)
    {
    case MeshType::QUAD4:
    case MeshType::QUAD8:
    case MeshType::QUAD9:
    {
        QPointGaussLegendre qp1d(1, GetQpOrder());
        qp1d.Create1DGaussPoint();
        _nQpPoints = qp1d.GetQpPointsNum() * qp1d.GetQpPointsNum();
        _QpCoords.resize(_nQpPoints * (GetDim() + 1), 0.0);
        int k = 0;
        // for this regular type mesh , it is just simple tensor product
        for (int j = 0; j < qp1d.GetQpPointsNum(); j++){
            for (int i = 0; i < qp1d.GetQpPointsNum(); i++){
                k += 1;
                (*this)(k, 1) = qp1d(i + 1, 1);
                (*this)(k, 2) = qp1d(j + 1, 1);
                (*this)(k, 0) = qp1d(i + 1, 0) * qp1d(j + 1, 0);
            }
        }
        return;
    }
    case MeshType::TRI3:
    case MeshType::TRI6:
    {
        // for the details, one is referred to:
        // http://www.ce.memphis.edu/7111/notes/class_notes/chapter_03d_slides.pdf
        // Thanks for Xiaoyuan's comments, now, the weights are already divided by 2, the correct one.
        switch (GetQpOrder())
        {
        case 0:
        case 1:
        {
            _nQpPoints = 1;
            _QpCoords.resize(_nQpPoints * (GetDim() + 1), 0.0);
            (*this)(1, 1) = 1.0 / 3.0;
            (*this)(1, 2) = 1.0 / 3.0;
            (*this)(1, 0) = 0.5;
            return;
        }
        case 2:
        {
            _nQpPoints = 3;
            _QpCoords.resize(_nQpPoints * (GetDim() + 1), 0.0);
            (*this)(1, 1) = 2.0 / 3.0;
            (*this)(1, 2) = 1.0 / 6.0;
            (*this)(1, 0) = 1.0 / 6.0;

            (*this)(2, 1) = 1.0 / 6.0;
            (*this)(2, 2) = 2.0 / 3.0;
            (*this)(2, 0) = 1.0 / 6.0;

            (*this)(3, 1) = 1.0 / 6.0;
            (*this)(3, 2) = 1.0 / 6.0;
            (*this)(3, 0) = 1.0 / 6.0;
            return;
        }
        case 3:
        {
            _nQpPoints = 4;
            _QpCoords.resize(_nQpPoints * (GetDim() + 1), 0.0);

            (*this)(1, 1) = 1.550510257216821e-01;
            (*this)(1, 2) = 1.785587282636164e-01;
            (*this)(1, 0) = 1.590206908719885e-01;

            (*this)(2, 1) = 6.449489742783178e-01;
            (*this)(2, 2) = 7.503111022260811e-02;
            (*this)(2, 0) = 9.097930912801141e-02;

            (*this)(3, 1) = 1.550510257216821e-01;
            (*this)(3, 2) = 6.663902460147013e-01;
            (*this)(3, 0) = 1.590206908719885e-01;

            (*this)(4, 1) = 6.449489742783178e-01;
            (*this)(4, 2) = 2.800199154990740e-01;
            (*this)(4, 0) = 9.097930912801141e-02;

            return;
        }
        case 4:
        {
            _nQpPoints = 6;
            _QpCoords.resize(_nQpPoints * (GetDim() + 1), 0.0);

            (*this)(1, 1) = 0.0915762135;
            (*this)(1, 2) = 0.8168475730;
            (*this)(1, 0) = 0.1099517437*0.5;

            (*this)(2, 1) = 0.0915762135;
            (*this)(2, 2) = 0.0915762135;
            (*this)(2, 0) = 0.1099517437*0.5;

            (*this)(3, 1) = 0.8168475730;
            (*this)(3, 2) = 0.0915762135;
            (*this)(3, 0) = 0.1099517437*0.5;

            (*this)(4, 1) = 0.4459484909;
            (*this)(4, 2) = 0.1081030182;
            (*this)(4, 0) = 0.2233815897*0.5;

            (*this)(5, 1) = 0.4459484909;
            (*this)(5, 2) = 0.4459484909;
            (*this)(5, 0) = 0.2233815897*0.5;

            (*this)(6, 1) = 0.1081030182;
            (*this)(6, 2) = 0.4459484909;
            (*this)(6, 0) = 0.2233815897*0.5;

            return;
        }
        case 5:
        {
            _nQpPoints = 7;
            _QpCoords.resize(_nQpPoints * (GetDim() + 1), 0.0);

            (*this)(1, 1) = 1.0 / 3.0;
            (*this)(1, 2) = 1.0 / 3.0;
            (*this)(1, 0) = 0.2250000000*0.5;

            (*this)(2, 1) = 0.1012865073;
            (*this)(2, 2) = 0.7974269854;
            (*this)(2, 0) = 0.1259391805*0.5;

            (*this)(3, 1) = 0.1012865073;
            (*this)(3, 2) = 0.1012865073;
            (*this)(3, 0) = 0.1259391805*0.5;

            (*this)(4, 1) = 0.7974269854;
            (*this)(4, 2) = 0.1012865073;
            (*this)(4, 0) = 0.1259391805*0.5;

            (*this)(5, 1) = 0.0597158718;
            (*this)(5, 2) = 0.4701420641;
            (*this)(5, 0) = 0.1323941528*0.5;

            (*this)(6, 1) = 0.4701420641;
            (*this)(6, 2) = 0.4701420641;
            (*this)(6, 0) = 0.1323941528*0.5;

            (*this)(7, 1) = 0.4701420641;
            (*this)(7, 2) = 0.0597158718;
            (*this)(7, 0) = 0.1323941528*0.5;

            return;
        }
        default:
        {
            string msg = "invalid gauss integration order(="+to_string(GetQpOrder())+ ") for 2D case";
            MessagePrinter::PrintErrorTxt(msg);
            MessagePrinter::AsFem_Exit();
            break;
        }
        }
        return;
    }
    default:
    {
        MessagePrinter::PrintErrorTxt("unsupported mesh type for 2D gauss-legendre integration");
        MessagePrinter::AsFem_Exit();
        break;
    }
    }
}
//****************************************************************
//*** For 3D Gauss-Legendre point generation
//****************************************************************
void QPointGaussLegendre::Create3DGaussPoint(MeshType meshtype){
    switch (meshtype){
    case MeshType::HEX8:
    case MeshType::HEX20:
    case MeshType::HEX27:{
        QPointGaussLegendre qp1d(1, GetQpOrder());
        qp1d.Create1DGaussPoint();
        _nQpPoints = qp1d.GetQpPointsNum() * qp1d.GetQpPointsNum() * qp1d.GetQpPointsNum();
        _QpCoords.resize(_nQpPoints * (GetDim() + 1), 0.0);
        int l = 0;
        // for this regular type mesh , it is just simple tensor product
        for (int k = 0; k < qp1d.GetQpPointsNum(); k++){
            for (int j = 0; j < qp1d.GetQpPointsNum(); j++){
                for (int i = 0; i < qp1d.GetQpPointsNum(); i++){
                    l += 1;
                    (*this)(l, 1) = qp1d(i + 1, 1);
                    (*this)(l, 2) = qp1d(j + 1, 1);
                    (*this)(l, 3) = qp1d(k + 1, 1);
                    (*this)(l, 0) = qp1d(i + 1, 0) * qp1d(j + 1, 0) * qp1d(k + 1, 0);
                }
            }
        }
        return;
        break;
    }
    case MeshType::TET4:
    case MeshType::TET10:{
        switch (GetQpOrder()){
        case 0:
        case 1:{
            _nQpPoints = 1;
            _QpCoords.resize(_nQpPoints * (GetDim() + 1), 0.0);
            (*this)(1, 1) = 0.25;
            (*this)(1, 2) = 0.25;
            (*this)(1, 3) = 0.25;
            (*this)(1, 0) = 1.0 / 6.0;
            return;
        }
        case 2:{
            _nQpPoints = 4;
            _QpCoords.resize(_nQpPoints * (GetDim() + 1), 0.0);
            const double a = 0.585410196624969;
            const double b = 0.138196601125011;

            (*this)(1, 1) = a;
            (*this)(1, 2) = b;
            (*this)(1, 3) = b;
            (*this)(1, 0) = 1.0 / 24.0;

            (*this)(2, 1) = b;
            (*this)(2, 2) = a;
            (*this)(2, 3) = b;
            (*this)(2, 0) = 1.0 / 24.0;

            (*this)(3, 1) = b;
            (*this)(3, 2) = b;
            (*this)(3, 3) = a;
            (*this)(3, 0) = 1.0 / 24.0;

            (*this)(4, 1) = b;
            (*this)(4, 2) = b;
            (*this)(4, 3) = b;
            (*this)(4, 0) = 1.0 / 24.0;
            return;
        }
        case 3:{
            _nQpPoints = 5;
            _QpCoords.resize(_nQpPoints * (GetDim() + 1), 0.0);

            (*this)(1, 1) = 0.25;
            (*this)(1, 2) = 0.25;
            (*this)(1, 3) = 0.25;
            (*this)(1, 0) = -2.0 / 15.0;

            (*this)(2, 1) = 0.5;
            (*this)(2, 2) = 1.0 / 6.0;
            (*this)(2, 3) = 1.0 / 6.0;
            (*this)(2, 0) = 0.075;

            (*this)(3, 1) = 1.0 / 6.0;
            (*this)(3, 2) = 0.5;
            (*this)(3, 3) = 1.0 / 6.0;
            (*this)(3, 0) = 0.075;

            (*this)(4, 1) = 1.0 / 6.0;
            (*this)(4, 2) = 1.0 / 6.0;
            (*this)(4, 3) = 0.5;
            (*this)(4, 0) = 0.075;

            (*this)(5, 1) = 1.0 / 6.0;
            (*this)(5, 2) = 1.0 / 6.0;
            (*this)(5, 3) = 1.0 / 6.0;
            (*this)(5, 0) = 0.075;

            return;
        }
        default:{
            string msg="invalid gauss integration order("+to_string(GetQpOrder())+") for 3D case";
            MessagePrinter::PrintErrorTxt(msg);
            MessagePrinter::AsFem_Exit();
            break;
        }
        }
        return;
    }
    default:{
        MessagePrinter::PrintErrorTxt("unsupported mesh type for 3D gauss-legendre integration");
        MessagePrinter::AsFem_Exit();
        break;
    }
    }
}
