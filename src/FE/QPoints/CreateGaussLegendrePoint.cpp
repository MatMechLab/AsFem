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
//+++ Purpose: generate the GaussâLegendre points for the integration
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/QPointGaussLegendre.h"

void QPointGaussLegendre::CreateQpoints(MeshType meshtype)
{
    if (!_HasDim && !_HasOrder && !_HasSettings){
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
void QPointGaussLegendre::Create1DGaussPoint()
{
    switch (GetQpOrder())
    {
    case 0:
    case 1:
    {
        _nQpPoints = 1;
        _QpCoords.resize(2, 0.0);
        (*this)(1, 0) = 2.0;
        (*this)(1, 1) = 0.0;
        return;
    }
    case 2:
    case 3:
    {
        _nQpPoints = 2;
        _QpCoords.resize(2 * (1 + 1), 0.0);
        (*this)(1, 0) = 1.0;
        (*this)(1, 1) = -sqrt(1.0 / 3.0);
        (*this)(2, 0) = 1.0;
        (*this)(2, 1) = sqrt(1.0 / 3.0);
        return;
    }
    case 4:
    case 5:
    {
        _nQpPoints = 3;
        _QpCoords.resize(3 * (1 + 1), 0.0);
        (*this)(1, 0) = 5.0 / 9.0;
        (*this)(1, 1) = -sqrt(0.6);
        (*this)(2, 0) = 8.0 / 9.0;
        (*this)(2, 1) = 0.0;

        (*this)(3, 0) = 5.0 / 9.0;
        (*this)(3, 1) = sqrt(0.6);
        return;
    }
    case 6:
    case 7:
    {
        _nQpPoints = 4;
        _QpCoords.resize(4 * (GetDim() + 1), 0.0);
        const static double t = sqrt(4.8);

        (*this)(1, 1) = -sqrt((3.0 + t) / 7.0);
        (*this)(2, 1) = -sqrt((3.0 - t) / 7.0);
        (*this)(3, 1) = sqrt((3.0 - t) / 7.0);
        (*this)(4, 1) = sqrt((3.0 + t) / 7.0);

        const static double w = 1.0 / 3.0 / t;
        (*this)(1, 0) = 0.5 - w;
        (*this)(2, 0) = 0.5 + w;
        (*this)(3, 0) = 0.5 + w;
        (*this)(4, 0) = 0.5 - w;
        return;
    }
    case 8:
    case 9:
    {
        _nQpPoints = 5;
        _QpCoords.resize(5 * (GetDim() + 1), 0.0);
        const static double t = sqrt(1120.0);

        (*this)(1, 1) = -sqrt((70.0 + t) / 126.0);
        (*this)(2, 1) = -sqrt((70.0 - t) / 126.0);
        (*this)(3, 1) = 0.0;
        (*this)(4, 1) = sqrt((70.0 - t) / 126.0);
        (*this)(5, 1) = sqrt((70.0 + t) / 126.0);

        (*this)(1, 0) = (21.0 * t + 117.60) / (t * (70.0 + t));
        (*this)(2, 0) = (21.0 * t - 117.60) / (t * (70.0 - t));
        (*this)(3, 0) = 2.0 * (1.0 - (*this)(1, 0) - (*this)(2, 0));
        (*this)(4, 0) = (*this)(2, 0);
        (*this)(5, 0) = (*this)(1, 0);
        return;
    }
    case 10:
    case 11:
    {
        _nQpPoints = 6;
        _QpCoords.resize(6 * (GetDim() + 1), 0.0);
        (*this)(1, 1) = double(-9.3246951420315202781230155449399e-01L);
        (*this)(2, 1) = double(-6.6120938646626451366139959501991e-01L);
        (*this)(3, 1) = double(-2.3861918608319690863050172168071e-01L);
        (*this)(4, 1) = -(*this)(1, 1);
        (*this)(5, 1) = -(*this)(2, 1);
        (*this)(6, 1) = -(*this)(3, 1);

        (*this)(1, 0) = double(1.7132449237917034504029614217273e-01L);
        (*this)(2, 0) = double(3.6076157304813860756983351383772e-01L);
        (*this)(3, 0) = double(4.6791393457269104738987034398955e-01L);

        (*this)(4, 0) = -(*this)(1, 0);
        (*this)(5, 0) = -(*this)(2, 0);
        (*this)(6, 0) = -(*this)(3, 0);
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
        for (int j = 0; j < qp1d.GetQpPointsNum(); j++)
        {
            for (int i = 0; i < qp1d.GetQpPointsNum(); i++)
            {
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
        switch (GetQpOrder())
        {
        case 0:
        case 1:
        {
            _nQpPoints = 1;
            _QpCoords.resize(_nQpPoints * (GetDim() + 1), 0.0);
            (*this)(1, 1) = 1.0 / 3.0;
            (*this)(1, 2) = 1.0 / 3.0;
            (*this)(1, 0) = 1.0;
            return;
        }
        case 2:
        {
            _nQpPoints = 3;
            _QpCoords.resize(_nQpPoints * (GetDim() + 1), 0.0);
            (*this)(1, 1) = 2.0 / 3.0;
            (*this)(1, 2) = 1.0 / 6.0;
            (*this)(1, 0) = 1.0 / 3.0;

            (*this)(2, 1) = 1.0 / 6.0;
            (*this)(2, 2) = 2.0 / 3.0;
            (*this)(2, 0) = 1.0 / 3.0;

            (*this)(3, 1) = 1.0 / 6.0;
            (*this)(3, 2) = 1.0 / 6.0;
            (*this)(3, 0) = 1.0 / 3.0;
            return;
        }
        case 3:
        {
            _nQpPoints = 4;
            _QpCoords.resize(_nQpPoints * (GetDim() + 1), 0.0);

            (*this)(1, 1) = double(1.5505102572168219018027159252941e-01L);
            (*this)(1, 2) = double(1.7855872826361642311703513337422e-01L);
            (*this)(1, 0) = double(1.5902069087198858469718450103758e-01L);

            (*this)(2, 1) = double(6.4494897427831780981972840747059e-01L);
            (*this)(2, 2) = double(7.5031110222608118177475598324603e-02L);
            (*this)(2, 0) = double(9.0979309128011415302815498962418e-02L);

            (*this)(3, 1) = double(1.5505102572168219018027159252941e-01L);
            (*this)(3, 2) = double(6.6639024601470138670269327409637e-01L);
            (*this)(3, 0) = double(1.5902069087198858469718450103758e-01L);

            (*this)(4, 1) = double(6.4494897427831780981972840747059e-01L);
            (*this)(4, 2) = double(2.8001991549907407200279599420481e-01L);
            (*this)(4, 0) = double(9.0979309128011415302815498962418e-02L);

            return;
        }
        case 4:
        {
            _nQpPoints = 6;
            _QpCoords.resize(_nQpPoints * (GetDim() + 1), 0.0);

            (*this)(1, 1) = 0.0915762135;
            (*this)(1, 2) = 0.8168475730;
            (*this)(1, 0) = 0.1099517437;

            (*this)(2, 1) = 0.0915762135;
            (*this)(2, 2) = 0.0915762135;
            (*this)(2, 0) = 0.1099517437;

            (*this)(3, 1) = 0.8168475730;
            (*this)(3, 2) = 0.0915762135;
            (*this)(3, 0) = 0.1099517437;

            (*this)(4, 1) = 0.4459484909;
            (*this)(4, 2) = 0.1081030182;
            (*this)(4, 0) = 0.2233815897;

            (*this)(5, 1) = 0.4459484909;
            (*this)(5, 2) = 0.4459484909;
            (*this)(5, 0) = 0.2233815897;

            (*this)(6, 1) = 0.1081030182;
            (*this)(6, 2) = 0.4459484909;
            (*this)(6, 0) = 0.2233815897;

            return;
        }
        case 5:
        {
            _nQpPoints = 7;
            _QpCoords.resize(_nQpPoints * (GetDim() + 1), 0.0);

            (*this)(1, 1) = 1.0 / 3.0;
            (*this)(1, 2) = 1.0 / 3.0;
            (*this)(1, 0) = 0.2250000000;

            (*this)(2, 1) = 0.1012865073;
            (*this)(2, 2) = 0.7974269854;
            (*this)(2, 0) = 0.1259391805;

            (*this)(3, 1) = 0.1012865073;
            (*this)(3, 2) = 0.1012865073;
            (*this)(3, 0) = 0.1259391805;

            (*this)(4, 1) = 0.7974269854;
            (*this)(4, 2) = 0.1012865073;
            (*this)(4, 0) = 0.1259391805;

            (*this)(5, 1) = 0.0597158718;
            (*this)(5, 2) = 0.4701420641;
            (*this)(5, 0) = 0.1323941528;

            (*this)(6, 1) = 0.4701420641;
            (*this)(6, 2) = 0.4701420641;
            (*this)(6, 0) = 0.1323941528;

            (*this)(7, 1) = 0.4701420641;
            (*this)(7, 2) = 0.0597158718;
            (*this)(7, 0) = 0.1323941528;

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
                    // cout<<"gpInd="<<l
                    //     <<": xi="<<(*this)(l,1)
                    //     <<", eta="<<(*this)(l,2)
                    //     <<": zeta="<<(*this)(l,3)<<endl;
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