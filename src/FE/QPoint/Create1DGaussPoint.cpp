//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "FE/QPoint.h"

void QPoint::Create1DGaussPoint(){
    switch (GetQpOrder()){
        case 0:
        case 1:
        {
            _nQpPoints=1;
            _qp_coords.resize(2,0.0);
            (*this)(1,0)=2.0;
            (*this)(1,1)=0.0;
            return;
        }
        case 2:
        case 3:
        {
            _nQpPoints=2;
             _qp_coords.resize(2*(1+1),0.0);
            (*this)(1,0)=1.0;
            (*this)(1,1)=-sqrt(1.0/3.0);
            (*this)(2,0)=1.0;
            (*this)(2,1)= sqrt(1.0/3.0);
            return;
        }
        case 4:
        case 5:
        {
            _nQpPoints=3;
            _qp_coords.resize(3*(1+1),0.0);
            (*this)(1,0)= 5.0/9.0;
            (*this)(1,1)=-sqrt(0.6);

            (*this)(2,0)= 8.0/9.0;
            (*this)(2,1)= 0.0;

            (*this)(3,0)= 5.0/9.0;
            (*this)(3,1)= sqrt(0.6);
            return;
        }
        case 6:
        case 7:
        {
            _nQpPoints=4;
            _qp_coords.resize(4*(GetDim()+1),0.0);
            const static double t=sqrt(4.8);

            (*this)(1,1)=-sqrt((3.0+t)/7.0);
            (*this)(2,1)=-sqrt((3.0-t)/7.0);
            (*this)(3,1)= sqrt((3.0-t)/7.0);
            (*this)(4,1)= sqrt((3.0+t)/7.0);

            const static double w=1.0/3.0/t;
            (*this)(1,0)=0.5-w;
            (*this)(2,0)=0.5+w;
            (*this)(3,0)=0.5+w;
            (*this)(4,0)=0.5-w;
            return;
        }
        case 8:
        case 9:
        {
            _nQpPoints=5;
            _qp_coords.resize(5*(GetDim()+1),0.0);
            const static double t=sqrt(1120.0);

            (*this)(1,1)=-sqrt((70.0+t)/126.0);
            (*this)(2,1)=-sqrt((70.0-t)/126.0);
            (*this)(3,1)= 0.0;
            (*this)(4,1)= sqrt((70.0-t)/126.0);
            (*this)(5,1)= sqrt((70.0+t)/126.0);

            (*this)(1,0)=(21.0*t+117.60)/(t*(70.0+t));
            (*this)(2,0)=(21.0*t-117.60)/(t*(70.0-t));
            (*this)(3,0)=2.0*(1.0-(*this)(1,0)-(*this)(2,0));
            (*this)(4,0)=(*this)(2,0);
            (*this)(5,0)=(*this)(1,0);
            return;
        }
        case 10:
        case 11:
        {
            _nQpPoints=6;
            _qp_coords.resize(6*(GetDim()+1),0.0);
            (*this)(1,1)=double(-9.3246951420315202781230155449399e-01L);
            (*this)(2,1)=double(-6.6120938646626451366139959501991e-01L);
            (*this)(3,1)=double(-2.3861918608319690863050172168071e-01L);
            (*this)(4,1)=-(*this)(1,1);
            (*this)(5,1)=-(*this)(2,1);
            (*this)(6,1)=-(*this)(3,1);

            (*this)(1,0)=double(1.7132449237917034504029614217273e-01L);
            (*this)(2,0)=double(3.6076157304813860756983351383772e-01L);
            (*this)(3,0)=double(4.6791393457269104738987034398955e-01L);

            (*this)(4,0)=-(*this)(1,0);
            (*this)(5,0)=-(*this)(2,0);
            (*this)(6,0)=-(*this)(3,0);
            return;
        }
        default:
            PetscPrintf(PETSC_COMM_WORLD,"*** Error: invalid gauss integration order(=%3d) for 1D case    !!!   ***\n",GetQpOrder());
            Msg_AsFem_Exit();
            break;
    }
}

//****************************************************
void QPoint::Create1DGaussLobattoPoint(){
    // For Gauss-Lobatto type integration points generatioin
    // taken from http://mathworld.wolfram.com/LobattoQuadrature.html
    switch (GetQpOrder()){
        case 0:
        case 1:
        {
            _nQpPoints=2;
            _qp_coords.resize(_nQpPoints*(GetDim()+1),0.0);

            (*this)(1,1)=-1.0;
            (*this)(2,1)= 1.0;

            (*this)(1,0)= 1.0;
            (*this)(2,0)= 1.0;
            return;
        }
        case 2:
        case 3:
        {
            _nQpPoints=3;
            _qp_coords.resize(_nQpPoints*(GetDim()+1),0.0);

            (*this)(1,1)=-1.0;
            (*this)(2,1)= 0.0;
            (*this)(3,1)= 1.0;

            (*this)(1,0)=1.0/3.0;
            (*this)(2,0)=4.0/3.0;
            (*this)(3,0)=1.0/3.0;
            return;
        }
        case 4:
        case 5:
        {
            _nQpPoints=4;
            _qp_coords.resize(_nQpPoints*(GetDim()+1),0.0);

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
        case 7:
        {
            _nQpPoints=5;
            _qp_coords.resize(_nQpPoints*(GetDim()+1),0.0);

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
        case 9:
        {
            _nQpPoints=6;
            _qp_coords.resize(_nQpPoints*(GetDim()+1),0.0);

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
            PetscPrintf(PETSC_COMM_WORLD,"*** Error: invalid gauss lobatto integration order(=%3d) for 1D case! ***\n",GetQpOrder());
            Msg_AsFem_Exit();
    }
}