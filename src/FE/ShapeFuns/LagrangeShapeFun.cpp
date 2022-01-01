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
//+++ Date   : 2020.11.29
//+++ Purpose: implement the general FE shape functions for FEM calculation
//+++          in this code AsFem offer you:
//+++            1) lagrange shape function in 1d case,i.e. edge2,edge3,edge4  
//+++            2) lagrange shape function for 2d case, i.e. quad4,8,9 and tri3,6 mesh
//+++            3) lagrange shape function for 3d case, i.e. hex8,20,27 and tet4,10 mesh
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/LagrangeShapeFun.h"

LagrangeShapeFun::LagrangeShapeFun(){
    _nOrder=0;_nFuns=0;
    _nDim=0;
    _MeshType=MeshType::EDGE2;
    _IsCartesianDeriv=true;// d/dX,d/dY; false-->d/dxi,d/deta...
    _DetJac=0.0;
    _dxdxi=0.0;_dydxi=0.0;_dzdxi=0.0;
    _dxdeta=0.0;_dydeta=0.0;_dzdeta=0.0;
    _dxdzeta=0.0;_dydzeta=0.0;_dzdzeta=0.0;
    _XJac[0][0]=0.0;_XJac[0][1]=0.0;_XJac[0][2]=0.0;
    _XJac[1][0]=0.0;_XJac[1][1]=0.0;_XJac[1][2]=0.0;
    _XJac[2][0]=0.0;_XJac[2][1]=0.0;_XJac[2][2]=0.0;
    _Jac[0][0]=0.0;_Jac[0][1]=0.0;_Jac[0][2]=0.0;
    _Jac[1][0]=0.0;_Jac[1][1]=0.0;_Jac[1][2]=0.0;
    _Jac[2][0]=0.0;_Jac[2][1]=0.0;_Jac[2][2]=0.0;
    _values.clear();
    _nValues=0;
    _HasDim=false;_HasOrder=false;_HasMeshType=false;

    _shape_value.clear();
    _shape_grad.clear();
}

LagrangeShapeFun::LagrangeShapeFun(int dim,MeshType meshtype){
    _nOrder=0;_nFuns=0;
    _nDim=dim;
    _MeshType=meshtype;
    _IsCartesianDeriv=true;// d/dX,d/dY; false-->d/dxi,d/deta...
    _DetJac=0.0;
    _dxdxi=0.0;_dydxi=0.0;_dzdxi=0.0;
    _dxdeta=0.0;_dydeta=0.0;_dzdeta=0.0;
    _dxdzeta=0.0;_dydzeta=0.0;_dzdzeta=0.0;
    _XJac[0][0]=0.0;_XJac[0][1]=0.0;_XJac[0][2]=0.0;
    _XJac[1][0]=0.0;_XJac[1][1]=0.0;_XJac[1][2]=0.0;
    _XJac[2][0]=0.0;_XJac[2][1]=0.0;_XJac[2][2]=0.0;
    _Jac[0][0]=0.0;_Jac[0][1]=0.0;_Jac[0][2]=0.0;
    _Jac[1][0]=0.0;_Jac[1][1]=0.0;_Jac[1][2]=0.0;
    _Jac[2][0]=0.0;_Jac[2][1]=0.0;_Jac[2][2]=0.0;
    _values.clear();
    _nValues=0;
    _HasDim=true;_HasOrder=false;_HasMeshType=true;

    _shape_value.clear();
    _shape_grad.clear();
}
//**********************************************************
//*** pre-calculation and allocation for memory
//**********************************************************
void LagrangeShapeFun::PreCalc(){
    if(_HasDim&&_HasMeshType){
        switch (GetMeshType()){
            case MeshType::EDGE2:
                _nFuns=2;
                _nOrder=1;
                break;
            case MeshType::EDGE3:
                _nFuns=3;
                _nOrder=2;
                break;
            case MeshType::EDGE4:
                _nFuns=4;
                _nOrder=3;
                break;
            case MeshType::QUAD4:
                _nFuns=4;
                _nOrder=1;
                break;
            case MeshType::QUAD8:
                _nFuns=8;
                _nOrder=2;
                break;
            case MeshType::QUAD9:
                _nFuns=9;
                _nOrder=2;
                break;
            case MeshType::TRI3:
                _nFuns=3;
                _nOrder=1;
                break;
            case MeshType::TRI6:
                _nFuns=6;
                _nOrder=2;
                break;
            case MeshType::HEX8:
                _nFuns=8;
                _nOrder=1;
                break;
            case MeshType::HEX20:
                _nFuns=20;
                _nOrder=2;
                break;
            case MeshType::HEX27:
                _nFuns=27;
                _nOrder=2;
                break;
            case MeshType::TET4:
                _nFuns=4;
                _nOrder=1;
                break;
            case MeshType::TET10:
                _nFuns=10;
                _nOrder=2;
                break;
            default:
                MessagePrinter::PrintErrorTxt("unsupported mesh type in AsFem, currently we only support 1D: edge2, edge3, edge4 shape function, 2D: quad4, quad8, quad9 and tri3, tri6 shape function, 3D: hex8, hex20, hex27 and tet4, tet10 shape function");
                MessagePrinter::AsFem_Exit();
                break;
        }
        _nValues=(1+GetDim())*_nFuns;
        _values.reserve(_nValues);
        for(int i=0;i<_nValues;++i) _values.push_back(0.0);


        _shape_value.reserve(_nFuns);
        _shape_grad.reserve(_nFuns);
        for(int i=0;i<_nFuns;++i){
            _shape_value.push_back(0.0);
            _shape_grad.push_back(Vector3d(0.0));
        }
        for(int i=0;i<_nFuns;++i){
            _shape_grad[i].setZero();
        }
    }
    else{
        MessagePrinter::PrintErrorTxt("information for shape function is not complete");
        MessagePrinter::AsFem_Exit();
    }
}