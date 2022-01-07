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
//+++ Date   : 2020.11.29
//+++ Purpose: implement the general FE shape functions for FEM calculation
//+++          in this code AsFem offer you:
//+++            1) lagrange shape function in 1d case,i.e. edge2,edge3,edge4  
//+++            2) lagrange shape function for 2d case, i.e. quad4,8,9 and tri3,6 mesh
//+++            3) lagrange shape function for 3d case, i.e. hex8,20,27 and tet4,10 mesh
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>

#include "Utils/Vector3d.h"
#include "Utils/MessagePrinter.h"

#include "Mesh/MeshType.h"
#include "Mesh/Nodes.h"


#include "FE/Lagrange1DShapeFun.h"
#include "FE/Lagrange2DShapeFun.h"
#include "FE/Lagrange3DShapeFun.h"

using namespace std;


class LagrangeShapeFun: public Lagrange1DShapeFun,
                        public Lagrange2DShapeFun,
                        public Lagrange3DShapeFun
{
public:
    LagrangeShapeFun();
    LagrangeShapeFun(int dim,MeshType meshtype);
    void PreCalc();
    void Calc(const double &xi,const Nodes &nodes,const bool &flag);// for 1D case
    void Calc(const double &xi,const double &eta,const Nodes &nodes,const bool &flag);// for 2D case
    void Calc(const double &xi,const double &eta,const double &zeta,const Nodes &nodes,const bool &flag); // for 3D case

    inline double& shape_value(const int &i){
        return _shape_value[i-1];
    }
    inline double shape_value(const int &i) const{
        return _shape_value[i-1];
    }
    inline Vector3d& shape_grad(const int &i){
        return _shape_grad[i-1];
    }
    inline Vector3d shape_grad(const int &i) const{
        return _shape_grad[i-1];
    }

// settings
    void SetDim(int dim) {_nDim=dim;}
    void SetShapeFunNum(int nnodes) {_nFuns=nnodes;}
    void SetShapeFunType(MeshType meshtype);

    // get information
    inline int GetShapeFunNums() const {return _nFuns;}
    inline int GetDim() const {return _nDim;}
    inline int GetOrder() const {return _nOrder;}
    inline double GetDetJac() const {return _DetJac;}
    inline MeshType GetMeshType() const {return _MeshType;}

    // operator overload
    inline double  operator()(const int &i,const int &j) const {return _values[(i-1)*(_nDim+1)+j];}
    inline double& operator()(const int &i,const int &j) {return _values[(i-1)*(_nDim+1)+j];}


protected:
    int _nOrder,_nFuns;
    int _nDim;
    MeshType _MeshType;
    bool _IsCartesianDeriv;
    double _DetJac;
    double _dxdxi,_dydxi,_dzdxi;
    double _dxdeta,_dydeta,_dzdeta;
    double _dxdzeta,_dydzeta,_dzdzeta;
    double _XJac[3][3],_Jac[3][3];
    vector<double> _values;
    int _nValues;
    bool _HasDim,_HasOrder,_HasMeshType;

    vector<double> _shape_value;
    vector<Vector3d> _shape_grad;
};
