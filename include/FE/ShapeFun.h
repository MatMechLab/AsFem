//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_SHAPEFUN_H
#define ASFEM_SHAPEFUN_H


#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>

#include "petsc.h"

//******************************************
#include "MessagePrinter/MessagePrinter.h"
#include "Mesh/MeshType.h"
#include "Mesh/Nodes.h"

#include "Utils/MathUtils.h"
#include "Utils/Vector3d.h"




using namespace std;



class ShapeFun
{
public:
    ShapeFun();
    ShapeFun(int dim,MeshType meshtype);
    void PreCalc();
    void Calc(const double &xi,const Nodes &nodes,const bool &flag);// for 1D case
    void Calc(const double &xi,const double &eta,const Nodes &nodes,const bool &flag);// for 2D case
    void Calc(const double &xi,const double &eta,const double &zeta,const Nodes &nodes,const bool &flag);


    inline double shape_value(const int &i) const{
        return _shape_value[i-1];
    }
    inline Vector3d shape_grad(const int &i) const{
        return _shape_grad[i-1];
    }



    // friend double operator*(const Eigen::Vector3d &a,const Eigen::Vector3d &b){
    //     return a(0)*b(0)+a(1)*b(1)+a(2)*b(2);
    // }


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


private:
    void Compute1DLagrangeShapeFun(const double &xi,const Nodes &nodes,bool flag=true);
    void Compute2DLagrangeShapeFun(const double &xi,const double &eta,const Nodes &nodes,bool flag=true);
    void Compute3DLagrangeShapeFun(const double &xi,const double &eta,const double &zeta,const Nodes &nodes,bool flag=true);


private:
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



#endif //ASFEM_SHAPEFUN_H