//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_VECTORXD_H
#define ASFEM_VECTORXD_H

#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;


class VectorXd{
public:
    VectorXd(){
        _vals.clear();_M=0;
    }
    VectorXd(const int &m){
        _vals.resize(m,0.0);_M=m;
    }
    void Resize(const int &m){
        _vals.resize(m);_M=m;
    }
    inline int GetM()const{return _M;}
    void Clean(){_vals.clear();}
    double* GetDataPtr(){
        return _vals.data();
    }
    //*****************************************
    //*** Operator overload
    //*****************************************
    inline double& operator()(const int &i){
        return _vals[i-1];
    }
    inline double operator()(const int &i)const{
        return _vals[i-1];
    }
    //*****************************************
    //*** For basic mathematic operator
    //*****************************************
    inline VectorXd& operator=(const double &val){
        for(int i=0;i<_M;++i) _vals[i]=val;
        return *this;
    }
    inline VectorXd operator=(const double &val)const{
        VectorXd temp(GetM());
        for(int i=0;i<_M;++i) temp._vals[i]=val;
        return temp;
    }

    void setZero(){
        for(int i=0;i<_M;++i) _vals[i]=0.0;
    }

private:
    vector<double> _vals;
    int _M;
};


#endif //ASFEM_VECTORXD_H