//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_MATRIXXD_H
#define ASFEM_MATRIXXD_H

#include <iostream>
#include <iomanip>
#include <vector>

#include "Eigen/Eigen"

using namespace std;


class MatrixXd{
public:
    MatrixXd(){
        _vals.clear();_M=0,_N=0;
    }
    MatrixXd(const int &m,const int &n){
        _vals.resize(m*n,0.0);_M=m;_N=n;_MN=m*n;
    }
    void Resize(const int &m,const int &n){
        _vals.resize(m,n);_M=m;_N=n;_MN=m*n;
    }
    double* GetDataPtr(){
        return _vals.data();
    }
    inline int GetM()const{return _M;}
    inline int GetN()const{return _N;}
    void Clean(){_vals.clear();}
    //*****************************************
    //*** Operator overload
    //*****************************************
    inline double& operator()(const int &i,const int &j){
        return _vals[(i-1)*_N+j-1];
    }
    inline double operator()(const int &i,const int &j)const{
        return _vals[(i-1)*_N+j-1];
    }
    //*****************************************
    //*** For basic mathematic operator
    //*****************************************
    inline MatrixXd& operator=(const double &val){
        for(int i=0;i<_MN;++i) _vals[i]=val;
        return *this;
    }
    inline MatrixXd operator=(const double &val)const{
        MatrixXd temp(GetM(),GetN());
        for(int i=0;i<_MN;++i) temp._vals[i]=val;
        return temp;
    }

    void setZero(){
        for(int i=0;i<_MN;++i) _vals[i]=0.0;
    }
    inline MatrixXd Inverse()const{
        Eigen::MatrixXd Mat(_M,_N),MatInv(_M,_N);
        MatrixXd temp(_M,_N);

        for(int i=1;i<=_M;i++){
            for(int j=1;j<=_N;j++){
                Mat.coeffRef(i-1,j-1)=(*this)(i,j);
            }
        }
        MatInv=Mat.inverse();
        for(int i=1;i<=_M;i++){
            for(int j=1;j<=_N;j++){
                temp(i,j)=MatInv.coeff(i-1,j-1);
            }
        }
        return temp;
    }
    inline double Det()const{
        Eigen::MatrixXd Mat(_M,_N);
        for(int i=1;i<=_M;i++){
            for(int j=1;j<=_N;j++){
                Mat.coeffRef(i-1,j-1)=(*this)(i,j);
            }
        }
        return Mat.determinant();
    }

private:
    vector<double> _vals;
    int _M,_N,_MN;
};



#endif //ASFEM_MATRIXXD_H