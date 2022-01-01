//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author  : Yang Bai
//+++ Date    : 2020.10.18
//+++ Reviewer: Xiaoyuan @ 2021.08.20
//+++ Purpose : Define the general Matrix  in AsFem
//+++           we mainly use this for the calculation of jacobian
//+++           If one wants to use Eigen's MatrixXd, please use
//+++           Eigen::MatrixXd, which is different with ours !!!
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Utils/MatrixXd.h"

MatrixXd::MatrixXd(){
    _vals.clear();_M=0,_N=0;
}
MatrixXd::MatrixXd(const MatrixXd &a){
    _vals.resize(a._M*a._N,0.0);
    _M=a._M;
    _N=a._N;
    _MN=a._MN;
    for(int i=0;i<a._MN;i++) _vals[i]=a._vals[i];
}
MatrixXd::MatrixXd(const int &m,const int &n){
    _vals.resize(m*n,0.0);_M=m;_N=n;_MN=m*n;
}
MatrixXd::MatrixXd(const int &m,const int &n,const double &val){
    _vals.resize(m*n,val);_M=m;_N=n;_MN=m*n;
}
