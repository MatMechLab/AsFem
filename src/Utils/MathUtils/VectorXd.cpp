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
//+++ Date   : 2020.10.18
//+++ Purpose: Define the general Vector array in AsFem
//+++          we mainly use this for the calculation of residual
//+++          If one wants to use Eigen's VectorXd, please use
//+++          Eigen::VectorXd, which is different with ours !!!
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Utils/VectorXd.h"

VectorXd::VectorXd(){
    _vals.clear();_M=0;
}
VectorXd::VectorXd(const VectorXd &a){
    _vals.resize(a._M);
    _M=a._M;
    for(int i=0;i<a._M;i++) _vals[i]=a._vals[i];
}
VectorXd::VectorXd(const int &m){
    _vals.resize(m,0.0);_M=m;
}
VectorXd::VectorXd(const int &m,const double &val){
    _vals.resize(m,val);_M=m;
}