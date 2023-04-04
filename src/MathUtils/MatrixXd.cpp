//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author   : Yang Bai
//+++ Date     : 2020.10.18
//+++ Reviewer : Xiaoyuan @ 2021.08.20
//+++ Purpose  : Define the general Matrix  in AsFem
//+++            we mainly use this for the calculation of jacobian
//+++            If one wants to use Eigen's MatrixXd, please use
//+++            Eigen::MatrixXd, which is different with ours !!!
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MathUtils/MatrixXd.h"

MatrixXd::MatrixXd(){
    m_m=m_n=m_mn=0;
    m_vals.clear();
}
MatrixXd::MatrixXd(const MatrixXd &a){
    m_vals=a.m_vals;
    m_m=a.m_m;m_n=a.m_n;
    m_mn=m_m*m_n;
}
MatrixXd::MatrixXd(const int &m,const int &n){
    m_m=m;m_n=n;
    m_mn=m*n;
    m_vals.resize(m_mn,0.0);
}
MatrixXd::MatrixXd(const int &m,const int &n,const double &val){
    m_m=m;m_n=n;
    m_mn=m*n;
    m_vals.resize(m_mn,val);
}

void MatrixXd::solve(const VectorXd &b,VectorXd &x) const{
    if(b.getM()!=getM()){
        MessagePrinter::printErrorTxt("size of rhs vector b is not equal to the row number of your matrix, can\'t execute the solve function");
        MessagePrinter::exitAsFem();
    }
    if(x.getM()!=getN()){
        MessagePrinter::printErrorTxt("size of solution vector x is not equal to the column number of your matrix, can\'t execute the solve function");
        MessagePrinter::exitAsFem();
    }
    Eigen::MatrixXd A(getM(),getN());
    Eigen::VectorXd B(getM()),X(getN());
    for(int i=1;i<=getM();i++){
        for(int j=1;j<=getN();j++){
            A.coeffRef(i-1,j-1)=(*this)(i,j);
        }
        B.coeffRef(i-1)=b(i);
    }
    X=A.fullPivLu().solve(B);
    for(int i=1;i<=getM();i++){
        x(i)=X.coeff(i-1);
    }
}
VectorXd MatrixXd::solve(const VectorXd &b) const{
    if(b.getM()!=getM()){
        MessagePrinter::printErrorTxt("size of rhs vector b is not equal to the row number of your matrix, can\'t execute the solve function");
        MessagePrinter::exitAsFem();
    }
    VectorXd x(getM(),0.0);
    Eigen::MatrixXd A(getM(),getN());
    Eigen::VectorXd B(getM()),X(getN());
    for(int i=1;i<=getM();i++){
        for(int j=1;j<=getN();j++){
            A.coeffRef(i-1,j-1)=(*this)(i,j);
        }
        B.coeffRef(i-1)=b(i);
    }
    X=A.fullPivLu().solve(B);
    for(int i=1;i<=getM();i++){
        x(i)=X.coeff(i-1);
    }
    return x;
}