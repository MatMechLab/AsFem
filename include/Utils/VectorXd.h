//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
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

#pragma once

#include <iostream>
#include <iomanip>
#include <vector>

#include "Utils/MessagePrinter.h"

using namespace std;

class VectorXd{
public:
    VectorXd();
    VectorXd(const VectorXd &a);
    VectorXd(const int &m);
    VectorXd(const int &m,const double &val);
    void Resize(const int &m){
        _vals.resize(m,0.0);_M=m;
    }
    void Resize(const int &m,const double &val){
        _vals.resize(m,val);_M=m;
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
    //*** for =
    inline VectorXd& operator=(const double &val){
        for(int i=0;i<_M;++i) _vals[i]=val;
        return *this;
    }
    inline VectorXd& operator=(const VectorXd &a){
        if(_M!=a.GetM()){
            MessagePrinter::PrintErrorTxt("a=b cant be applied for two vectors with different size");
            MessagePrinter::AsFem_Exit();
        }
        else{
            for(int i=0;i<_M;++i) _vals[i]=a._vals[i];
            return *this;
        }
        return *this;
    }
    //*** for +
    inline VectorXd operator+(const double &val){
        VectorXd temp(_M);
        for(int i=0;i<_M;++i) temp._vals[i]=_vals[i]+val;
        return temp;
    }
    inline VectorXd operator+(const VectorXd &a){
        if(_M!=a.GetM()){
            MessagePrinter::PrintErrorTxt("a+b cant be applied for two vectors with different size");
            MessagePrinter::AsFem_Exit();
        }
        else{
            VectorXd temp(_M);
            for(int i=0;i<_M;++i) temp._vals[i]=_vals[i]+a._vals[i];
            return temp;
        }
    }
    //*** for +=
    inline VectorXd& operator+=(const double &val){
        for(int i=0;i<_M;++i) _vals[i]+=val;
        return *this;
    }
    inline VectorXd& operator+=(const VectorXd &a){
        if(_M!=a.GetM()){
            MessagePrinter::PrintErrorTxt("a+b cant be applied for two vectors with different size");
            MessagePrinter::AsFem_Exit();
        }
        else{
            for(int i=0;i<_M;++i) _vals[i]+=a._vals[i];
            return *this;
        }
    }
    //************************
    //*** for -
    inline VectorXd operator-(const double &val){
        VectorXd temp(_M);
        for(int i=0;i<_M;++i) temp._vals[i]=_vals[i]-val;
        return temp;
    }
    inline VectorXd operator-(const VectorXd &a){
        if(_M!=a.GetM()){
            MessagePrinter::PrintErrorTxt("a+b cant be applied for two vectors with different size");
            MessagePrinter::AsFem_Exit();
        }
        else{
            VectorXd temp(_M);
            for(int i=0;i<_M;++i) temp._vals[i]=_vals[i]-a._vals[i];
            return temp;
        }
    }
    //*** for -=
    inline VectorXd& operator-=(const double &val){
        for(int i=0;i<_M;++i) _vals[i]-=val;
        return *this;
    }
    inline VectorXd& operator-=(const VectorXd &a){
        if(_M!=a.GetM()){
            MessagePrinter::PrintErrorTxt("a-b cant be applied for two vectors with different size");
            MessagePrinter::AsFem_Exit();
        }
        else{
            for(int i=0;i<_M;++i) _vals[i]-=a._vals[i];
            return *this;
        }
    }
    //***********************************************
    //*** for *
    inline VectorXd operator*(const double &val){
        VectorXd temp(_M);
        for(int i=0;i<_M;++i) temp._vals[i]=_vals[i]*val;
        return temp;
    }
    //*** for *=
    inline VectorXd& operator*=(const double &val){
        for(int i=0;i<_M;++i) _vals[i]*=val;
        return *this;
    }
    //**********************************************
    //*** for /
    inline VectorXd operator/(const double &val){
        VectorXd temp(_M);
        for(int i=0;i<_M;++i) temp._vals[i]=_vals[i]/val;
        return temp;
    }
    //*** for /=
    inline VectorXd& operator/=(const double &val){
        for(int i=0;i<_M;++i) _vals[i]/=val;
        return *this;
    }
    //***********************************************
    void setZero(){
        for(int i=0;i<_M;++i) _vals[i]=0.0;
    }
    void setRandom(){
        srand(time(0));
        for(int i=0;i<_M;++i) _vals[i]=static_cast<double>(1.0*rand()/RAND_MAX);
    }

private:
    vector<double> _vals;
    int _M;
};