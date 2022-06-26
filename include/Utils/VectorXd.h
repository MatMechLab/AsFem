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

/**
 * This class implement the vector with dynamic size, which different from the vector3d(fixed size, mainly used in shape function calculation)
 */
class VectorXd{
public:
    /**
     * construtor for different purpose
     */
    VectorXd();
    VectorXd(const VectorXd &a);
    VectorXd(const int &m);
    VectorXd(const int &m,const double &val);
    
    /**
     * resize the vector
     * @param m the size of the vector
     */
    void Resize(const int &m){
        _vals.resize(m,0.0);_M=m;
    }
    /**
     * resize the vector with initial value
     * @param m the size of the vector
     * @param val the intial value of the resized vector
     */
    void Resize(const int &m,const double &val){
        _vals.resize(m,val);_M=m;
    }
    /**
     * return the size of the vector
     */
    inline int GetM()const{return _M;}
    
    /**
     * clean the whole vector
     */
    void Clean(){_vals.clear();_M=0;}
    
    /**
     * get the pointer of the vector's data
     */
    double* GetDataPtr(){
        return _vals.data();
    }
    //*****************************************
    //*** Operator overload
    //*****************************************
    /**
     * () operator for the element access of the vector
     * @param i the index of the single element
     */
    inline double& operator()(const int &i){
        return _vals[i-1];
    }
    /**
     * const () operator for the elemenet access of the vector 
     * @param i the index of the single element
     */
    inline double operator()(const int &i)const{
        return _vals[i-1];
    }
    //*****************************************
    //*** For basic mathematic operator
    //*****************************************
    //*** for =
    /**
     * '=' operator for scalar
     * @param val right hand side scalar
     */
    inline VectorXd& operator=(const double &val){
        fill(_vals.begin(),_vals.end(),val);
        return *this;
    }
    /**
     * '=' operator for vector
     * @param a the right hand side vector
     */
    inline VectorXd& operator=(const VectorXd &a){
        if(_M==0){
            _M=a.GetM();_vals.resize(_M,0.0);
            for(int i=0;i<_M;++i) _vals[i]=a._vals[i];
            return *this;
        }
        else{
            if(_M==a.GetM()){
                for(int i=0;i<_M;++i) _vals[i]=a._vals[i];
                return *this;
            }
            else{
                MessagePrinter::PrintErrorTxt("a=b cant be applied for two vectors with different size");
                MessagePrinter::AsFem_Exit();
            }
        }
        return *this;
    }
    //*** for +
    /**
     * '+' operator for scalar
     * @param val the right hand side scalar
     */
    inline VectorXd operator+(const double &val)const{
        VectorXd temp(_M);
        for(int i=0;i<_M;++i) temp._vals[i]=_vals[i]+val;
        return temp;
    }
    /**
     * '+' for vector
     * @param a the right hand side vector
     */
    inline VectorXd operator+(const VectorXd &a)const{
        VectorXd temp(_M);
        if(_M==a.GetM()){
            for(int i=0;i<_M;++i) temp._vals[i]=_vals[i]+a._vals[i];
            return temp;
        }
        else{
            MessagePrinter::PrintErrorTxt("a+b cant be applied for two vectors with different size");
            MessagePrinter::AsFem_Exit();
        }
        return temp;
    }
    //*** for +=
    /**
     * '+=' operator for scalar
     * @param val the right hand side scalar
     */
    inline VectorXd& operator+=(const double &val){
        for(int i=0;i<_M;++i) _vals[i]+=val;
        return *this;
    }
    /**
     * '+=' for vector
     * @param a the right hand side vector
     */
    inline VectorXd& operator+=(const VectorXd &a){
        if(_M==a.GetM()){
            for(int i=0;i<_M;++i) _vals[i]+=a._vals[i];
            return *this;
        }
        else{
            MessagePrinter::PrintErrorTxt("a+b cant be applied for two vectors with different size");
            MessagePrinter::AsFem_Exit();
        }
        return *this;
    }
    //************************
    //*** for -
    /**
     * '-' operator for scalar
     * @param val right hand side scalar value
     */
    inline VectorXd operator-(const double &val)const{
        VectorXd temp(_M);
        for(int i=0;i<_M;++i) temp._vals[i]=_vals[i]-val;
        return temp;
    }
    /**
     * '-' operator for vector
     * @param a right hand side vector
     */
    inline VectorXd operator-(const VectorXd &a)const{
        VectorXd temp(_M);
        if(_M==a.GetM()){
            VectorXd temp(_M);
            for(int i=0;i<_M;++i) temp._vals[i]=_vals[i]-a._vals[i];
            return temp;
        }
        else{
            MessagePrinter::PrintErrorTxt("a+b cant be applied for two vectors with different size");
            MessagePrinter::AsFem_Exit();
        }
        return temp;
    }
    //*** for -=
    /**
     * '-' operator for scalar
     * @param val right hand side scalar
     */
    inline VectorXd& operator-=(const double &val){
        for(int i=0;i<_M;++i) _vals[i]-=val;
        return *this;
    }
    /**
     * '-=' operator for vector
     * @param a right hand side vector
     */
    inline VectorXd& operator-=(const VectorXd &a){
        if(_M==a.GetM()){
            for(int i=0;i<_M;++i) _vals[i]-=a._vals[i];
            return *this;
        }
        else{
            MessagePrinter::PrintErrorTxt("a-b cant be applied for two vectors with different size");
            MessagePrinter::AsFem_Exit();
        }
        return *this;
    }
    //***********************************************
    //*** for *
    /**
     * '*' operator for scalar
     * @param val right hand side scalar
     */
    inline VectorXd operator*(const double &val)const{
        VectorXd temp(_M);
        for(int i=0;i<_M;++i) temp._vals[i]=_vals[i]*val;
        return temp;
    }
    //*** for *=
    /**
     * '*=' operator for scalar
     * @param val the right hand side scalar
     */
    inline VectorXd& operator*=(const double &val){
        for(int i=0;i<_M;++i) _vals[i]*=val;
        return *this;
    }
    //**********************************************
    //*** for /
    /**
     * '/' operator for scalar
     * @param val right hand side scalar
     */
    inline VectorXd operator/(const double &val)const{
        if(abs(val)<1.0e-16){
            MessagePrinter::PrintErrorTxt("x/0 is not acceptable for '/' operator");
            MessagePrinter::AsFem_Exit();
        }
        VectorXd temp(_M);
        for(int i=0;i<_M;++i) temp._vals[i]=_vals[i]/val;
        return temp;
    }
    //*** for /=
    /**
     * '/=' operator for scalar
     * @param val right hand side scalar
     */
    inline VectorXd& operator/=(const double &val){
        if(abs(val)<1.0e-16){
            MessagePrinter::PrintErrorTxt("x/0 is not acceptable for '/' operator");
            MessagePrinter::AsFem_Exit();
        }
        for(int i=0;i<_M;++i) _vals[i]/=val;
        return *this;
    }
    //***********************************************
    /**
     * set vector's value to be zero
     */
    void setZero(){
        fill(_vals.begin(),_vals.end(),0.0);
    }

    /**
     * set vector's components to be random values
     */
    void setRandom(){
        srand(time(0));
        for(int i=0;i<_M;++i) _vals[i]=static_cast<double>(1.0*rand()/RAND_MAX);
    }

private:
    vector<double> _vals;/**< the double array for vector's components*/
    int _M;/**< the length of the vector*/
};
