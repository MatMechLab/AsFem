//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.10.18
//+++ Purpose: Define the general Vector array in AsFem.
//+++          We mainly use this for the calculation of residual.
//+++          If one wants to use Eigen's VectorXd, please use
//+++          Eigen::VectorXd, which is different with ours !!!
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <cmath>
#include <algorithm>
#include "Utils/MessagePrinter.h"

using std::sqrt;
using std::abs;
using std::fill;

/**
 * This class implements the vector with dynamic size, which is different from the vector3d(fixed size, mainly used in shape function calculation)
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
    void resize(const int &m){
        m_vals.resize(m,0.0);m_m=m;
    }
    /**
     * resize the vector with initial value
     * @param m the size of the vector
     * @param val the intial value of the resized vector
     */
    void resize(const int &m,const double &val){
        m_vals.resize(m,val);m_m=m;
    }
    /**
     * return the size of the vector
     */
    inline int getM()const{return m_m;}
    
    /**
     * clean the whole vector
     */
    void clean(){m_vals.clear();}
    
    /**
     * get the pointer of current vector's data
     */
    double* getDataPtr(){
        return m_vals.data();
    }
    //*****************************************
    //*** Operator overload
    //*****************************************
    /**
     * () operator for the element access of the vector
     * @param i the index of the single element
     */
    inline double& operator()(const int &i){
        if(i<1||i>m_m){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range(m="+to_string(m_m)+")");
            MessagePrinter::exitAsFem();
        }
        return m_vals[i-1];
    }
    /**
     * const () operator for the elemenet access of the vector 
     * @param i the index of the single element
     */
    inline double operator()(const int &i)const{
        if(i<1||i>m_m){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range(m="+to_string(m_m)+")");
            MessagePrinter::exitAsFem();
        }
        return m_vals[i-1];
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
        for(int i=0;i<m_m;++i) m_vals[i]=val;
        return *this;
    }
    /**
     * '=' operator for vector
     * @param a the right hand side vector
     */
    inline VectorXd& operator=(const VectorXd &a){
        if(m_m==0){
            m_m=a.getM();m_vals.resize(m_m,0.0);
            for(int i=0;i<m_m;++i) m_vals[i]=a.m_vals[i];
            return *this;
        }
        else{
            if(m_m==a.getM()){
                for(int i=0;i<m_m;++i) m_vals[i]=a.m_vals[i];
                return *this;
            }
            else{
                MessagePrinter::printErrorTxt("a=b can\'t be applied to two vectors with different size");
                MessagePrinter::exitAsFem();
            }
        }
        return *this;
    }
    //*** for +
    /**
     * '+' operator for scalar
     * @param val the right hand side scalar
     */
    inline VectorXd operator+(const double &val){
        VectorXd temp(m_m);
        for(int i=0;i<m_m;++i) temp.m_vals[i]=m_vals[i]+val;
        return temp;
    }
    /**
     * '+' for vector
     * @param a the right hand side vector
     */
    inline VectorXd operator+(const VectorXd &a){
        VectorXd temp(m_m);
        if(m_m==a.getM()){
            for(int i=0;i<m_m;++i) temp.m_vals[i]=m_vals[i]+a.m_vals[i];
            return temp;
        }
        else{
            MessagePrinter::printErrorTxt("a+b can\'t be applied for two vectors with different size");
            MessagePrinter::exitAsFem();
        }
        return temp;
    }
    //*** for +=
    /**
     * '+=' operator for scalar
     * @param val the right hand side scalar
     */
    inline VectorXd& operator+=(const double &val){
        for(int i=0;i<m_m;++i) m_vals[i]+=val;
        return *this;
    }
    /**
     * '+=' for vector
     * @param a the right hand side vector
     */
    inline VectorXd& operator+=(const VectorXd &a){
        if(m_m==a.getM()){
            for(int i=0;i<m_m;++i) m_vals[i]+=a.m_vals[i];
            return *this;
        }
        else{
            MessagePrinter::printErrorTxt("a+b can\'t be applied for two vectors with different size");
            MessagePrinter::exitAsFem();
        }
        return *this;
    }
    //************************
    //*** for -
    /**
     * '-' operator for scalar
     * @param val right hand side scalar value
     */
    inline VectorXd operator-(const double &val){
        VectorXd temp(m_m);
        for(int i=0;i<m_m;++i) temp.m_vals[i]=m_vals[i]-val;
        return temp;
    }
    /**
     * '-' operator for vector
     * @param a right hand side vector
     */
    inline VectorXd operator-(const VectorXd &a){
        VectorXd temp(m_m);
        if(m_m==a.getM()){
            VectorXd temp(m_m);
            for(int i=0;i<m_m;++i) temp.m_vals[i]=m_vals[i]-a.m_vals[i];
            return temp;
        }
        else{
            MessagePrinter::printErrorTxt("a+b can\'t be applied to two vectors with different size");
            MessagePrinter::exitAsFem();
        }
        return temp;
    }
    //*** for -=
    /**
     * '-' operator for scalar
     * @param val right hand side scalar
     */
    inline VectorXd& operator-=(const double &val){
        for(int i=0;i<m_m;++i) m_vals[i]-=val;
        return *this;
    }
    /**
     * '-=' operator for vector
     * @param a right hand side vector
     */
    inline VectorXd& operator-=(const VectorXd &a){
        if(m_m==a.getM()){
            for(int i=0;i<m_m;++i) m_vals[i]-=a.m_vals[i];
            return *this;
        }
        else{
            MessagePrinter::printErrorTxt("a-b can\'t be applied to two vectors with different size");
            MessagePrinter::exitAsFem();
        }
        return *this;
    }
    //***********************************************
    //*** for *
    /**
     * '*' operator for scalar
     * @param val right hand side scalar
     */
    inline VectorXd operator*(const double &val){
        VectorXd temp(m_m);
        for(int i=0;i<m_m;++i) temp.m_vals[i]=m_vals[i]*val;
        return temp;
    }
    //*** for *=
    /**
     * '*=' operator for scalar
     * @param val the right hand side scalar
     */
    inline VectorXd& operator*=(const double &val){
        for(int i=0;i<m_m;++i) m_vals[i]*=val;
        return *this;
    }
    //**********************************************
    //*** for /
    /**
     * '/' operator for scalar
     * @param val right hand side scalar
     */
    inline VectorXd operator/(const double &val){
        VectorXd temp(m_m);
        if(abs(val)<1.0e-16){
            MessagePrinter::printErrorTxt("val="+to_string(val)+" is singular for / operator in VectorXd");
            MessagePrinter::exitAsFem();
        }
        for(int i=0;i<m_m;++i) temp.m_vals[i]=m_vals[i]/val;
        return temp;
    }
    //*** for /=
    /**
     * '/=' operator for scalar
     * @param val right hand side scalar
     */
    inline VectorXd& operator/=(const double &val){
        if(abs(val)<1.0e-16){
            MessagePrinter::printErrorTxt("val="+to_string(val)+" is singular for /= operator in VectorXd");
            MessagePrinter::exitAsFem();
        }
        for(int i=0;i<m_m;++i) m_vals[i]/=val;
        return *this;
    }
    //***********************************************
    /**
     * set vector's value to be zero
     */
    void setToZero(){
        fill(m_vals.begin(),m_vals.end(),0.0);
    }
    /**
     * set vector's components to be random values
     */
    void setToRandom(){
        srand(time(0));
        for(int i=0;i<m_m;++i) m_vals[i]=static_cast<double>(1.0*rand()/RAND_MAX);
    }
    /**
     * get the L2 norm of current vector
     */
    inline double norm()const{
        double sum=0.0;
        for(int i=0;i<m_m;++i) sum+=static_cast<double>(m_vals[i]*m_vals[i]);
        return sqrt(sum);
    }
    /**
     * get the squared L2 norm of current vector
     */
    inline double normsq()const{
        double sum=0.0;
        for(int i=0;i<m_m;++i) sum+=static_cast<double>(m_vals[i]*m_vals[i]);
        return sum;
    }

private:
    vector<double> m_vals;/**< the double array for vector's components*/
    int m_m;/**< the length of the vector*/
};