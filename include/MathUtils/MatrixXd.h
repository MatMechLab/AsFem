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

#pragma once

#include "Eigen/Eigen"

#include "Utils/MessagePrinter.h"
#include "MathUtils/VectorXd.h"

using std::fill;

/**
 * This class implements the general matrix calculation, it shoul be noted that the vector
 * should be one of the special case of this class
 */
class MatrixXd{
public:
    /**
     * Constructor function for different purpose
     */
    MatrixXd();
    MatrixXd(const MatrixXd &a);
    MatrixXd(const int &m,const int &n);
    MatrixXd(const int &m,const int &n,const double &val);
    
    /**
     * Resize the matrix, the memory is reallocated, the matrix is set to zero by default
     * @param m integer, the size of the 1st dimention
     * @param n integra, the size of the 2nd dimention
     */
    void resize(const int &m,const int &n){
        m_vals.resize(m*n,0.0);m_m=m;m_n=n;m_mn=m*n;
    }

    /**
     * Resize the matrix with a initial value, the memory is reallocated
     * @param m integer, the size of the 1st dimention
     * @param n integra, the size of the 2nd dimention
     * @param val the initial value for the resized matrix
     */
    void resize(const int &m,const int &n,const double &val){
        m_vals.resize(m*n,val);m_m=m;m_n=n;m_mn=m*n;
    }

    /**
     * Return the pointer of the matrix's data (its a vector<double> type)
     */
    double* getDataPtr(){
        return m_vals.data();
    }

    /**
     * Return the size of the 1st dimension
     */
    inline int getM()const{return m_m;}
    /**
     * Return the size of the 2nd dimension
     */
    inline int getN()const{return m_n;}
    /**
     * Clean the whole matrix data
     */
    void clean(){m_vals.clear();}
    //*****************************************
    //*** Operator overload
    //*****************************************
    /**
     * The () operator for the data access
     * @param i the index of the 1st dimension, it should start from 1, not 0!!!
     * @param j the index of the 2nd dimension, it should start from 1, not 0!!!
     */
    inline double& operator()(const int &i,const int &j){
        if(i<1||i>m_m){
            MessagePrinter::printErrorTxt("i= "+to_string(i)+" is out of range(m="+to_string(m_m)+") in MatrixXd.h");
            MessagePrinter::exitAsFem();
        }
        if(j<1||i>m_n){
            MessagePrinter::printErrorTxt("j= "+to_string(j)+" is out of range(m="+to_string(m_n)+") in MatrixXd.h");
            MessagePrinter::exitAsFem();
        }
        return m_vals[(i-1)*m_n+j-1];
    }
    /**
     * The () operator for the data access with constant reference(not editable!)
     * @param i the index of the 1st dimension, it should start from 1, not 0!!!
     * @param j the index of the 2nd dimension, it should start from 1, not 0!!!
     */
    inline double operator()(const int &i,const int &j)const{
        if(i<1||i>m_m){
            MessagePrinter::printErrorTxt("i= "+to_string(i)+" is out of range(m="+to_string(m_m)+") in MatrixXd.h");
            MessagePrinter::exitAsFem();
        }
        if(j<1||i>m_n){
            MessagePrinter::printErrorTxt("j= "+to_string(j)+" is out of range(m="+to_string(m_n)+") in MatrixXd.h");
            MessagePrinter::exitAsFem();
        }
        return m_vals[(i-1)*m_n+j-1];
    }
    /**
     * The [] operator, the data of our matrix is just simple vector in 1D
     * @param i the index of the data vector element, it should start from 1, not 0!!!
     */
    inline double& operator[](const int &i){
        if(i<1||i>m_mn){
            MessagePrinter::printErrorTxt("i= "+to_string(i)+" is out of range(mn="+to_string(m_mn)+") in MatrixXd.h");
            MessagePrinter::exitAsFem();
        }
        return m_vals[i-1];
    }
    /**
     * The [] operator with constant reference, the data of our matrix is just simple vector in 1D
     * @param i the index of the data vector element, it should start from 1, not 0!!!
     */
    inline double operator[](const int &i)const{
        if(i<1||i>m_mn){
            MessagePrinter::printErrorTxt("i= "+to_string(i)+" is out of range(mn="+to_string(m_mn)+") in MatrixXd.h");
            MessagePrinter::exitAsFem();
        }
        return m_vals[i-1];
    }
    //*****************************************
    //*** For basic mathematic operator
    //*****************************************
    //*** for =
    /**
     * The '=' for equal operator
     * @param val the double type value to set up the whole matrix
     */
    inline MatrixXd& operator=(const double &val){
        fill(m_vals.begin(),m_vals.end(),val);
        return *this;
    }
    /**
     * The '=' for equal operator between two (same) matrix
     * @param a the right-hand side matrix with the same dimensions as current one
     */
    inline MatrixXd& operator=(const MatrixXd &a){
        if(m_m==0&&m_n==0){
            m_m=a.getM();m_n=a.getN();
            m_mn=m_m*m_n;m_vals.resize(m_mn,0.0);
            for(int i=0;i<m_mn;++i) m_vals[i]=a.m_vals[i];
            return *this;
        }
        else{
            if(m_m==a.getM()&&m_n==a.getN()){
                for(int i=0;i<m_mn;++i) m_vals[i]=a.m_vals[i];
                return *this;
            }
            else{
                MessagePrinter::printErrorTxt("a=b can\'t be applied to two matrix with different size");
                MessagePrinter::exitAsFem();
            }
        }
        return *this;
    }
    //****************************
    //*** for +
    /**
     * The '+' operator between matrix and scalar
     * @param val the right-hand side scalar (double type)
     */
    inline MatrixXd operator+(const double &val)const{
        MatrixXd temp(m_m,m_n);
        for(int i=0;i<m_mn;++i) temp.m_vals[i]=m_vals[i]+val;
        return temp;
    }
    /**
     * The '+' operator between matrix and matrix
     * @param a the right-hand side matrix (the dimensions should be the same)
     */
    inline MatrixXd operator+(const MatrixXd &a)const{
        MatrixXd temp(m_m,m_n);
        if(m_m==a.getM()&&m_n==a.getN()){
            for(int i=0;i<m_mn;++i) temp.m_vals[i]=m_vals[i]+a.m_vals[i];
            return temp;
        }
        else{
            MessagePrinter::printErrorTxt("a+b can\'t be applied to two matrix with different size");
            MessagePrinter::exitAsFem();
        }
        return temp;
    }
    //*** for +=
    /**
     * The '+=' operator between matrix and scalar
     * @param val the right-hand side scalar (double type)
     */
    inline MatrixXd& operator+=(const double &val){
        for(int i=0;i<m_mn;++i) m_vals[i]=m_vals[i]+val;
        return *this;
    }
    /**
     * The '+=' operator between matrix and matrix
     * @param a the right-hand side matrix (the dimensions should be the same)
     */
    inline MatrixXd& operator+=(const MatrixXd &a){
        MatrixXd temp(m_m,m_n);
        if(m_m==a.getM()&&m_n==a.getN()){
            for(int i=0;i<m_mn;++i) m_vals[i]=m_vals[i]+a.m_vals[i];
            return *this;
        }
        else{
            MessagePrinter::printErrorTxt("a+b can\'t be applied to two matrix with different size");
            MessagePrinter::exitAsFem();
        }
        return *this;
    }
    //****************************
    //*** for -
    /**
     * The '-' operator between matrix and scalar
     * @param val the right-hand side scalar (double type)
     */
    inline MatrixXd operator-(const double &val)const{
        MatrixXd temp(m_m,m_n);
        for(int i=0;i<m_mn;++i) temp.m_vals[i]=m_vals[i]-val;
        return temp;
    }
    /**
     * The '-' operator between matrix and matrix
     * @param a the right-hand side matrix (the dimensions should be the same)
     */
    inline MatrixXd operator-(const MatrixXd &a)const{
        MatrixXd temp(m_m,m_n);
        if(m_m==a.getM()&&m_n==a.getN()){
            for(int i=0;i<m_mn;++i) temp.m_vals[i]=m_vals[i]-a.m_vals[i];
            return temp;
        }
        else{
            MessagePrinter::printErrorTxt("a-b can\'t be applied to two matrix with different size");
            MessagePrinter::exitAsFem();
        }
        return temp;
    }
    //*** for -=
    /**
     * The '-=' operator between matrix and scalar
     * @param val the right-hand side scalar (double type)
     */
    inline MatrixXd& operator-=(const double &val){
        for(int i=0;i<m_mn;++i) m_vals[i]=m_vals[i]-val;
        return *this;
    }
    /**
     * The '-=' operator between matrix and matrix
     * @param a the right-hand side matrix (the dimensions should be the same)
     */
    inline MatrixXd& operator-=(const MatrixXd &a){
        MatrixXd temp(m_m,m_n);
        if(m_m==a.getM()&&m_n==a.getN()){
            for(int i=0;i<m_mn;++i) m_vals[i]=m_vals[i]-a.m_vals[i];
            return *this;
        }
        else{
            MessagePrinter::printErrorTxt("a-b can\'t be applied to two matrix with different size");
            MessagePrinter::exitAsFem();
        }
        return *this;
    }
    //****************************
    //*** for *
    /**
     * The '*' operator between matrix and scalar
     * @param val the right-hand side scalar (double type)
     */
    inline MatrixXd operator*(const double &val)const{
        MatrixXd temp(m_m,m_n);
        for(int i=0;i<m_mn;++i) temp.m_vals[i]=m_vals[i]*val;
        return temp;
    }
    /**
     * The '*' operator between matrix and vector
     * @param a the right-hand side vector (the dimensions should be the same)
     */
    inline VectorXd operator*(const VectorXd &a)const{
        VectorXd temp(m_m,0.0);
        if(m_n!=a.getM()){
            MessagePrinter::printErrorTxt("A*b should be applied to A matrix with the same cols as b vector!");
            MessagePrinter::exitAsFem();
        }
        else{
            for(int i=1;i<=m_m;i++){
                temp(i)=0.0;
                for(int j=1;j<=m_n;j++){
                    temp(i)+=(*this)(i,j)*a(j);
                }
            }
            return temp;
        }
        return temp;
    }
    /**
     * The '*' operator between matrix and matrix
     * @param a the right-hand side matrix (the dimensions should be the same)
     */
    inline MatrixXd operator*(const MatrixXd &a)const{
        MatrixXd temp(m_m,a.getN());
        if(m_n!=a.getM()){
            MessagePrinter::printErrorTxt("A*B should be applied to A matrix with the same cols as the rows of B matrix!");
            MessagePrinter::exitAsFem();
        }
        else{
            for(int i=1;i<=m_m;i++){
                for(int j=1;j<=a.getN();j++){
                    temp(i,j)=0.0;
                    for(int k=1;k<=a.getM();k++){
                        temp(i,j)+=(*this)(i,k)*a(k,j);
                    }
                }
            }
            return temp;
        }
        return temp;
    }
    //*** for *=
    /**
     * The '*=' operator between matrix and scalar
     * @param val the right-hand side scalar (double type)
     */
    inline MatrixXd& operator*=(const double &val){
        for(int i=0;i<m_mn;++i) m_vals[i]=m_vals[i]*val;
        return *this;
    }
    //****************************
    //*** for /
    /**
     * The '/' operator between matrix and scalar
     * @param val the right-hand side scalar (double type)
     */
    inline MatrixXd operator/(const double &val)const{
        if(abs(val)<1.0e-15){
            MessagePrinter::printErrorTxt("val="+to_string(val)+" is singular for / operator in MatrixXd");
            MessagePrinter::exitAsFem();
        }
        MatrixXd temp(m_m,m_n);
        for(int i=0;i<m_mn;++i) temp.m_vals[i]=m_vals[i]/val;
        return temp;
    }
    //*** for /=
    /**
     * The '/=' operator between matrix and scalar
     * @param val the right-hand side scalar (double type)
     */
    inline MatrixXd& operator/=(const double &val){
        if(abs(val)<1.0e-15){
            MessagePrinter::printErrorTxt("val="+to_string(val)+" is singular for / operator in MatrixXd");
            MessagePrinter::exitAsFem();
        }
        for(int i=0;i<m_mn;++i) m_vals[i]=m_vals[i]/val;
        return *this;
    }
    /**
     * This function will set the whole matrix to zero
     */
    void setToZero(){
        fill(m_vals.begin(),m_vals.end(),0.0);
    }
    /**
     * This function will set each element of the matrix to be random value
     */
    void setToRandom(){
        srand(time(0));
        for(int i=0;i<m_mn;++i) m_vals[i]=static_cast<double>(1.0*rand()/RAND_MAX);
    }
    /**
     * This function return the inverse matrix of current one, it should be noted
     * this function will not change the value of current matrix
     */
    inline MatrixXd inverse()const{
        if(m_m!=m_n){
            MessagePrinter::printErrorTxt("the inverse operation only works for square matrix");
            MessagePrinter::exitAsFem();
        }
        Eigen::MatrixXd Mat(m_m,m_n),MatInv(m_m,m_n);
        MatrixXd temp(m_m,m_n);
        for(int i=1;i<=m_m;i++){
            for(int j=1;j<=m_n;j++){
                Mat.coeffRef(i-1,j-1)=(*this)(i,j);
            }
        }
        MatInv=Mat.inverse();
        for(int i=1;i<=m_m;i++){
            for(int j=1;j<=m_n;j++){
                temp(i,j)=MatInv.coeff(i-1,j-1);
            }
        }
        return temp;
    }
    /**
     * This function return the determinant of the current matrix
     */
    inline double det()const{
        Eigen::MatrixXd Mat(m_m,m_n);
        for(int i=1;i<=m_m;i++){
            for(int j=1;j<=m_n;j++){
                Mat.coeffRef(i-1,j-1)=(*this)(i,j);
            }
        }
        return Mat.determinant();
    }
    /**
     * This function return the transporse matrix of current one,
     * the current matrix will not be changed
     */
    inline MatrixXd transpose() const{
        MatrixXd temp(m_n,m_m);
        for(int i=1;i<=m_m;i++){
            for(int j=1;j<=m_n;j++){
                temp(j,i)=(*this)(i,j);
            }
        }
        return temp;
    }
    /**
     * This function return the transposed matrix, important, the current matrix will
     * be transposed, if you don't want to use this one, then please call transpose()
     */
    inline void transposed(){
        MatrixXd temp=(*this).transpose();
        (*this).resize(temp.getM(),temp.getN());
        (*this)=temp;
    }

    /**
     * This function solve the system equation Ax=b
     * @param x input vector which stores the solution
     * @param b input vector which serves as the right hand side term
    */
    void solve(const VectorXd &b,VectorXd &x) const;
    /**
     * This function solve the system equation Ax=b, with a given b vector and return the solution x
     * @param b input vector which serves as the right hand side term
    */
    VectorXd solve(const VectorXd &b) const;

private:
    vector<double> m_vals;/**< double type vector to store the matrix element*/
    int m_m; /**< the integer variable for the 1st dimension of the matrix*/
    int m_n; /**< the integer variable for the 2nd dimension of the matrix*/
    int m_mn;/**< the integer variable for the total length of the matrix*/
};