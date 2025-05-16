//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.10.17
//+++ Update : 2022.07.24->re-write & re-design @Yang Bai
//+++ Purpose: Implement rank-4 tensor class for the common
//+++          tensor manipulation in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <cmath>

#include "MathUtils/Vector3d.h"
#include "MathUtils/Rank2Tensor.h"

#include "Utils/MessagePrinter.h"

class Rank2Tensor;

using std::fill;
using std::sqrt;
using std::abs;

/**
 * This class implement the general manipulation for rank-2 tensor.
 */
class Rank4Tensor{
public:
    /**
     * enum for different initializing methods
     */
    enum InitMethod{
        ZERO,
        IDENTITY,
        IDENTITY4,
        IDENTITY4TRANS,
        IDENTITY4SYMMETRIC,
        RANDOM
    };
public:
    /**
     * constructor
     */
    Rank4Tensor();
    Rank4Tensor(const double &val);
    Rank4Tensor(const Rank4Tensor &a);
    Rank4Tensor(const InitMethod &method);
    ~Rank4Tensor();

    //**********************************************************************
    //*** for operator override
    //**********************************************************************
    //*******************************
    //*** for ()  operator
    //*******************************
    /**
     * () operator for rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param i i index, start from 1 instead of 0 !!!
     * @param j j index, start from 1 instead of 0 !!!
     * @param k k index, start from 1 instead of 0 !!!
     * @param l l index, start from 1 instead of 0 !!!
     */
    inline double operator()(const int &i,const int &j,const int &k,const int &l) const{
        if(i<1||i>3 || j<1||j>3 || k<1||k>3 || l<1||l>3){
            MessagePrinter::printErrorTxt("your i or j or k or l is out of range when you access a rank-4 tensor");
            MessagePrinter::exitAsFem();
        }
        return m_vals[i-1][j-1][k-1][l-1];
    }
    /**
     * () operator for rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param i i index, start from 1 instead of 0 !!!
     * @param j j index, start from 1 instead of 0 !!!
     * @param k k index, start from 1 instead of 0 !!!
     * @param l l index, start from 1 instead of 0 !!!
     */
    inline double& operator()(const int &i,const int &j,const int &k,const int &l){
        if(i<1||i>3 || j<1||j>3 || k<1||k>3 || l<1||l>3){
            MessagePrinter::printErrorTxt("your i or j or k or l is out of range when you access a rank-4 tensor");
            MessagePrinter::exitAsFem();
        }
        return m_vals[i-1][j-1][k-1][l-1];
    }
    /**
     * get the voigt component of current rank-4 tensor
     * @param i i-index for 1st dimension
     * @param j j-index for 2nd dimension
     */
    double getVoigtComponent(const int &i,const int &j)const;

    /**
     * return the reference via the voigt component of current rank-4 tensor
     * @param i i-index for 1st dimension
     * @param j j-index for 2nd dimension
     */
    double& voigtComponent(const int &i,const int &j);
    /**
     * return the C_iJkL*N,J*N,L value for jacobian matrix calculation
     * @param i the 1st dimension index
     * @param j the 2nd dimension index
     * @param grad_test the test shape function's gradient
     * @param grad_trial the trial shape function's gradient
     */
    inline double getIKComponent(const int &i,const int &k,const Vector3d &grad_test,const Vector3d &grad_trial)const{
        if(i<1||i>3){
            MessagePrinter::printErrorTxt("your i(="+to_string(i)+") is out of range when you call getIKComponent");
            MessagePrinter::exitAsFem();
        }
        if(k<1||k>3){
            MessagePrinter::printErrorTxt("your k(="+to_string(k)+") is out of range when you call getIKComponent");
            MessagePrinter::exitAsFem();
        }
        return ( (*this)(i,1,k,1)*grad_test(1)
                +(*this)(i,2,k,1)*grad_test(2)
                +(*this)(i,3,k,1)*grad_test(3))*grad_trial(1)
              +( (*this)(i,1,k,2)*grad_test(1)
                +(*this)(i,2,k,2)*grad_test(2)
                +(*this)(i,3,k,2)*grad_test(3))*grad_trial(2)
              +( (*this)(i,1,k,3)*grad_test(1)
                +(*this)(i,2,k,3)*grad_test(2)
                +(*this)(i,3,k,3)*grad_test(3))*grad_trial(3);
    }
    //*******************************
    //*** for =  operator
    //*******************************
    /**
     * = operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side scalar
     */
    inline Rank4Tensor& operator=(const double &a){
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                for(int k=0;k<3;k++){
                    for(int l=0;l<3;l++){
                        m_vals[i][j][k][l]=a;
                    }
                }
            }
        }
        return *this;
    }
    /**
     * = operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand rank-4 tensor
     */
    inline Rank4Tensor& operator=(const Rank4Tensor &a){
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                for(int k=0;k<3;k++){
                    for(int l=0;l<3;l++){
                        m_vals[i][j][k][l]=a.m_vals[i][j][k][l];
                    }
                }
            }
        }
        return *this;
    }
    //*******************************
    //*** for + operator
    //*******************************
    /**
     * + operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side scalar
     */
    inline Rank4Tensor operator+(const double &a) const{
        Rank4Tensor temp(0.0);
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                for(int k=0;k<3;k++){
                    for(int l=0;l<3;l++){
                        temp.m_vals[i][j][k][l]=m_vals[i][j][k][l]+a;
                    }
                }
            }
        }
        return temp;
    }
    /**
     * + operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand rank-4 tensor
     */
    inline Rank4Tensor operator+(const Rank4Tensor &a) const{
        Rank4Tensor temp(0.0);
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                for(int k=0;k<3;k++){
                    for(int l=0;l<3;l++){
                        temp.m_vals[i][j][k][l]=m_vals[i][j][k][l]+a.m_vals[i][j][k][l];
                    }
                }
            }
        }
        return temp;
    }
    //*******************************
    //*** for += operator
    //*******************************
    /**
     * += operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side scalar
     */
    inline Rank4Tensor& operator+=(const double &a){
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                for(int k=0;k<3;k++){
                    for(int l=0;l<3;l++){
                        m_vals[i][j][k][l]=m_vals[i][j][k][l]+a;
                    }
                }
            }
        }
        return *this;
    }
    /**
     * += operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side rank-4 tensor
     */
    inline Rank4Tensor& operator+=(const Rank4Tensor &a){
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                for(int k=0;k<3;k++){
                    for(int l=0;l<3;l++){
                        m_vals[i][j][k][l]=m_vals[i][j][k][l]+a.m_vals[i][j][k][l];
                    }
                }
            }
        }
        return *this;
    }
    //*******************************
    //*** for - operator
    //*******************************
    /**
     * - operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side scalar
     */
    inline Rank4Tensor operator-(const double &a) const{
        Rank4Tensor temp(0.0);
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                for(int k=0;k<3;k++){
                    for(int l=0;l<3;l++){
                        temp.m_vals[i][j][k][l]=m_vals[i][j][k][l]-a;
                    }
                }
            }
        }
        return temp;
    }
    /**
     * - operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand rank-4 tensor
     */
    inline Rank4Tensor operator-(const Rank4Tensor &a) const{
        Rank4Tensor temp(0.0);
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                for(int k=0;k<3;k++){
                    for(int l=0;l<3;l++){
                        temp.m_vals[i][j][k][l]=m_vals[i][j][k][l]-a.m_vals[i][j][k][l];
                    }
                }
            }
        }
        return temp;
    }
    //*******************************
    //*** for -= operator
    //*******************************
    /**
     * -= operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side scalar
     */
    inline Rank4Tensor& operator-=(const double &a){
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                for(int k=0;k<3;k++){
                    for(int l=0;l<3;l++){
                        m_vals[i][j][k][l]=m_vals[i][j][k][l]-a;
                    }
                }
            }
        }
        return *this;
    }
    /**
     * -= operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side rank-4 tensor
     */
    inline Rank4Tensor& operator-=(const Rank4Tensor &a){
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                for(int k=0;k<3;k++){
                    for(int l=0;l<3;l++){
                        m_vals[i][j][k][l]=m_vals[i][j][k][l]-a.m_vals[i][j][k][l];
                    }
                }
            }
        }
        return *this;
    }
    //*******************************
    //*** for * operator
    //*******************************
    /**
     * * operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side scalar
     */
    inline Rank4Tensor operator*(const double &a) const{
        Rank4Tensor temp(0.0);
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                for(int k=0;k<3;k++){
                    for(int l=0;l<3;l++){
                        temp.m_vals[i][j][k][l]=m_vals[i][j][k][l]*a;
                    }
                }
            }
        }
        return temp;
    }
    /**
     * * operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}=\mathbb{C}_{ijkm}\mathbf{A}_{ml}\f$
     * @param a right hand side ran-2 tensor
     */
    Rank4Tensor operator*(const Rank2Tensor &a) const;
    /**
     * double dot : between a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$ and rank-2 tensor
     * @param a right hand side rank-2 tensor
     */
    Rank2Tensor doubledot(const Rank2Tensor &a) const;
    /**
     * double dot : between a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$ and another rank-4 tensor
     * return \f$\mathbb{C}_{ijkl}=\mathbb{A}_{ijmn}\mathbb{B}_{mnkl}\f$
     * @param a right hand side rank-4 tensor
     */
    Rank4Tensor doubledot(const Rank4Tensor &a) const;
    //*******************************
    //*** for *= operator
    //*******************************
    /**
     * *= operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side scalar
     */
    inline Rank4Tensor& operator*=(const double &a){
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                for(int k=0;k<3;k++){
                    for(int l=0;l<3;l++){
                        m_vals[i][j][k][l]=m_vals[i][j][k][l]*a;
                    }
                }
            }
        }
        return *this;
    }
    //*******************************
    //*** for / operator
    //*******************************
    /**
     * / operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side scalar
     */
    inline Rank4Tensor operator/(const double &a) const{
        if(abs(a)<1.0e-16){
            MessagePrinter::printErrorTxt("a="+to_string(a)+" is singular for / operator in rank-4 tensor");
            MessagePrinter::exitAsFem();
        }
        Rank4Tensor temp(0.0);
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                for(int k=0;k<3;k++){
                    for(int l=0;l<3;l++){
                        temp.m_vals[i][j][k][l]=m_vals[i][j][k][l]/a;
                    }
                }
            }
        }
        return temp;
    }
    //*******************************
    //*** for /= operator
    //*******************************
    /**
     * /= operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side scalar
     */
    inline Rank4Tensor& operator/=(const double &a){
        if(abs(a)<1.0e-16){
            MessagePrinter::printErrorTxt("a="+to_string(a)+" is singular for /= operator in rank-4 tensor");
            MessagePrinter::exitAsFem();
        }
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                for(int k=0;k<3;k++){
                    for(int l=0;l<3;l++){
                        m_vals[i][j][k][l]=m_vals[i][j][k][l]/a;
                    }
                }
            }
        }
        return *this;
    }
    //*******************************************************************
    //*** for left hand side manipulation
    //*******************************************************************
    /**
     * * operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param lhs left hand side scalar
     * @param a right hand side rank-4 tensor
     */
    friend Rank4Tensor operator*(const double &lhs,const Rank4Tensor &a);
    /**
     * * operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param lhs left hand side rank-2 tensor
     * @param a right hand side rank-4 tensor
     */
    friend Rank4Tensor operator*(const Rank2Tensor &lhs,const Rank4Tensor &a);
    //**********************************************************************
    //*** for general settings
    //**********************************************************************
    /**
     * set current rank-4 tensor to 0
     */
    inline void setToZeros(){
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                for(int k=0;k<3;k++){
                    for(int l=0;l<3;l++){
                        m_vals[i][j][k][l]=0.0;
                    }
                }
            }
        }
    }
    /**
     * set current rank-4 tensor to idendity
     */
    inline void setToIdentity(){
        setToZeros();
        for(int i=1;i<=3;i++) (*this)(i,i,i,i)=1.0;
    }
    /**
     * set current rank-4 tensor to rank-4 identity
     */
    inline void setToIdentity4(){
        // maps a rank2 tensor to itself(no symmetric consideriation here),i.e. Iden4:rank2=rank2
        setToZeros();
        for(int i=1;i<=3;i++){
            for(int j=1;j<=3;j++){
                for(int k=1;k<=3;k++){
                    for(int l=1;l<=3;l++){
                        (*this)(i,j,k,l)=1.0*((i==k)&&(j==l));
                    }
                }
            }
        }
    }
    /**
     * set current rank-4 tensor to rank-4 transposed identity
     */
    inline void setIdentity4Transpose(){
        // maps a rank-2 tensor to its transpose, A^{T}=I4:A
        setToZeros();
        for(int i=1;i<=3;i++){
            for(int j=1;j<=3;j++){
                for(int k=1;k<=3;k++){
                    for(int l=1;l<=3;l++){
                        (*this)(i,j,k,l)=1.0*((j==k)&&(i==l));
                    }
                }
            }
        }
    }
    /**
     * set current rank-4 tensor to symmetric identity rank-4 tensor
     */
    inline void setToIdentity4Symmetric(){
        // symmetric fourth-order tensor
        setToZeros();
        for(int i=1;i<=3;i++){
            for(int j=1;j<=3;j++){
                for(int k=1;k<=3;k++){
                    for(int l=1;l<=3;l++){
                        (*this)(i,j,k,l)=0.5*((i==k)&&(j==l))
                                        +0.5*((i==l)&&(j==k));
                    }
                }
            }
        }
    }
    /**
     * set current rank-4 tensor to random values
     */
    inline void setToRandom(){
        srand(time(0));
        for(int i=1;i<=3;i++){
            for(int j=1;j<=3;j++){
                for(int k=1;k<=3;k++){
                    for(int l=1;l<=3;l++){
                        (*this)(i,j,k,l)=static_cast<double>(1.0*rand()/RAND_MAX);
                    }
                }
            }
        }
    }
    /**
     * set rank-4 tensor from lamme constant and shear modulus
     * @param Lame first Lame constant
     * @param G shear modulus
     */
    void setFromLameAndG(const double &Lame,const double &G);
    /**
     * set rank-4 tensor from Youngs modulus and poisson ratio
     * @param E Youngs modulus
     * @param Nu Poisson ratio
     */
    void setFromEAndNu(const double &E,const double &Nu);
    /**
     * set rank-4 tensor from bulk modulus and shear modulus
     * @param E Bulk modulus
     * @param Nu Shear modulus
     */
    void setFromKAndG(const double &K,const double &G);
    /**
     * set rank-4 tensor from 9-component
     * @param vec 9-component for rank-4 tensor
     */
    void setFromSymmetric9(const vector<double> &vec);
    /**
     * set rank-4 tensor to Orthotropic
     * @param vec 9-component for rank-4 tensor
     */
    void setToOrthotropic(const vector<double> &vec);
    //**********************************************************************
    //*** for advanced mathematic manipulation
    //**********************************************************************
    /**
     * rotate the current a rank-4 tensor by a rank-2 rotation tensor, this will return a new and rotated 
     * rank-4 tensor, the original one will not be changed
     * @param rotate right hand side rank-2 tensor
     */
    Rank4Tensor rotate(const Rank2Tensor &rotate) const;
    /**
     * rotate the current a rank-4 tensor by a rank-2 rotation tensor, current tensor value will be changed
     * @param rot right hand side rank-2 tensor
     */
    void rotated(const Rank2Tensor &rot);
    /**
     * pushforward current rank-4 tensor by a rank-2 rotation tensor F
     * @param F right hand side rank-2 tensor
     */
    Rank4Tensor pushForward(const Rank2Tensor &F) const;
    /**
     * conjugate pushforward current rank-4 tensor by a rank-2 rotation tensor F,
     * the conjugate tensor is given by:
     * C_ijkl=F_im*C_mjnl*F_kn
     * @param F right hand side rank-2 tensor
     */
    Rank4Tensor conjPushForward(const Rank2Tensor &F) const;

private:
    double m_vals[3][3][3][3];/**< the tensor's matrix */

};