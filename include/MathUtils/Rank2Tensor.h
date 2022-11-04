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
//+++ Date   : 2020.10.17
//+++ Update : 2022.07.24->re-write & re-design @Yang Bai
//+++ Purpose: Implement rank-2 tensor class for the common
//+++          tensor manipulation in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <cmath>

#include "MathUtils/Vector3d.h"
#include "MathUtils/Rank4Tensor.h"

#include "Utils/MessagePrinter.h"

class Rank4Tensor;

using std::fill;
using std::sqrt;
using std::abs;

/**
 * This class implement the general manipulation for rank-2 tensor.
 */
class Rank2Tensor{
public:
    /**
     * different initial method for rank-2 tensor
     */
    enum InitMethod{
        ZERO,
        IDENTITY,
        RANDOM
    };
public:
    /**
     * constructor
     */
    Rank2Tensor();
    Rank2Tensor(const double &val);
    Rank2Tensor(const Rank2Tensor &a);
    Rank2Tensor(const InitMethod &initmethod);
    ~Rank2Tensor();
    //**********************************************************************
    //*** for row, col and other elements access
    //**********************************************************************
    /**
     * get the ith row of the rank-2 tensor
     * @param i i-th row number, start from 1 to 3
     */
    inline Vector3d getIthRow(const int &i)const{
        Vector3d temp(0.0);
        temp(1)=(*this)(i,1);
        temp(2)=(*this)(i,2);
        temp(3)=(*this)(i,3);
        return temp;
    }
    /**
     * get the ith column of the rank-2 tensor
     * @param i i-th col number, start from 1 to 3
     */
    inline Vector3d getIthCol(const int &i)const{
        Vector3d temp(0.0);
        temp(1)=(*this)(1,i);
        temp(2)=(*this)(2,i);
        temp(3)=(*this)(3,i);
        return temp;
    }
    //**********************************************************************
    //*** for operator override
    //**********************************************************************
    //*******************************
    //*** for ()  operator
    //*******************************
    /** for index based access(start from 1, instead of zero !!!)
     * @param i i index of the rank-2 tensor, start from 1
     * @param j j index of the rank-2 tensor, start from 1
     */
    inline double operator()(const int &i,const int &j) const{
        if(i<1||i>3 || j<1||j>3 ){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" or j="+to_string(j)+" is out of range when you call a rank-2 tensor");
            MessagePrinter::exitAsFem();
        }
        return m_vals[(i-1)*N+j-1];
    }
    /** for index based access(start from 1, instead of zero !!!)
     * @param i i index of the rank-2 tensor, start from 1
     * @param j j index of the rank-2 tensor, start from 1
     */
    inline double& operator()(const int &i,const int &j){
        if(i<1||i>3 || j<1||j>3 ){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" or j="+to_string(j)+" is out of range when you call a rank-2 tensor");
            MessagePrinter::exitAsFem();
        }
        return m_vals[(i-1)*N+j-1];
    }
    //*******************************
    //*** for []  operator
    //*******************************
    /**
     * [] operator for element access
     * @param i global index, range from 1~9
     */
    inline double operator[](const int &i) const{
        if(i<1||i>N2){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range when you call a rank-2 tensor");
            MessagePrinter::exitAsFem();
        }
        return m_vals[i-1];
    }
    /**
     * [] operator for element access
     * @param i global index, range from 1~9
     */
    inline double& operator[](const int &i){
        if(i<1||i>N2){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range when you call a rank-2 tensor");
            MessagePrinter::exitAsFem();
        }
        return m_vals[i-1];
    }
    //*******************************
    //*** for =  operator
    //*******************************
    /**
     * '=' operator for scalar
     * @param a right hand side scalar
     */
    inline Rank2Tensor& operator=(const double &a){
        for(int i=0;i<N2;i++) m_vals[i]=a;
        return *this;
    }
    /**
     * '=' operator for rank-2 tensor
     * @param a the right hand side rank-2 tensor
     */
    inline Rank2Tensor& operator=(const Rank2Tensor &a){
        for(int i=0;i<N2;i++) m_vals[i]=a.m_vals[i];
        return *this;
    }
    //*******************************
    //*** for +  operator
    //*******************************
    /**
     * '+' operator for scalar
     * @param a right hand side scalar
     */
    inline Rank2Tensor operator+(const double &a) const{
        Rank2Tensor temp(0.0);
        for(int i=0;i<N2;i++) temp.m_vals[i]=m_vals[i]+a;
        return temp;
    }
    /**
     * '+' operator for rank-2 tensor
     * @param a right hand side rank-2 tensor
     */
    inline Rank2Tensor operator+(const Rank2Tensor &a) const{
        Rank2Tensor temp(0.0);
        for(int i=0;i<N2;i++) temp.m_vals[i]=m_vals[i]+a.m_vals[i];
        return temp;
    }
    //*******************************
    //*** for +=  operator
    //*******************************
    /**
     * '+=' operator for scalar
     * @param a right hand side scalar
     */
    inline Rank2Tensor& operator+=(const double &a) {
        for(int i=0;i<N2;i++) m_vals[i]+=a;
        return *this;
    }
    /**
     * '+=' for rank-2 tensor
     * @param a right hand side tensor
     */
    inline Rank2Tensor& operator+=(const Rank2Tensor &a){
        for(int i=0;i<N2;i++) m_vals[i]+=a.m_vals[i];
        return *this;
    }
    /**
     * add scalar value to diagnal element, the current tensor=old-tensor+I*a, where I is the identity tensor
     * @param val the scalar value to be added
     */
    inline Rank2Tensor& addIa(const double &val){
        (*this)(1,1)+=val;
        (*this)(2,2)+=val;
        (*this)(3,3)+=val;
        return *this;
    }
    //*******************************
    //*** for -  operator
    //*******************************
    /**
     * '-' operator for scalar
     * @param a right hand side scalar
     */
    inline Rank2Tensor operator-(const double &a) const{
        Rank2Tensor temp(0.0);
        for(int i=0;i<N2;i++) temp.m_vals[i]=m_vals[i]-a;
        return temp;
    }
    /**
     * '-' operator for rank-2 tensor
     * @param a right hand side rank-2 tensor
     */
    inline Rank2Tensor operator-(const Rank2Tensor &a) const{
        Rank2Tensor temp(0.0);
        for(int i=0;i<N2;i++) temp.m_vals[i]=m_vals[i]-a.m_vals[i];
        return temp;
    }
    //*******************************
    //*** for -=  operator
    //*******************************
    /**
     * '-=' operator for scalar
     * @param a right hand side scalar
     */
    inline Rank2Tensor& operator-=(const double &a) {
        for(int i=0;i<N2;i++) m_vals[i]-=a;
        return *this;
    }
    /**
     * '-=' for rank-2 tensor
     * @param a right hand side tensor
     */
    inline Rank2Tensor& operator-=(const Rank2Tensor &a){
        for(int i=0;i<N2;i++) m_vals[i]-=a.m_vals[i];
        return *this;
    }
    //*******************************
    //*** for *  operator
    //*******************************
    /**
     * '*' for scalar
     * @param a right hand side scalar
     */
    inline Rank2Tensor operator*(const double &a) const{
        Rank2Tensor temp(0.0);
        for(int i=0;i<N2;++i) temp.m_vals[i]=m_vals[i]*a;
        return temp;
    }
    /**
     * '*' operator for vector3d
     * @param a right hand side vector3d
     */
    inline Vector3d operator*(const Vector3d &a) const{
        Vector3d temp(0.0);
        for(int i=1;i<=N;i++){
            temp(i)=(*this)(i,1)*a(1)+(*this)(i,2)*a(2)+(*this)(i,3)*a(3);
        }
        return temp;
    }
    /**
     * '*' operator for rank-2 tensor, return \f$a_{ij}=b_{ik}c_{kj}\f$
     * @param a the right hand side rank-2 tensor
     */
    inline Rank2Tensor operator*(const Rank2Tensor &a) const{
        // return A*B(still rank-2 tensor)
        Rank2Tensor temp(0.0);
        for(int i=1;i<=N;i++){
            for(int j=1;j<=N;j++){
                temp(i,j)=(*this)(i,1)*a(1,j)+(*this)(i,2)*a(2,j)+(*this)(i,3)*a(3,j);
            }
        }
        return temp;
    }
    /**
     * double dot(:, the double contruction) between two rank-2 tensor, the result is \f$\sum a_{ij}b_{ij}\f$, which is a scalar
     * @param a the right hand side rank-2 tensor
     */
    inline double doubledot(const Rank2Tensor &a) const{
        // return A:B calculation
        double sum=0.0;
        for(int i=1;i<=N;i++){
            for(int j=1;j<=N;j++){
                // You may see A:B=A_ijB_ji in other books/literature, here we use A_ijB_ij
                // in Rank4Tensor, we follow the same definition!
                sum+=(*this)(i,j)*a(i,j);// use this to get the positive definite case!!!
            }
        }
        return sum;
    }
    /**
     * double dot(:, the double contruction) between rank-2 and rank-4 tensor, 
     * the result is \f$\sum a_{ij}c_{ijkl}=b_{kl}\f$, which is a scalar
     * @param a the right hand side rank-2 tensor
     */
    Rank2Tensor doubledot(const Rank4Tensor &a) const;
    //*******************************
    //*** for *=  operator
    //*******************************
    /**
     * '*=' operator for scalar
     * @param a right hand side scalar
     */
    inline Rank2Tensor& operator*=(const double &a) {
        for(int i=0;i<N2;i++) m_vals[i]*=a;
        return *this;
    }
    /**
     * '*=' for rank-2 tensor
     * @param a right hand side rank-2 tensor
     */
    inline Rank2Tensor& operator*=(const Rank2Tensor &a){
        Rank2Tensor temp(0.0);
        temp=(*this)*a;
        (*this)=temp;
        return *this;
    }
    //*******************************
    //*** for /  operator
    //*******************************
    /**
     * '/' for scalar
     * @param a right hand side scalar
     */
    inline Rank2Tensor operator/(const double &a) const{
        if(abs(a)<1.0e-16){
            MessagePrinter::printErrorTxt("a="+to_string(a)+" is singular for / operator in rank-2 tensor");
            MessagePrinter::exitAsFem();
        }
        Rank2Tensor temp(0.0);
        for(int i=0;i<N2;i++) temp.m_vals[i]=m_vals[i]/a;
        return temp;
    }
    //*******************************
    //*** for /=  operator
    //*******************************
    /**
     * '/=' operator for scalar
     * @param a right hand side scalar
     */
    inline Rank2Tensor& operator/=(const double &a){
        if(abs(a)<1.0e-16){
            MessagePrinter::printErrorTxt("a="+to_string(a)+" is singular for /= operator in rank-2 tensor");
            MessagePrinter::exitAsFem();
        }
        for(int i=0;i<N2;i++) m_vals[i]/=a;
        return *this;
    }
    //*******************************************************************
    //*** for left hand side manipulation
    //*******************************************************************
    /**
     * '*' for left hand side scalar
     * @param lhs left hand side scalar value
     * @param a right hand side rank-2 tensor
     */
    friend Rank2Tensor operator*(const double &lhs,const Rank2Tensor &a);
    //*** for left hand vector times rank-2 tensor
    /**
     * '*' operator for left hand side vector3d
     * @param lhs the left hand side vector3d
     * @param a the right hand side rank-2 tensor
     */
    friend Vector3d operator*(const Vector3d &lhs,const Rank2Tensor &a);
    //*******************************************************************
    //*** for advanced math operators
    //*******************************************************************
    /**
     * get the exponetial formula of current rank-2 tensor
    */
    Rank2Tensor exp()const{
        Rank2Tensor I;
        I.setToIdentity();
        return I
              +(*this)
              +(*this)*(*this)*(1.0/(1.0*2.0))
              +(*this)*(*this)*(*this)*(1.0/(1.0*2.0*3.0))
              +(*this)*(*this)*(*this)*(*this)*(1.0/(1.0*2.0*3.0*4.0))
              +(*this)*(*this)*(*this)*(*this)*(*this)*(1.0/(1.0*2.0*3.0*4.0*5.0))
              +(*this)*(*this)*(*this)*(*this)*(*this)*(*this)*(1.0/(1.0*2.0*3.0*4.0*5.0*6.0));
    }
    /**
     * get the exponential of input rank-2 tensor
     * @param a the given rank-2 tensor
    */
    friend Rank2Tensor exp(const Rank2Tensor &a);
    /**
     * get the input rank-2 tensor's (a*Rank-2) exponential's first order derivative w.r.t. scalar a
     * @param a the scalar factor
     * @param b the rank-2 tensor
    */
    friend Rank2Tensor dexp(const double &a,const Rank2Tensor &b);
    //*******************************************************************
    //*** for higher order tensor calculation
    //*******************************************************************
    /**
     * Otime(\f$\otimes\f$) between two rank-2 tensor, which will return \f$c_{ijkl}=a_{ij}b_{kl}\f$
     * @param a right hand side rank-2 tensor
     */
    Rank4Tensor otimes(const Rank2Tensor &a) const;
    /**
     * Odot (\f$\odot\f$) between two rank-2 tensor, which will return \f$c_{ijkl}=\frac{1}{2}(a_{ik}b_{jl}+a_{il}b_{jk})\f$
     * @param a right hand side rank-2 tensor
     */
    Rank4Tensor odot(const Rank2Tensor &a) const;

    /**
     * IJlk times (\f$IJ\otimes lk\f$) operator
     * @param a right hand side rank-2 tensor
     */
    Rank4Tensor ijXlk(const Rank2Tensor &a) const;

    /**
     * IkJl times (\f$Ik\otimes Jl\f$)
     * @param a right hand side rank-2 tensor
     */
    Rank4Tensor ikXjl(const Rank2Tensor &a) const;
    /**
     * IklJ times (\f$Ik\otimes lJ\f$)
     * @param a right hand side rank-2 tensor
     */
    Rank4Tensor ikXlj(const Rank2Tensor &a) const;

    /**
     * IlJk times (\f$Il\otimes Jk\f$)
     * @param a right hand side rank-2 tensor
     */
    Rank4Tensor ilXjk(const Rank2Tensor &a) const;
    /**
     * IlkJ times (\f$Il\otimes kJ\f$)
     * @param a right hand side rank-2 tensor
     */
    Rank4Tensor ilXkj(const Rank2Tensor &a) const;

    /**
     * JiKl times (\f$Ji\otimes Kl\f$)
     * @param a right hand side rank-2 tensor
     */
    Rank4Tensor jiXkl(const Rank2Tensor &a) const;
    /**
     * Jilk times (\f$Ji\otimes lK\f$)
     * @param a right hand side rank-2 tensor
     */
    Rank4Tensor jiXlk(const Rank2Tensor &a) const;

    /**
     * JkIl times (\f$Jk\otimes Il\f$)
     * @param a right hand side rank-2 tensor
     */
    Rank4Tensor jkXil(const Rank2Tensor &a) const;
    /**
     * JklI times (\f$Jk\otimes lI\f$)
     * @param a right hand side rank-2 tensor
     */
    Rank4Tensor jkXli(const Rank2Tensor &a) const;

    /**
     * JlIk times (\f$Jl\otimes Ik\f$)
     * @param a right hand side rank-2 tensor
     */
    Rank4Tensor jlXik(const Rank2Tensor &a) const;
    /**
     * JlKi times (\f$Jl\otimes Ki\f$)
     * @param a right hand side rank-2 tensor
     */
    Rank4Tensor jlXki(const Rank2Tensor &a) const;
    
    /**
     * kiJL times (\f$ki\otimes JL\f$)
     * @param a right hand side rank-2 tensor
     */
    Rank4Tensor kiXjl(const Rank2Tensor &a) const;
    /**
     * kiLJ times (\f$ki\otimes LJ\f$)
     * @param a right hand side rank-2 tensor
     */
    Rank4Tensor kiXlj(const Rank2Tensor &a) const;

    /**
     * kjil times (\f$kj\otimes IL\f$)
     * @param a right hand side rank-2 tensor
     */
    Rank4Tensor kjXil(const Rank2Tensor &a) const;
    /**
     * kjLI times (\f$kj\otimes LI\f$)
     * @param a right hand side rank-2 tensor
     */
    Rank4Tensor kjXli(const Rank2Tensor &a) const;
    
    /**
     * klIJ times (\f$kl\otimes IJ\f$)
     * @param a right hand side rank-2 tensor
     */
    Rank4Tensor klXij(const Rank2Tensor &a) const;
    /**
     * klJI times (\f$kl\otimes JI\f$)
     * @param a right hand side rank-2 tensor
     */
    Rank4Tensor klXji(const Rank2Tensor &a) const;

    /**
     * liKJ times (\f$li\otimes KJ\f$)
     * @param a right hand side rank-2 tensor
     */
    Rank4Tensor liXkj(const Rank2Tensor &a) const;
    /**
     * liJK times (\f$li\otimes JK\f$)
     * @param a right hand side rank-2 tensor
     */
    Rank4Tensor liXjk(const Rank2Tensor &a) const;

    /**
     * ljIK times (\f$lj\otimes IK\f$)
     * @param a right hand side rank-2 tensor
     */
    Rank4Tensor ljXik(const Rank2Tensor &a) const;
    /**
     * ljKI times (\f$lj\otimes KI\f$)
     * @param a right hand side rank-2 tensor
     */
    Rank4Tensor ljXki(const Rank2Tensor &a) const;

    /**
     * lkIJ times (\f$lk\otimes IJ\f$)
     * @param a right hand side rank-2 tensor
     */
    Rank4Tensor lkXij(const Rank2Tensor &a) const;
    /**
     * lkJI times (\f$lk\otimes JI\f$)
     * @param a right hand side rank-2 tensor
     */
    Rank4Tensor lkXji(const Rank2Tensor &a) const;
    //**********************************************************************
    //*** for general settings
    //**********************************************************************
    /**
     * set the elements of current rank-2 tensor to be zero
     */
    inline void setToZeros(){
        for(int i=0;i<N2;i++) m_vals[i]=0.0;
    }
    /**
     * set current rannk-2 tensor to be an identitiy tensor, where \f$a_{ij}=\delta_{ij}\f$.
     */
    inline void setToIdentity(){
        for(int i=1;i<=N;++i){
            for(int j=1;j<=N;++j){
                if(i==j){
                    (*this)(i,j)=1.0;
                }
                else{
                    (*this)(i,j)=0.0;
                }
            }
        }
    }
    /**
     * set current rank-2 tensor to be a random one
     */
    inline void setToRandom(){
        srand(time(0));
        for(int i=1;i<=N;++i){
            for(int j=1;j<=N;++j){
                (*this)(i,j)=static_cast<double>(1.0*rand()/RAND_MAX);
            }
        }
    }
    /**
     * cross dot \f$\otimes\f$ for two vector(from STL), set current one to \f$c_{ij}=a_{i}\times b_{j}\f$.
     * @param a vector<double> for 1st dimension
     * @param b vector<double> for 2nd dimension
     */
    inline void setFromVectorDyad(const vector<double> &a,const vector<double> &b){
        for(int i=1;i<=N;i++){
            for(int j=1;j<=N;j++){
                (*this)(i,j)=a[i-1]*b[j-1];
            }
        }
    }
    /**
     * fill up current rank-2 tensor from the dyad of two vector
     * @param a the first vector
     * @param b the second vector
     */ 
    inline void setFromVectorDyad(const Vector3d &a,const Vector3d &b){
        for(int i=1;i<=N;i++){
            for(int j=1;j<=N;j++){
                (*this)(i,j)=a(i)*b(j);
            }
        }
    }
    /**
     * fill up current rank-2 tensor from the displacement gradient
     * @param gradUx the gradient of Ux, namely, \f$\nabla u_{x}\f$
     */ 
    void setFromGradU(const Vector3d &gradUx);
    /**
     * fill up current rank-2 tensor from the gradient of ux and uy
     * @param gradUx the gradient of Ux, namely, \f$\nabla u_{x}\f$
     * @param gradUy the gradient of Uy, namely, \f$\nabla u_{y}\f$
     */
    void setFromGradU(const Vector3d &gradUx,const Vector3d &gradUy);
    /**
     * fill up current rank-2 tensor from the gradient of ux and uy
     * @param gradUx the gradient of Ux, namely, \f$\nabla u_{x}\f$
     * @param gradUy the gradient of Uy, namely, \f$\nabla u_{y}\f$
     */
    void setFromGradU(const Vector3d &gradUx,const Vector3d &gradUy,const Vector3d &gradUz);
    /**
     * fill up current rank-2 tensor with eular angle, which should be used for the rotation tensor
     * @param theta1 the first eular angle, the angle must be degree
     * @param theta2 the second eular angle, the angle must be degree
     * @param theta3 the third eular angle, the angle must be degree
     */
    void setRotationTensorFromEulerAngle(const double &theta1,const double &theta2,const double &theta3);

    //**********************************************************************
    //*** for some common mathematic manipulations
    //**********************************************************************
    /**
     * return the trace of a rank-2 tensor, result is \f$\sum a_{ii}\f$
     */
    inline double trace() const{
        return (*this)(1,1)+(*this)(2,2)+(*this)(3,3);
    }
    /**
     * return the determinant of the current rank-2 tensor
     */
    inline double det() const{
        // taken from http://mathworld.wolfram.com/Determinant.html
        // Eq.10
        return (*this)(1,1)*(*this)(2,2)*(*this)(3,3)
              -(*this)(1,1)*(*this)(2,3)*(*this)(3,2)
              -(*this)(1,2)*(*this)(2,1)*(*this)(3,3)
              +(*this)(1,2)*(*this)(2,3)*(*this)(3,1)
              +(*this)(1,3)*(*this)(2,1)*(*this)(3,2)
              -(*this)(1,3)*(*this)(2,2)*(*this)(3,1);
    }
    /**
     * return the \f$L_{2}\f$ norm of current rank-2 tensor, result is \f$\sqrt{\sum a_{ij}^{2}}\f$
     */
    inline double norm() const{
        double sum=0.0;
        for(int i=0;i<N2;i++) sum+=static_cast<double>(m_vals[i]*m_vals[i]);
        return sqrt(sum);
    }
    
    /**
     * return the \f$L_{2}\f$ norm^2 of current rank-2 tensor, result is \f$\sqrt{\sum a_{ij}^{2}}\f$
     */
    inline double normsq() const{
        double sum=0.0;
        for(int i=0;i<N2;i++) sum+=static_cast<double>(m_vals[i]*m_vals[i]);
        return sum;
    }
    //*** for the different invariants of stress(strain)
    /**
     * return the first invariant of current rank-2 tensor, namely, \f$I_{1}\f$.
     */
    inline double firstInvariant() const{
        return trace();
    }
    /**
     * return the second invariant of current tensor, namely, \f$I_{2}\f$.
     */
    inline double secondInvariant() const{
        double trAA=((*this)*(*this)).trace();
        double trA=trace();
        return 0.5*(trA*trA-trAA);
    }
    /**
     * return the third invariant of current tensor, namely, \f$I_{3}\f$.
     */
    inline double thirdInvariant() const{
        return det();
    }
    //*** for inverse
    /**
     * return the inverse tensor of current one, namely, \f$\mathbf{A}^{-1}\f$.
     */
    inline Rank2Tensor inverse() const{
        double J=det();
        if(abs(J)<1.0e-16){
            MessagePrinter::printErrorTxt("inverse operation failed for a singular rank-2 tensor !");
            MessagePrinter::exitAsFem();
        }
        Rank2Tensor inv(0.0);
        // taken from wiki:
        //   https://en.wikipedia.org/wiki/Invertible_matrix
        double A= (*this)(2,2)*(*this)(3,3)-(*this)(2,3)*(*this)(3,2);
        double D=-(*this)(1,2)*(*this)(3,3)+(*this)(1,3)*(*this)(3,2);
        double G= (*this)(1,2)*(*this)(2,3)-(*this)(1,3)*(*this)(2,2);
        inv(1,1)=A/J;inv(1,2)=D/J;inv(1,3)=G/J;

        double B=-(*this)(2,1)*(*this)(3,3)+(*this)(2,3)*(*this)(3,1);
        double E= (*this)(1,1)*(*this)(3,3)-(*this)(1,3)*(*this)(3,1);
        double H=-(*this)(1,1)*(*this)(2,3)+(*this)(1,3)*(*this)(2,1);
        inv(2,1)=B/J;inv(2,2)=E/J;inv(2,3)=H/J;

        double C= (*this)(2,1)*(*this)(3,2)-(*this)(2,2)*(*this)(3,1);
        double F=-(*this)(1,1)*(*this)(3,2)+(*this)(1,2)*(*this)(3,1);
        double I= (*this)(1,1)*(*this)(2,2)-(*this)(1,2)*(*this)(2,1);
        inv(3,1)=C/J;inv(3,2)=F/J;inv(3,3)=I/J;
        return inv;
    }
    /**
     * return the transpose tensor of current one, namely, \f$\mathbf{A}^{T}\f$.
     */
    inline Rank2Tensor transpose() const{
        Rank2Tensor temp(0.0);
        for(int i=1;i<=N;i++){
            for(int j=1;j<=N;j++){
                temp(i,j)=(*this)(j,i);
            }
        }
        return temp;
    }
    /**
     * transpose current tensor, and overwrite the original one
     */
    inline void transposed(){
        Rank2Tensor temp(0.0);
        temp=(*this).transpose();
        (*this)=temp;
    }
    /**
     * return the deviatoric part of current rank-2 tensor
     */
    inline Rank2Tensor dev()const{
        Rank2Tensor temp(0.0),I(0.0);
        I.setToIdentity();
        temp=(*this)-I*(this->trace()/3.0);
        return temp;
    }
    //**********************************************************************
    //*** for decomposition
    //**********************************************************************
    /**
     * calculate the eigen value and eigen vector for current rank-2 tensor
     * @param eigval the double array, which stores the eigen value
     * @param eigvec the rank-2 tensor, where each column store the related eigen vector
     */
    void calcEigenValueAndEigenVectors(double (&eigval)[3],Rank2Tensor &eigvec) const;
    /**
     * calculate the positive projection tensor(a rank-4 tensor), this algorithm is taken from Miehe's paper, for the details, please see the cpp file
     * @param eigval the double array, which stores the eigen value
     * @param eigvec the rank-2 tensor, whoses column stores the related eigen vector
     */
    Rank4Tensor calcPositiveProjTensor(double (&eigval)[3],Rank2Tensor &eigvec) const;
    /**
     * calculate the positive projection tensor (rank-4 tensor) based on current rank-2 tensor
     */
    Rank4Tensor getPositiveProjectionTensor() const;
    //**********************************************************************
    //*** for decomposition
    //**********************************************************************
    /**
     * print all the elements to the terminal
     */
    inline void print() const{
        PetscPrintf(PETSC_COMM_WORLD,"*** %14.6e ,%14.6e ,%14.6e ***\n",(*this)(1,1),(*this)(1,2),(*this)(1,3));
        PetscPrintf(PETSC_COMM_WORLD,"*** %14.6e ,%14.6e ,%14.6e ***\n",(*this)(2,1),(*this)(2,2),(*this)(2,3));
        PetscPrintf(PETSC_COMM_WORLD,"*** %14.6e ,%14.6e ,%14.6e ***\n",(*this)(3,1),(*this)(3,2),(*this)(3,3));
    }



private:
    const int N=3;/**< the dimension of current rank-2 tensor */
    const int N2=9;/**< the total length of current rank-2 tensor */
    vector<double> m_vals;/**< vector for the elements of rank-2 tensor */
};