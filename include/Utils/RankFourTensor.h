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
//+++ Date   : 2020.10.17
//+++ Purpose: Implement rank-4 tensor class for some common
//+++          tensor calculation for solid-mechanics problem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#pragma once

#include <iostream>
#include <iomanip>
#include <cmath>

#include "petsc.h"

//****************************
#include "Utils/MessagePrinter.h"
#include "Utils/RankTwoTensor.h"

using namespace std;

/**
 * forward declaration
 */
class RankTwoTensor;

/**
 * this class implement the general operators/calculations in Rank-4 tensor
 */
class RankFourTensor{
public:
    /**
     * constructor for different purpose
     */
    RankFourTensor();
    RankFourTensor(const double &val);
    RankFourTensor(const RankFourTensor &a);
    
    /**
     * enum for different initializing
     */
    enum InitMethod{
        InitZero,
        InitIdentity,
        InitIdentity4,
        InitIdentitySymmetric4,
        InitRandom
    };
    RankFourTensor(const InitMethod &method);

    /**
     * get i-k components for jacobian matrix from a rank-4 tensor
     * @param i integer for row-index
     * @param j integer for col-index
     * @param grad_test gradient of the test function
     * @param grad_trial gradient of the tril function
     */
    inline double GetIKjlComponent(const int &i,const int &k,
                                   const Vector3d &grad_test,
                                   const Vector3d &grad_trial)const{
        // function to get K_ik=C_ijkl*N,j*N,l(for assemble local K by using the rank-4 tensor!!!)
        // where j is the index of trial fun
        //       i is the index of test fun
        if(i<1||i>3 || k<1 || k>3){
            MessagePrinter::PrintErrorTxt(" your i or k is out of range for GetIKjlComponent function");
            MessagePrinter::AsFem_Exit();
        }
        return ((*this)(i,1,k,1)*grad_trial(1)
               +(*this)(i,1,k,2)*grad_trial(2)
               +(*this)(i,1,k,3)*grad_trial(3))*grad_test(1)
              +((*this)(i,2,k,1)*grad_trial(1)
               +(*this)(i,2,k,2)*grad_trial(2)
               +(*this)(i,2,k,3)*grad_trial(3))*grad_test(2)
              +((*this)(i,3,k,1)*grad_trial(1)
               +(*this)(i,3,k,2)*grad_trial(2)
               +(*this)(i,3,k,3)*grad_trial(3))*grad_test(3);
    }
    //********************************************
    //*** For operator overload
    //********************************************
    // for index based access
    /**
     * () operator for rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param i i index
     * @param j j index
     * @param k k index
     * @param l l index
     */
    inline double operator()(const int &i,const int &j,const int &k,const int &l) const{
        if(i<1||i>3 || j<1||j>3 || k<1||k>3 || l<1||l>3){
            MessagePrinter::PrintErrorTxt("your i or j or k or l is out of range when you call a rank-4 tensor");
            MessagePrinter::AsFem_Exit();
        }
        return _vals[(((i-1)*_N+j-1)*_N+k-1)*_N+l-1];
    }
    /**
     * () operator for rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param i i index
     * @param j j index
     * @param k k index
     * @param l l index
     */
    inline double& operator()(const int &i,const int &j,const int &k,const int &l){
        if(i<1||i>3 || j<1||j>3 || k<1||k>3 || l<1||l>3){
            MessagePrinter::PrintErrorTxt("your i or j or k or l is out of range when you call a rank-4 tensor");
            MessagePrinter::AsFem_Exit();
        }
        return _vals[(((i-1)*_N+j-1)*_N+k-1)*_N+l-1];
    }
    /**
     * get the voigt component of a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param i i index
     * @param j j index
     */
    double VoigtIJcomponent(const int &i,const int &j) const;
    // for components based access
    /**
     * [] operator to get the component of a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param i i index
     */
    inline double  operator[](const int &i) const{
        return _vals[i-1];
    }
    /**
     * [] operator to get the component of a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param i i index
     */
    inline double& operator[](const int &i){
        return _vals[i-1];
    }
    inline double GetIthVoigtComponent(const int &i)const{
        // TODO: implement efficient voigt access!!!
        return _vals[i];
    }
    //*** for =
    /**
     * = operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side scalar
     */
    inline RankFourTensor& operator=(const double &a){
        for(int i=0;i<_N4;++i) _vals[i]=a;
        return *this;
    }
    /**
     * = operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand rank-4 tensor
     */
    inline RankFourTensor& operator=(const RankFourTensor &a){
        for(int i=0;i<_N4;++i) _vals[i]=a._vals[i];
        return *this;
    }
    //*** for +
    /**
     * + operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side scalar
     */
    inline RankFourTensor operator+(const double &a) const{
        RankFourTensor temp(0.0);
        for(int i=0;i<_N4;++i) temp._vals[i]=_vals[i]+a;
        return temp;
    }
    /**
     * + operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand rank-4 tensor
     */
    inline RankFourTensor operator+(const RankFourTensor &a) const{
        RankFourTensor temp(0.0);
        for(int i=0;i<_N4;++i) temp._vals[i]=_vals[i]+a._vals[i];
        return temp;
    }
    //**** for +=
    /**
     * += operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side scalar
     */
    inline RankFourTensor& operator+=(const double &a){
        for(int i=0;i<_N4;++i) _vals[i]=_vals[i]+a;
        return *this;
    }
    /**
     * += operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side rank-4 tensor
     */
    inline RankFourTensor& operator+=(const RankFourTensor &a){
        for(int i=0;i<_N4;++i) _vals[i]=_vals[i]+a._vals[i];
        return *this;
    }
    //*** for -
    /**
     * - operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side scalar
     */
    inline RankFourTensor operator-(const double &a) const{
        RankFourTensor temp(0.0);
        for(int i=0;i<_N4;++i) temp._vals[i]=_vals[i]-a;
        return temp;
    }
    /**
     * -= operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side rank-4 tensor
     */
    inline RankFourTensor operator-(const RankFourTensor &a) const{
        RankFourTensor temp(0.0);
        for(int i=0;i<_N4;++i) temp._vals[i]=_vals[i]-a._vals[i];
        return temp;
    }
    //**** for -=
    /**
     * -= operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side scalar
     */
    inline RankFourTensor& operator-=(const double &a){
        for(int i=0;i<_N4;++i) _vals[i]=_vals[i]-a;
        return *this;
    }
    /**
     * -= operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side rank-4 tensor
     */
    inline RankFourTensor& operator-=(const RankFourTensor &a){
        for(int i=0;i<_N4;++i) _vals[i]=_vals[i]-a._vals[i];
        return *this;
    }
    //*** for *
    /**
     * * operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side scalar
     */
    inline RankFourTensor operator*(const double &a) const{
        RankFourTensor temp(0.0);
        for(int i=0;i<_N4;++i) temp._vals[i]=_vals[i]*a;
        return temp;
    }
    /**
     * * operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side ran-2 tensor
     */
    RankFourTensor operator*(const RankTwoTensor &a) const;
    /**
     * * operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param lhs left hand side scalar
     * @param a right hand side rank-4 tensor
     */
    friend RankFourTensor operator*(const double &lhs,const RankFourTensor &a);
    /**
     * * operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param lhs left hand side rank-2 tensor
     * @param a right hand side rank-4 tensor
     */
    friend RankFourTensor operator*(const RankTwoTensor &lhs,const RankFourTensor &a);

    //**** for *=
    /**
     * *= operator to a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$
     * @param a right hand side scalar
     */
    inline RankFourTensor& operator*=(const double &a){
        for(int i=0;i<_N4;++i) _vals[i]=_vals[i]*a;
        return *this;
    }
    //*** for double dot operator
    /**
     * double dot : between a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$ and rank-2 tensor
     * @param a right hand side rank-2 tensor
     */
    RankTwoTensor DoubleDot(const RankTwoTensor &a) const;
    /**
     * double dot : between a rank-4 tensor \f$\mathbb{C}_{ijkl}\f$ and rank-4 tensor
     * @param a right hand side rank-2 tensor
     */
    inline RankFourTensor DoubleDot(const RankFourTensor &a) const{
        // C_ijkl=A_ijmn*B_mnkl
        RankFourTensor temp(0.0);
        for(int i=1;i<=_N;++i){
            for(int j=1;j<=_N;++j){
                for(int k=1;k<=_N;++k){
                    for(int l=1;l<=_N;++l){
                        for(int m=1;m<=_N;++m){
                            for(int n=1;n<=_N;++n){
                                temp(i,j,k,l)+=(*this)(i,j,m,n)*a(m,n,k,l);
                            }
                        }
                    }
                }
            }
        }
        return temp;
    }
    //*** for rotated rank-4 tensor
    /**
     * rotat the current a rank-4 tensor by a rank-2 rotation tensor
     * @param rotate right hand side rank-2 tensor
     */
    RankFourTensor Rotate(const RankTwoTensor &rotate) const;
    //*** for push forward and pull back operator
    /**
     * push current rank-4 tensor by a rank-2 rotation tensor F
     * @param F right hand side rank-2 tensor
     */
    RankFourTensor PushByF(const RankTwoTensor &F) const;
    /**
     * pushforward current rank-4 tensor by a rank-2 rotation tensor F
     * @param F right hand side rank-2 tensor
     */
    RankFourTensor PushForward(const RankTwoTensor &F) const;
    //*** for conjugate dot Cijkl=Aijml*Smk
    /**
     * conjugate dot between current rank-4 tensor and a rank-2 tensor
     * @param a right hand side rank-2 tensor
     */
    RankFourTensor ConjDot(const RankTwoTensor &a) const;
    //**********************************************
    //*** some setting functions
    //**********************************************
    /**
     * set current rank-4 tensor to 0
     */
    inline void SetToZeros(){
        for(int i=0;i<_N4;i++) _vals[i]=0.0;
    }
    /**
     * set current rank-4 tensor to idendity
     */
    inline void SetToIdentity(){
        SetToZeros();
        for(int i=1;i<=_N;++i) (*this)(i,i,i,i)=1.0;
    }
    /**
     * set current rank-4 tensor to rank-4 identity
     */
    inline void SetToIdentity4(){
        // maps a rank2 tensor to itself(no symmetric consideriation here),i.e. Iden4:rank2=rank2
        SetToZeros();
        for(int i=1;i<=_N;++i){
            for(int j=1;j<=_N;++j){
                for(int k=1;k<=_N;++k){
                    for(int l=1;l<=_N;++l){
                        (*this)(i,j,k,l)=1.0*((i==k)&&(j==l));
                    }
                }
            }
        }
    }
    /**
     * set current rank-4 tensor to rank-4 transposed identity
     */
    inline void SetIdentity4Transpose(){
        // maps a rank-2 tensor to its transpose, A^{T}=I4:A
        SetToZeros();
        for(int i=1;i<=_N;i++){
            for(int j=1;j<=_N;j++){
                for(int k=1;k<=_N;k++){
                    for(int l=1;l<=_N;l++){
                        (*this)(i,j,k,l)=1.0*((j==k)&&(i==l));
                    }
                }
            }
        }
    }
    /**
     * set current rank-4 tensor to symmetric identity rank-4 tensor
     */
    inline void SetToIdentitySymmetric4(){
        // symmetric fourth-order tensor
        SetToZeros();
        for(int i=1;i<=_N;++i){
            for(int j=1;j<=_N;++j){
                for(int k=1;k<=_N;++k){
                    for(int l=1;l<=_N;++l){
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
    inline void SetToRandom(){
        srand(time(0));
        for(int i=1;i<=_N;++i){
            for(int j=1;j<=_N;++j){
                for(int k=1;k<=_N;++k){
                    for(int l=1;l<=_N;++l){
                        (*this)(i,j,k,l)=static_cast<double>(1.0*rand()/RAND_MAX);
                    }
                }
            }
        }
    }
    //**********************************************
    //*** some fill-in functions
    //**********************************************
    /**
     * eps operator in 2d
     * @param i i index
     * @param j j index
     */
    inline int Eps(const int &i,const int &j) const{
        if(i==1&&j==2) {
            return 1;
        }
        else if(i==2&&j==1){
            return -1;
        }
        return 0;
    }
    /**
     * eps operator in 3d
     * @param i i index
     * @param j j index
     * @param k k index
     */
    inline int Eps(const int &i,const int &j,const int &k) const{
        if(i==1&&j>1&&k>1){
            return Eps(j-1,k-1);
        }
        else if(j==1&&i>1&&k>1){
            return -Eps(i-1,k-1);
        }
        else if(k==1&&i>1&&j>1){
            return Eps(i-1,j-1);
        }
        return 0;
    }
    //*** fill method
    /**
     * set rank-4 tensor from lamme constant and shear modulus
     * @param Lame first Lame constant
     * @param G shear modulus
     */
    void SetFromLameandG(const double &Lame,const double &G);
    /**
     * set rank-4 tensor from Youngs modulus and poisson ratio
     * @param E Youngs modulus
     * @param Nu Poisson ratio
     */
    void SetFromEandNu(const double &E,const double &Nu);
    /**
     * set rank-4 tensor from bulk modulus and shear modulus
     * @param E Bulk modulus
     * @param Nu Shear modulus
     */
    void SetFromKandG(const double &K,const double &G);
    /**
     * set rank-4 tensor from 9-component
     * @param vec 9-component for rank-4 tensor
     */
    void SetFromSymmetric9(const vector<double> &vec);
    /**
     * set rank-4 tensor to Orthotropic
     * @param vec 9-component for rank-4 tensor
     */
    void SetToOrthotropic(const vector<double> &vec);

    //*****************************************
    //*** Print rank-4 tensor
    //*****************************************
    /**
     * print out the rank-4 tensor
     */
    inline void Print() const{
        for(int i=1;i<=_N;++i){
            for(int j=1;j<=_N;++j){
                PetscPrintf(PETSC_COMM_WORLD,"*** i=%d,j=%d: %14.6e  %14.6e  %14.6e***\n",i,j,(*this)(i,j,1,1),(*this)(i,j,1,2),(*this)(i,j,1,3));
                PetscPrintf(PETSC_COMM_WORLD,"***          %14.6e  %14.6e  %14.6e***\n",(*this)(i,j,2,1),(*this)(i,j,2,2),(*this)(i,j,2,3));
                PetscPrintf(PETSC_COMM_WORLD,"***          %14.6e  %14.6e  %14.6e***\n",(*this)(i,j,3,1),(*this)(i,j,3,2),(*this)(i,j,3,3));
            }
        }
    }
    //*******************
    /**
     * print out the rank-4 tensor in voigt notation
     */
    inline void PrintVoigt() const{
        // taken from: https://en.wikipedia.org/wiki/Linear_elasticity
        PetscPrintf(PETSC_COMM_WORLD,"*** %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e***\n",(*this)(1,1,1,1),(*this)(1,1,2,2),(*this)(1,1,3,3),(*this)(1,1,2,3),(*this)(1,1,1,3),(*this)(1,1,1,2));
        PetscPrintf(PETSC_COMM_WORLD,"*** %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e***\n",(*this)(2,2,1,1),(*this)(2,2,2,2),(*this)(2,2,3,3),(*this)(2,2,2,3),(*this)(2,2,1,3),(*this)(2,2,1,2));
        PetscPrintf(PETSC_COMM_WORLD,"*** %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e***\n",(*this)(3,3,1,1),(*this)(3,3,2,2),(*this)(3,3,3,3),(*this)(3,3,2,3),(*this)(3,3,1,3),(*this)(3,3,1,2));
        PetscPrintf(PETSC_COMM_WORLD,"*** %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e***\n",(*this)(2,3,1,1),(*this)(2,3,2,2),(*this)(2,3,3,3),(*this)(2,3,2,3),(*this)(2,3,1,3),(*this)(2,3,1,2));
        PetscPrintf(PETSC_COMM_WORLD,"*** %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e***\n",(*this)(3,1,1,1),(*this)(3,1,2,2),(*this)(3,1,3,3),(*this)(3,1,2,3),(*this)(3,1,1,3),(*this)(3,1,1,2));
        PetscPrintf(PETSC_COMM_WORLD,"*** %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e***\n",(*this)(1,2,1,1),(*this)(1,2,2,2),(*this)(1,2,3,3),(*this)(1,2,2,3),(*this)(1,2,1,3),(*this)(1,2,1,2));
    }
private:
    double _vals[81];
    const int _N=3;
    const int _N2=9;
    const int _N4=81;
};
