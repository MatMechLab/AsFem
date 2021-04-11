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
//+++ Date   : 2020.10.17
//+++ Purpose: Implement rank-2 tensor class for some common
//+++          tensor calculation for solid-mechanics problem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <iomanip>
#include <vector>

#include "petsc.h"

//****************************
#include "Utils/MessagePrinter.h"
#include "Utils/Vector3d.h"
#include "Utils/RankFourTensor.h"

using namespace std;


class RankFourTensor;

class RankTwoTensor{
public:
    RankTwoTensor();
    RankTwoTensor(const double &val);
    RankTwoTensor(const RankTwoTensor &a);
    
    enum InitMethod{
        InitZero,
        InitIdentity,
        InitRandom
    };
    RankTwoTensor(const InitMethod &method);
    // this is quite helpful for deformation gradient calculation
    RankTwoTensor(const Vector3d &r1,
                  const Vector3d &r2);// for 2d case
    RankTwoTensor(const Vector3d &r1,
                  const Vector3d &r2,
                  const Vector3d &r3);// for 3d case
    // set from voigt notation
    RankTwoTensor(const double &v11,const double &v22,const double &v12);// for 2d voigt
    RankTwoTensor(const double &v11,const double &v22,const double &v33,
                  const double &v23,const double &v31,const double &v12);// for 3d voigt
    // set from full input
    RankTwoTensor(const double &v11,const double &v12,
                  const double &v21,const double &v22);// for 2d all elements
    RankTwoTensor(const double &v11,const double &v12,const double &v13,
                  const double &v21,const double &v22,const double &v23,
                  const double &v31,const double &v32,const double &v33);// for 3d all elements
    
    //*******************************************************************
    //*** For some basic operators of rank-2 tensor
    //*******************************************************************
    //*** for index based access(start from 1, instead of zero !!!)
    inline double operator()(const int &i,const int &j) const{
        return _vals[(i-1)*_N+j-1];
    }
    inline double& operator()(const int &i,const int &j){
        return _vals[(i-1)*_N+j-1];
    }
    //*** for component based access
    inline double operator[](const int &i) const{
        return _vals[i-1];
    }
    inline double& operator[](const int &i){
        return _vals[i-1];
    }
    //*** for columne and row based operator
    inline Vector3d IthRow(const int &i)const{
        Vector3d temp(0.0);
        temp(1)=(*this)(i,1);
        temp(2)=(*this)(i,2);
        temp(3)=(*this)(i,3);
        return temp;
    }
    inline Vector3d IthCol(const int &i)const{
        Vector3d temp(0.0);
        temp(1)=(*this)(1,i);
        temp(2)=(*this)(2,i);
        temp(3)=(*this)(3,i);
        return temp;
    }
    inline double GetIthVoigtComponent(const int &i)const{
        if(i==1){
            return (*this)(1,1);
        }
        else if(i==2){
            return (*this)(2,2);
        }
        else if(i==3){
            return (*this)(3,3);
        }
        else if(i==4){
            return(*this)(2,3);
        }
        else if(i==5){
            return (*this)(1,3);
        }
        else if(i==6){
            return (*this)(1,2);
        }
        return 0;
    }
    //**************************************************
    //*** for some basic mathematic operators
    //**************************************************
    //*** for =
    inline RankTwoTensor& operator=(const double &a){
        for(int i=0;i<_N2;++i) _vals[i]=a;
        return *this;
    }
    inline RankTwoTensor& operator=(const RankTwoTensor &a){
        for(int i=0;i<_N2;++i) _vals[i]=a._vals[i];
        return *this;
    }
    //*** for + operator
    inline RankTwoTensor operator+(const double &a) const{
        RankTwoTensor temp(0.0);
        for(int i=0;i<_N2;++i) temp._vals[i]=_vals[i]+a;
        return temp;
    }
    inline RankTwoTensor operator+(const RankTwoTensor &a) const{
        RankTwoTensor temp(0.0);
        for(int i=0;i<_N2;++i) temp._vals[i]=_vals[i]+a._vals[i];
        return temp;
    }
    //*** for += operator
    inline RankTwoTensor& operator+=(const double &a) {
        for(int i=0;i<_N2;++i) _vals[i]+=a;
        return *this;
    }
    inline RankTwoTensor& operator+=(const RankTwoTensor &a){
        for(int i=0;i<_N2;++i) _vals[i]+=a._vals[i];
        return *this;
    }
    //*** for - operator
    inline RankTwoTensor operator-(const double &a) const{
        RankTwoTensor temp(0.0);
        for(int i=0;i<_N2;++i) temp._vals[i]=_vals[i]-a;
        return temp;
    }
    inline RankTwoTensor operator-(const RankTwoTensor &a) const{
        RankTwoTensor temp(0.0);
        for(int i=0;i<_N2;++i) temp._vals[i]=_vals[i]-a._vals[i];
        return temp;
    }
    //*** for -= operator
    inline RankTwoTensor& operator-=(const double &a) {
        for(int i=0;i<_N2;++i) _vals[i]-=a;
        return *this;
    }
    inline RankTwoTensor& operator-=(const RankTwoTensor &a){
        for(int i=0;i<_N2;++i) _vals[i]-=a._vals[i];
        return *this;
    }
    //*** for * operator
    //**********************
    inline RankTwoTensor operator*(const double &a) const{
        RankTwoTensor temp(0.0);
        for(int i=0;i<_N2;++i) temp._vals[i]=_vals[i]*a;
        return temp;
    }
    //*** for left hand side scalar times rank-2 tensor
    friend RankTwoTensor operator*(const double &lhs,const RankTwoTensor &a);
    //*** for left hand vector times rank-2 tensor
    friend Vector3d operator*(const Vector3d &lhs,const RankTwoTensor &a);

    inline Vector3d operator*(const Vector3d &a) const{
        Vector3d temp(0.0);
        for(int i=1;i<=_N;++i){
            temp(i)=0.0;
            for(int j=1;j<=_N;++j){
                temp(i)+=(*this)(i,j)*a(j);
            }
        }
        return temp;
    }
    inline RankTwoTensor operator*(const RankTwoTensor &a) const{
        // return A*B(still rank-2 tensor)
        RankTwoTensor temp(0.0);
        for(int i=1;i<=_N;++i){
            for(int j=1;j<=_N;++j){
                temp(i,j)=0.0;
                for(int k=1;k<=_N;++k){
                    temp(i,j)+=(*this)(i,k)*a(k,j);
                }
            }
        }
        return temp;
    }
    //*** for *= operator
    inline RankTwoTensor& operator*=(const double &a) {
        for(int i=0;i<_N2;++i) _vals[i]*=a;
        return *this;
    }
    inline RankTwoTensor& operator*=(const RankTwoTensor &a){
        RankTwoTensor temp=(*this)*a;
        (*this)=temp;
        return *this;
    }
    //*** for cross-dot 
    inline void VectorCrossDot(const double (&a)[3],const double (&b)[3]){
        for(int i=1;i<=_N;++i){
            for(int j=1;j<=_N;++j){
                (*this)(i,j)=a[i-1]*b[j-1];
            }
        }
    }
    inline void VectorCrossDot(const vector<double> &a,const vector<double> &b){
        for(int i=1;i<=_N;++i){
            for(int j=1;j<=_N;++j){
                (*this)(i,j)=a[i-1]*b[j-1];
            }
        }
    }
    inline void VectorCrossDot(const Vector3d &a,const Vector3d &b){
        for(int i=1;i<=_N;++i){
            for(int j=1;j<=_N;++j){
                (*this)(i,j)=a(i)*b(j);
            }
        }
    }
    friend RankTwoTensor VecCrossDot(const Vector3d &a,const Vector3d &b);
    //*** for double dot operator
    inline double DoubleDot(const RankTwoTensor &a) const{
        // return A:B calculation
        double sum=0.0;
        for(int i=1;i<=_N;++i){
            for(int j=1;j<=_N;++j){
                // J.N. Reddy use A_ijB_ji
                sum+=(*this)(i,j)*a(i,j);// use this to get the positive definite scale!!!
            }
        }
        return sum;
    }
    //*** for /
    inline RankTwoTensor operator/(const double &a) const{
        RankTwoTensor temp(0.0);
        for(int i=0;i<_N2;++i) temp._vals[i]=_vals[i]/a;
        return temp;
    }
    //*** for /=
    inline RankTwoTensor& operator/=(const double &a){
        for(int i=0;i<_N2;++i) _vals[i]=_vals[i]/a;
        return *this;
    }
    //********************************************************
    //**** other mathematic related functions
    //********************************************************
    inline double Trace() const{
        return (*this)(1,1)+(*this)(2,2)+(*this)(3,3);
    }
    inline double Det() const{
        // taken from http://mathworld.wolfram.com/Determinant.html
        // Eq.10
        return (*this)(1,1)*(*this)(2,2)*(*this)(3,3)
              -(*this)(1,1)*(*this)(2,3)*(*this)(3,2)
              -(*this)(1,2)*(*this)(2,1)*(*this)(3,3)
              +(*this)(1,2)*(*this)(2,3)*(*this)(3,1)
              +(*this)(1,3)*(*this)(2,1)*(*this)(3,2)
              -(*this)(1,3)*(*this)(2,2)*(*this)(3,1);
    }
    inline double Norm() const{
        double sum=0.0;
        for(int i=0;i<_N2;i++) sum+=_vals[i]*_vals[i];
        return sqrt(sum);
    }
    //*** for the different invariants of stress(strain)
    inline double FirstInvariant() const{
        return Trace();
    }
    inline double SecondInvariant() const{
        double trAA=((*this)*(*this)).Trace();
        double trA=Trace();
        return 0.5*(trA*trA-trAA);
    }
    inline double ThirdInvariant() const{
        return Det();
    }
    //*** for inverse
    inline RankTwoTensor Inverse() const{
        double J=Det();
        if(abs(J)<1.0e-15){
            MessagePrinter::PrintErrorTxt("inverse operation failed for a singular rank-2 tensor !");
            MessagePrinter::AsFem_Exit();
        }
        RankTwoTensor inv(0.0);
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
    inline RankTwoTensor Transpose() const{
        RankTwoTensor temp(0.0);
        for(int i=1;i<=_N;++i){
            for(int j=1;j<=_N;++j){
                temp(i,j)=(*this)(j,i);
            }
        }
        return temp;
    }
    //********************************************************
    //**** For stress and strain decomposition
    //********************************************************
    void CalcEigenValueAndEigenVectors(double (&eigval)[3],RankTwoTensor &eigvec) const;
    RankFourTensor CalcPostiveProjTensor(double (&eigval)[3],RankTwoTensor &eigvec) const;
    RankFourTensor GetPostiveProjTensor() const;
    //*******************************************************************
    //*** some setting functions
    //*******************************************************************
    inline void SetToZeros(){
        for(int i=0;i<_N2;++i) _vals[i]=0.0;
    }
    inline void SetToIdentity(){
        for(int i=1;i<=_N;++i){
            for(int j=1;j<=_N;++j){
                if(i==j){
                    (*this)(i,j)=1.0;
                }
                else{
                    (*this)(i,j)=0.0;
                }
            }
        }
    }
    inline void SetToRandom(){
        srand(time(0));
        for(int i=1;i<=_N;++i){
            for(int j=1;j<=_N;++j){
                (*this)(i,j)=static_cast<double>(1.0*rand()/RAND_MAX);
            }
        }
    }
    // for deformation gradient based calculation(or similar calculation)
    void SetFromGradU(const Vector3d &gradUx);
    void SetFromGradU(const Vector3d &gradUx,const Vector3d &gradUy);
    void SetFromGradU(const Vector3d &gradUx,const Vector3d &gradUy,const Vector3d &gradUz);
    void SetFromVoigt(const double &v11,const double &v22,const double &v12);// for 2D
    void SetFromVoigt(const double &v11,const double &v22,const double &v33,
                      const double &v23,const double &v31,const double &v12);// for 3D
    //********************************************************
    //**** For rotation tensor
    //********************************************************
    void SetRotationTensorFromEulerAngle(const double &theta1,const double &theta2,const double &theta3);
    //*******************************************************************
    //*** some higher order tensor calculation
    //*******************************************************************
    // must be put into cpp file, otherwise the compiler will complain!!!
    RankFourTensor CrossDot(const RankTwoTensor &a) const;
    RankFourTensor ODot(const RankTwoTensor &a) const;

    RankFourTensor IJlkDot(const RankTwoTensor &a) const;

    RankFourTensor IklJDot(const RankTwoTensor &a) const;
    RankFourTensor IkJlDot(const RankTwoTensor &a) const;

    RankFourTensor IlJkDot(const RankTwoTensor &a) const;
    RankFourTensor IlkJDot(const RankTwoTensor &a) const;

    RankFourTensor JkIlDot(const RankTwoTensor &a) const;
    RankFourTensor JklIDot(const RankTwoTensor &a) const;
    
    RankFourTensor klJIDot(const RankTwoTensor &a) const;
    RankFourTensor lkJIDot(const RankTwoTensor &a) const;
    //*******************************************************************
    //*** Print functions
    //*******************************************************************
    inline void Print() const{
        PetscPrintf(PETSC_COMM_WORLD,"*** %14.6e ,%14.6e ,%14.6e ***\n",(*this)(1,1),(*this)(1,2),(*this)(1,3));
        PetscPrintf(PETSC_COMM_WORLD,"*** %14.6e ,%14.6e ,%14.6e ***\n",(*this)(2,1),(*this)(2,2),(*this)(2,3));
        PetscPrintf(PETSC_COMM_WORLD,"*** %14.6e ,%14.6e ,%14.6e ***\n",(*this)(3,1),(*this)(3,2),(*this)(3,3));
    }

private:
    const int _N=3;
    const int _N2=9;
    double _vals[9];
};