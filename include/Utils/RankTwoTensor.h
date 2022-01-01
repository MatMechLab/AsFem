//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
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

#include "Eigen/Eigen"

//****************************
#include "Utils/MessagePrinter.h"
#include "Utils/Vector3d.h"
#include "Utils/RankFourTensor.h"


using namespace std;


class RankFourTensor;

/**
 * The class implement the operator and calculation for rank-2 tensor
 */
class RankTwoTensor{
public:
    /**
     * Constructor for different purpose
     */
    RankTwoTensor();
    RankTwoTensor(const double &val);
    RankTwoTensor(const RankTwoTensor &a);
   
    /**
     * different intialization method for rank-2 tensor
     */ 
    enum InitMethod{
        InitZero,
        InitIdentity,
        InitRandom
    };
    /**
     * Constructor with given initializing method
     */
    RankTwoTensor(const InitMethod &method);
    // this is quite helpful for deformation gradient calculation
    /**
     * Constructor with given vector for 2d case
     */
    RankTwoTensor(const Vector3d &r1,
                  const Vector3d &r2);// for 2d case
    /**
     * Constructor with given vector for 3d case
     */
    RankTwoTensor(const Vector3d &r1,
                  const Vector3d &r2,
                  const Vector3d &r3);// for 3d case
    /**
     * Constructor with given voigt components in 2d case
     */
    RankTwoTensor(const double &v11,const double &v22,const double &v12);// for 2d voigt
    /**
     * Constructor with given voigt components in 3d case
     */
    RankTwoTensor(const double &v11,const double &v22,const double &v33,
                  const double &v23,const double &v31,const double &v12);// for 3d voigt
    // set from full input
    /**
     * Constructor with the full components (4) for 2d case
     */
    RankTwoTensor(const double &v11,const double &v12,
                  const double &v21,const double &v22);// for 2d all elements
    /**
     * Constructor with the full components (9) for 3d case
     */
    RankTwoTensor(const double &v11,const double &v12,const double &v13,
                  const double &v21,const double &v22,const double &v23,
                  const double &v31,const double &v32,const double &v33);// for 3d all elements
    
    //*******************************************************************
    //*** For some basic operators of rank-2 tensor
    //*******************************************************************
    /** for index based access(start from 1, instead of zero !!!)
     * @param i i index of the rank-2 tensor, start from 1
     * @param j j index of the rank-2 tensor, start from 1
     */
    inline double operator()(const int &i,const int &j) const{
        if(i<1||i>3 || j<1||j>3 ){
            MessagePrinter::PrintErrorTxt("your i or j is out of range when you call a rank-2 tensor");
            MessagePrinter::AsFem_Exit();
        }
        return _vals[(i-1)*_N+j-1];
    }
    /** for index based access(start from 1, instead of zero !!!)
     * @param i i index of the rank-2 tensor, start from 1
     * @param j j index of the rank-2 tensor, start from 1
     */
    inline double& operator()(const int &i,const int &j){
        if(i<1||i>3 || j<1||j>3 ){
            MessagePrinter::PrintErrorTxt("your i or j is out of range when you call a rank-2 tensor");
            MessagePrinter::AsFem_Exit();
        }
        return _vals[(i-1)*_N+j-1];
    }
    //*** for component based access
    /**
     * [] operator for element access
     * @param i global index, range from 1~9
     */
    inline double operator[](const int &i) const{
        return _vals[i-1];
    }
    /**
     * [] operator for element access
     * @param i global index, range from 1~9
     */
    inline double& operator[](const int &i){
        return _vals[i-1];
    }
    //*** for columne and row based operator
    /**
     * get the ith row of the rank-2 tensor
     * @param i i-th row number, start from 1 to 3
     */
    inline Vector3d IthRow(const int &i)const{
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
    inline Vector3d IthCol(const int &i)const{
        Vector3d temp(0.0);
        temp(1)=(*this)(1,i);
        temp(2)=(*this)(2,i);
        temp(3)=(*this)(3,i);
        return temp;
    }
    /**
     * get i-th voigt components
     * @param i the index of the voigt notation
     */
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
    /**
     * '=' operator for scalar
     * @param a right hand side scalar
     */
    inline RankTwoTensor& operator=(const double &a){
        for(int i=0;i<_N2;++i) _vals[i]=a;
        return *this;
    }
    /**
     * '=' operator for rank-2 tensor
     * @param a the right hand side rank-2 tensor
     */
    inline RankTwoTensor& operator=(const RankTwoTensor &a){
        for(int i=0;i<_N2;++i) _vals[i]=a._vals[i];
        return *this;
    }
    //*** for + operator
    /**
     * '+' operator for scalar
     * @param a right hand side scalar
     */
    inline RankTwoTensor operator+(const double &a) const{
        RankTwoTensor temp(0.0);
        for(int i=0;i<_N2;++i) temp._vals[i]=_vals[i]+a;
        return temp;
    }
    /**
     * '+' operator for rank-2 tensor
     * @param a right hand side rank-2 tensor
     */
    inline RankTwoTensor operator+(const RankTwoTensor &a) const{
        RankTwoTensor temp(0.0);
        for(int i=0;i<_N2;++i) temp._vals[i]=_vals[i]+a._vals[i];
        return temp;
    }
    //*** for += operator
    /**
     * '+=' operator for scalar
     * @param a right hand side scalar
     */
    inline RankTwoTensor& operator+=(const double &a) {
        for(int i=0;i<_N2;++i) _vals[i]+=a;
        return *this;
    }
    /**
     * '+=' for rank-2 tensor
     * @param a right hand side tensor
     */
    inline RankTwoTensor& operator+=(const RankTwoTensor &a){
        for(int i=0;i<_N2;++i) _vals[i]+=a._vals[i];
        return *this;
    }
    //*** for - operator
    /**
     * '-' for scalar
     * @param a right hand side scalar
     */
    inline RankTwoTensor operator-(const double &a) const{
        RankTwoTensor temp(0.0);
        for(int i=0;i<_N2;++i) temp._vals[i]=_vals[i]-a;
        return temp;
    }
    /**
     * '-' for rank-2 tensor
     * @param a right hand side rank-2 tensor
     */
    inline RankTwoTensor operator-(const RankTwoTensor &a) const{
        RankTwoTensor temp(0.0);
        for(int i=0;i<_N2;++i) temp._vals[i]=_vals[i]-a._vals[i];
        return temp;
    }
    //*** for -= operator
    /**
     * '-=' operator for scalar
     * @param a right hand side scalar
     */
    inline RankTwoTensor& operator-=(const double &a) {
        for(int i=0;i<_N2;++i) _vals[i]-=a;
        return *this;
    }
    /**
     * '-=' operator for rank-2 tensor
     * @param a right hand side rank-2 tensor
     */
    inline RankTwoTensor& operator-=(const RankTwoTensor &a){
        for(int i=0;i<_N2;++i) _vals[i]-=a._vals[i];
        return *this;
    }
    //*** for * operator
    //**********************
    /**
     * '*' for scalar
     * @param a right hand side scalar
     */
    inline RankTwoTensor operator*(const double &a) const{
        RankTwoTensor temp(0.0);
        for(int i=0;i<_N2;++i) temp._vals[i]=_vals[i]*a;
        return temp;
    }
    //*** for left hand side scalar times rank-2 tensor
    /**
     * '*' for left hand side scalar
     * @param lhs left hand side scalar value
     * @param a right hand side rank-2 tensor
     */
    friend RankTwoTensor operator*(const double &lhs,const RankTwoTensor &a);
    //*** for left hand vector times rank-2 tensor
    /**
     * '*' operator for left hand side vector3d
     * @param lhs the left hand side vector3d
     * @param a the right hand side rank-2 tensor
     */
    friend Vector3d operator*(const Vector3d &lhs,const RankTwoTensor &a);

    /**
     * '*' operator for vector3d
     * @param a right hand side vector3d
     */
    inline Vector3d operator*(const Vector3d &a) const{
        Vector3d temp(0.0);
        for(int i=1;i<=_N;++i){
            temp(i)=(*this)(i,1)*a(1)+(*this)(i,2)*a(2)+(*this)(i,3)*a(3);
        }
        return temp;
    }
    /**
     * '*' operator for rank-2 tensor, return \f$a_{ij}=b_{ik}c_{kj}\f$
     * @param a the right hand side rank-2 tensor
     */
    inline RankTwoTensor operator*(const RankTwoTensor &a) const{
        // return A*B(still rank-2 tensor)
        RankTwoTensor temp(0.0);
        for(int i=1;i<=_N;++i){
            for(int j=1;j<=_N;++j){
                temp(i,j)=(*this)(i,1)*a(1,j)+(*this)(i,2)*a(2,j)+(*this)(i,3)*a(3,j);
            }
        }
        return temp;
    }
    //*** for *= operator
    /**
     * '*=' operator for scalar
     * @param a right hand side scalar
     */
    inline RankTwoTensor& operator*=(const double &a) {
        for(int i=0;i<_N2;++i) _vals[i]*=a;
        return *this;
    }
    /**
     * '*=' for rank-2 tensor
     * @param a right hand side rank-2 tensor
     */
    inline RankTwoTensor& operator*=(const RankTwoTensor &a){
        RankTwoTensor temp=(*this)*a;
        (*this)=temp;
        return *this;
    }
    //*** for cross-dot
    /**
     * cross dot \f$\otimes\f$ for two double array, set current one to \f$c_{ij}=a_{i}b_{j}\f$.
     * @param a double array for 1st dimension
     * @param b double array for 2nd dimension
     */
    inline void VectorOTimes(const double (&a)[3],const double (&b)[3]){
        for(int i=1;i<=_N;++i){
            for(int j=1;j<=_N;++j){
                (*this)(i,j)=a[i-1]*b[j-1];
            }
        }
    }
    
    /**
     * cross dot \f$\otimes\f$ for two vector(from STL), set current one to \f$c_{ij}=a_{i}b_{j}\f$.
     * @param a vector<double> for 1st dimension
     * @param b vector<double> for 2nd dimension
     */
    inline void VectorOTimes(const vector<double> &a,const vector<double> &b){
        for(int i=1;i<=_N;++i){
            for(int j=1;j<=_N;++j){
                (*this)(i,j)=a[i-1]*b[j-1];
            }
        }
    }
    
    /**
     * cross dot \f$\otimes\f$ for two vector3d, set current one to \f$c_{ij}=a_{i}b_{j}\f$.
     * @param a vector3d for 1st dimension
     * @param b vector3d for 2nd dimension
     */
    inline void VectorOTimes(const Vector3d &a,const Vector3d &b){
        for(int i=1;i<=_N;++i){
            for(int j=1;j<=_N;++j){
                (*this)(i,j)=a(i)*b(j);
            }
        }
    }

    /**
     * cross dot \f$\otimes\f$ for two vector3d, return \f$c_{ij}=a_{i}b_{j}\f$.
     * @param a vector3d for 1st dimension
     * @param b vector3d for 2nd dimension
     * return a new rank-2 tensor
     */
    friend RankTwoTensor VecOTimes(const Vector3d &a,const Vector3d &b);
    //*** for double dot operator
    /**
     * double dot(:, the double contruction) between two rank-2 tensor, the result is \f$\sum a_{ij}b_{ij}\f$, which is a scalar
     * @param a the right hand side rank-2 tensor
     */
    inline double DoubleDot(const RankTwoTensor &a) const{
        // return A:B calculation
        double sum=0.0;
        for(int i=1;i<=_N;++i){
            for(int j=1;j<=_N;++j){
                // You may see A:B=A_ijB_ji in other books/literature, here we use A_ijB_ij
                // to keep the same, in Rank4Tensor, we follow the same definition!
                sum+=(*this)(i,j)*a(i,j);// use this to get the positive definite case!!!
            }
        }
        return sum;
    }
    //*** for /
    /**
     * '/' operator for scalar
     * @param a the right hand side scalar
     */
    inline RankTwoTensor operator/(const double &a) const{
        RankTwoTensor temp(0.0);
        if(abs(a)<1.0e-16){
            MessagePrinter::PrintErrorTxt("rank-2/0 is singular, the input scalar is zero, which is invalid for '/' operator of a rank-2 tensor");
            MessagePrinter::AsFem_Exit();
        }
        for(int i=0;i<_N2;++i) temp._vals[i]=_vals[i]/a;
        return temp;
    }
    //*** for /=
    /**
     * '/=' operator for scalar
     * @param right hand side scalar
     */
    inline RankTwoTensor& operator/=(const double &a){
        if(abs(a)<1.0e-16){
            MessagePrinter::PrintErrorTxt("rank-2/0 is singular, the input scalar is zero, which is invalid for '/=' operator of a rank-2 tensor");
            MessagePrinter::AsFem_Exit();
        }
        for(int i=0;i<_N2;++i) _vals[i]=_vals[i]/a;
        return *this;
    }
    //********************************************************
    //**** other mathematic related functions
    //********************************************************
    /**
     * return the trace of a rank-2 tensor, result is \f$\sum a_{ii}\f$
     */
    inline double Trace() const{
        return (*this)(1,1)+(*this)(2,2)+(*this)(3,3);
    }
    /**
     * return the determinant of the current rank-2 tensor
     */
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
    /**
     * return the \f$L_{2}\f$ norm of current rank-2 tensor, result is \f$\sqrt{\sum a_{ij}^{2}}\f$
     */
    inline double Norm() const{
        double sum=0.0;
        for(int i=0;i<_N2;i++) sum+=_vals[i]*_vals[i];
        return sqrt(sum);
    }
    
    /**
     * return the \f$L_{2}\f$ norm^2 of current rank-2 tensor, result is \f$\sqrt{\sum a_{ij}^{2}}\f$
     */
    inline double Norm2() const{
        double sum=0.0;
        for(int i=0;i<_N2;i++) sum+=_vals[i]*_vals[i];
        return sum;
    }

    //*** for the different invariants of stress(strain)
    /**
     * return the first invariant of current rank-2 tensor, namely, \f$I_{1}\f$.
     */
    inline double FirstInvariant() const{
        return Trace();
    }
    /**
     * return the second invariant of current tensor, namely, \f$I_{2}\f$.
     */
    inline double SecondInvariant() const{
        double trAA=((*this)*(*this)).Trace();
        double trA=Trace();
        return 0.5*(trA*trA-trAA);
    }
    /**
     * return the third invariant of current tensor, namely, \f$I_{3}\f$.
     */
    inline double ThirdInvariant() const{
        return Det();
    }
    //*** for inverse
    /**
     * return the inverse tensor of current one, namely, \f$\mathbf{A}^{-1}\f$.
     */
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
    /**
     * return the transpose tensor of current one, namely, \f$\mathbf{A}^{T}\f$.
     */
    inline RankTwoTensor Transpose() const{
        RankTwoTensor temp(0.0);
        for(int i=1;i<=_N;++i){
            for(int j=1;j<=_N;++j){
                temp(i,j)=(*this)(j,i);
            }
        }
        return temp;
    }
    /**
     * transpose current tensor, and overwrite the original one
     */
    inline void Transposed(){
        RankTwoTensor temp(0.0);
        temp=(*this).Transpose();
        (*this)=temp;
    }
    /**
     * return the deviatoric part of current rank-2 tensor
     */
    inline RankTwoTensor Dev()const{
        RankTwoTensor temp(0.0),I(0.0);
        I.SetToIdentity();
        temp=(*this)-I*(this->Trace()/3.0);
        return temp;
    }
    //********************************************************
    //**** For stress and strain decomposition
    //********************************************************
    /**
     * calculate the eigen value and eigen vector for current rank-2 tensor
     * @param eigval the double array, which stores the eigen value
     * @param eigvec the rank-2 tensor, where each column store the related eigen vector
     */
    void CalcEigenValueAndEigenVectors(double (&eigval)[3],RankTwoTensor &eigvec) const;
    /**
     * calculate the positive projection tensor(a rank-4 tensor), this algorithm is taken from Miehe's paper, for the details, please see the cpp file
     * @param eigval the double array, which stores the eigen value
     * @param eigvec the rank-2 tensor, whoses column stores the related eigen vector
     */
    RankFourTensor CalcPositiveProjTensor(double (&eigval)[3],RankTwoTensor &eigvec) const;
    /**
     * calculate the positive projection tensor (rank-4 tensor) based on current rank-2 tensor
     */
    RankFourTensor GetPositiveProjTensor() const;
    //*******************************************************************
    //*** some setting functions
    //*******************************************************************
    /**
     * set the elements of current rank-2 tensor to be zero
     */
    inline void SetToZeros(){
        for(int i=0;i<_N2;++i) _vals[i]=0.0;
    }
    /**
     * set current rannk-2 tensor to be an identitiy tensor, where \f$a_{ij}=\delta_{ij}\f$.
     */
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
    /**
     * set current rank-2 tensor to be a random one
     */
    inline void SetToRandom(){
        srand(time(0));
        for(int i=1;i<=_N;++i){
            for(int j=1;j<=_N;++j){
                (*this)(i,j)=static_cast<double>(1.0*rand()/RAND_MAX);
            }
        }
    }
    // for deformation gradient based calculation(or similar calculation)
    /**
     * fill up current rank-2 tensor from the displacement gradient
     * @param gradUx the gradient of Ux, namely, \f$\nabla u_{x}\f$
     */ 
    void SetFromGradU(const Vector3d &gradUx);
    /**
     * fill up current rank-2 tensor from the gradient of ux and uy
     * @param gradUx the gradient of Ux, namely, \f$\nabla u_{x}\f$
     * @param gradUy the gradient of Uy, namely, \f$\nabla u_{y}\f$
     */
    void SetFromGradU(const Vector3d &gradUx,const Vector3d &gradUy);
    /**
     * fill up current rank-2 tensor from the gradient of ux and uy
     * @param gradUx the gradient of Ux, namely, \f$\nabla u_{x}\f$
     * @param gradUy the gradient of Uy, namely, \f$\nabla u_{y}\f$
     */
    void SetFromGradU(const Vector3d &gradUx,const Vector3d &gradUy,const Vector3d &gradUz);
    /**
     * fill up current rank-2 tensor from voigt notation elements for 2d case 
     * @param v11 the voigt notation marked element
     */
    void SetFromVoigt(const double &v11,const double &v22,const double &v12);// for 2D
    /**
     * fill up current rank-2 tensor from voigt notation elements for 3d  case
     */ 
    void SetFromVoigt(const double &v11,const double &v22,const double &v33,
                      const double &v23,const double &v31,const double &v12);// for 3D
    //********************************************************
    //**** For rotation tensor
    //********************************************************
    /**
     * fill up current rank-2 tensor with eular angle, which should be used for the rotation tensor
     * @param theta1 the first eular angle, the angle must be degree
     * @param theta2 the second eular angle, the angle must be degree
     * @param theta3 the third eular angle, the angle must be degree
     */
    void SetRotationTensorFromEulerAngle(const double &theta1,const double &theta2,const double &theta3);
    //*******************************************************************
    //*** some higher order tensor calculation
    //*******************************************************************
    // must be put into cpp file, otherwise the compiler will complain!!!
    /**
     * Otime(\f$\otimes\f$) between two rank-2 tensor, which will return \f$c_{ijkl}=a_{ij}b_{kl}\f$
     */
    RankFourTensor OTimes(const RankTwoTensor &a) const;
    /**
     * Odot (\f$\odot\f$) between two rank-2 tensor, which will return \f$c_{ijkl}=\frac{1}{2}(a_{ik}b_{jl}+a_{il}b_{jk})\f$
     */
    RankFourTensor ODot(const RankTwoTensor &a) const;

    /**
     * IJlk times (\f$IJ\otimes lk\f$) operator
     */
    RankFourTensor IJlkTimes(const RankTwoTensor &a) const;

    /**
     * IklJ times (\f$Ik\otimes lJ\f$)
     */
    RankFourTensor IklJTimes(const RankTwoTensor &a) const;
    /**
     * IkJl times (\f$Ik\otimes Jl\f$)
     */
    RankFourTensor IkJlTimes(const RankTwoTensor &a) const;

    /**
     * IlJk times (\f$Il\otimes Jk\f$)
     */
    RankFourTensor IlJkTimes(const RankTwoTensor &a) const;
    /**
     * IlkJ times (\f$Il\otimes kJ\f$)
     */
    RankFourTensor IlkJTimes(const RankTwoTensor &a) const;
    
    /**
     * JkIl times (\f$Jk\otimes Il\f$)
     */
    RankFourTensor JkIlTimes(const RankTwoTensor &a) const;
    /**
     * JklI times (\f$Jk\otimes lI\f$)
     */
    RankFourTensor JklITimes(const RankTwoTensor &a) const;
    
    /**
     * klJI times (\f$kl\otimes JI\f$)
     */
    RankFourTensor klJITimes(const RankTwoTensor &a) const;
    /**
     * lkJI times (\f$lk\otimes JI\f$)
     */
    RankFourTensor lkJITimes(const RankTwoTensor &a) const;
    //*******************************************************************
    //*** Print functions
    //*******************************************************************
    /**
     * print all the elements to the terminal
     */
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
