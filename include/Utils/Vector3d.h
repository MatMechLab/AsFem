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
//+++ Purpose: Implement the vector (3 components) operator and 
//+++          calculation
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <cmath>

#include "petsc.h"


//*****************************************
//*** for AsFem's own header files
//*****************************************


using namespace std;

/**
 * This class defines the vector with only 3 componenets
 */
class Vector3d{
public:
    /**
     * Constructor for vector3d class
     */
    Vector3d();
    Vector3d(const double &val);
    Vector3d(const Vector3d &a);
    
    /**
     * () operator for the reference vector3d
     * @param i the i-th index, which must start from 1(\f$ 1\leq i\leq\f$)
     */
    inline double& operator()(const int &i){
        return _vals[i-1];
    }
    /**
     * () operator for the copy of vector3d
     * @param i the i-th index, i>=1 and i<=3
     */
    inline double operator()(const int &i)const{
        return _vals[i-1];
    }
    
    /**
     * this function will set the vector3d to be zero
     */
    void setZero(){
        _vals[0]=0.0;_vals[1]=0.0;_vals[2]=0.0;
    }
    /**
     * this function will fill the vector with 3 different random number
     */
    void setRandom(){
        srand(time(0));
        _vals[0]=static_cast<double>(1.0*rand()/RAND_MAX);
        _vals[1]=static_cast<double>(1.0*rand()/RAND_MAX);
        _vals[2]=static_cast<double>(1.0*rand()/RAND_MAX);
    }
    //***********************************************
    //*** for some basic math operators 
    //***********************************************
    //*** for =
    /**
     * '=' operator for vector3d with a scalar value
     * @param val the scalar value on the right side
     */
    inline Vector3d& operator=(const double &val){
        _vals[0]=val;_vals[1]=val;_vals[2]=val;
        return *this;
    }
    /**
     * '=' operator for vector3d with another vector3d
     * @param a the right-hand side vector3d
     */
    inline Vector3d& operator=(const Vector3d &a){
        _vals[0]=a._vals[0];_vals[1]=a._vals[1];_vals[2]=a._vals[2];
        return *this;
    }
    //*** for +
    /**
     * '+' operator for vector3d with a scalar
     * @param val scalar value on the right side
     */
    inline Vector3d operator+(const double &val){
        Vector3d temp(0.0);
        temp._vals[0]=_vals[0]+val;
        temp._vals[1]=_vals[1]+val;
        temp._vals[2]=_vals[2]+val;
        return temp;
    }
    /**
     * '+' operator for vector3d with another vector3d
     * @param a the right-hand side vector
     */
    inline Vector3d operator+(const Vector3d &a){
        Vector3d temp(0.0);
        temp._vals[0]=_vals[0]+a._vals[0];
        temp._vals[1]=_vals[1]+a._vals[1];
        temp._vals[2]=_vals[2]+a._vals[2];
        return temp;
    }
    //*** for +=
    /**
     * for '+=' operator of vector3d with a scalar
     * @param val the scalar value on the right hand side
     */
    inline Vector3d& operator+=(const double &val){
        _vals[0]=_vals[0]+val;
        _vals[1]=_vals[1]+val;
        _vals[2]=_vals[2]+val;
        return (*this);
    }
    /**
     * for '+=' operator of vector3d with another vector3d
     * @param a the right-hand side vector3d
     */
    inline Vector3d& operator+=(const Vector3d &a){
        _vals[0]=_vals[0]+a._vals[0];
        _vals[1]=_vals[1]+a._vals[1];
        _vals[2]=_vals[2]+a._vals[2];
        return (*this);
    }
    //***********************************
    //*** for -
    //***********************************
    /**
     * '-' operator for vector3d with a scalar
     * @param val the right-hand side scalar value
     */
    inline Vector3d operator-(const double &val){
        Vector3d temp(0.0);
        temp._vals[0]=_vals[0]-val;
        temp._vals[1]=_vals[1]-val;
        temp._vals[2]=_vals[2]-val;
        return temp;
    }
    /**
     * '-' operator for vector3d with another vector3d
     * @param a the right hand side vector3d
     */
    inline Vector3d operator-(const Vector3d &a){
        Vector3d temp(0.0);
        temp._vals[0]=_vals[0]-a._vals[0];
        temp._vals[1]=_vals[1]-a._vals[1];
        temp._vals[2]=_vals[2]-a._vals[2];
        return temp;
    }
    /**
     * for '-=' operator of vector3d with a scalar
     * @param val right hand side scalar
     */
    inline Vector3d& operator-=(const double &val){
        _vals[0]=_vals[0]-val;
        _vals[1]=_vals[1]-val;
        _vals[2]=_vals[2]-val;
        return (*this);
    }
    /**
     * '-=' operator of vector3d with another vector3d
     * @param a right hand side vector3d
     */
    inline Vector3d& operator-=(const Vector3d &a){
        _vals[0]=_vals[0]-a._vals[0];
        _vals[1]=_vals[1]-a._vals[1];
        _vals[2]=_vals[2]-a._vals[2];
        return (*this);
    }
    //*** for *
    /**
     * for '*' operator of vector3d with a scalar
     * @param val right hand side scalar
     */
    inline Vector3d operator*(const double &val){
        Vector3d temp(0.0);
        temp._vals[0]=_vals[0]*val;
        temp._vals[1]=_vals[1]*val;
        temp._vals[2]=_vals[2]*val;
        return temp;
    }
    // Please put all the friend funs to the cpp file !!!
    /**
     * the left hand side '*' of vector3d with a scalar
     * @param val the left hand side scalar
     * @param a the right hand side vector3d
     * this will return a new vector3d
     */
    friend Vector3d operator*(const double &val,const Vector3d &a);

    /**
     * '*' operator of vector3d with another vector3d
     * @param right hand side vector3d
     * this will return a scalar value \f$\sum_{i=1}{3}b_{i}*a_{i}\f$
     */
    inline double operator*(const Vector3d &a)const{
        return _vals[0]*a._vals[0]
              +_vals[1]*a._vals[1]
              +_vals[2]*a._vals[2];
    }
    /**
     * Save function as '*' operator between two vector3d
     * @param a the input vector3d
     */
    inline double ODot(const Vector3d &a)const{
        return _vals[0]*a._vals[0]
            +_vals[1]*a._vals[1]
            +_vals[2]*a._vals[2];
    }
    //*** for *=
    /**
     * '*=' operator of vector3d with a scalar
     * @param val the right hand side scalar
     */
    inline Vector3d& operator*=(const double &val){
        _vals[0]=_vals[0]*val;
        _vals[1]=_vals[1]*val;
        _vals[2]=_vals[2]*val;
        return (*this);
    }
    //*** for /
    /**
     * '/' operator of vector3d with a scalar
     * @param right hand side scalar
     */
    inline Vector3d operator/(const double &val){
        Vector3d temp(0.0);
        temp._vals[0]=_vals[0]/val;
        temp._vals[1]=_vals[1]/val;
        temp._vals[2]=_vals[2]/val;
        return temp;
    }
    //*** for /=
    /**
     * '/=' operator of vector3d with a scalar
     * @param val the right hand side scalar
     */
    inline Vector3d& operator/=(const double &val){
        _vals[0]=_vals[0]/val;
        _vals[1]=_vals[1]/val;
        _vals[2]=_vals[2]/val;
        return (*this);
    }
    //*** for different norm calculation
    /**
     * get the square of \f$L_{2}\f$ norm of vector3d, the result is \f$\sum_{i=1}{3}a_{i}*a_{i}\f$
     */
    inline double normsq()const{
        return _vals[0]*_vals[0]+_vals[1]*_vals[1]+_vals[2]*_vals[2];
    }
    /**
     * the \f$L_{2}\f$ norm of vector3d
     */
    inline double norm()const{
        return sqrt(_vals[0]*_vals[0]+_vals[1]*_vals[1]+_vals[2]*_vals[2]);
    }
    //*** for output operator
    /**
     * '<<' operator override for cout output
     */
    friend ostream& operator<<(ostream& os,const Vector3d &vec);
    
    /**
     * print the components of vector3d
     */
    inline void Print()const{
        PetscPrintf(PETSC_COMM_WORLD,"%14.6e ,%14.6e ,%14.6e\n",_vals[0],_vals[1],_vals[2]);
    }
private:
    double _vals[3];/**< the componenets of vector3d*/
};
