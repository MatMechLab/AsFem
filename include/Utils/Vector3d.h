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
//+++ Date   : 2020.10.17
//+++ Purpose: Implement the vector (3 components) operator and 
//+++          calculation
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <iomanip>
#include <cmath>

#include "petsc.h"


//*****************************************
//*** for AsFem's own header files
//*****************************************


using namespace std;



class Vector3d{
public:
    Vector3d(const double &val);
    inline double& operator()(const int &i){
        return _vals[i-1];
    }
    inline double operator()(const int &i)const{
        return _vals[i-1];
    }
    void setZero(){
        _vals[0]=0.0;_vals[1]=0.0;_vals[2]=0.0;
    }
    void setRandom(){
        srand(time(0));
        _vals[0]=static_cast<double>(1.0*random()/RAND_MAX);
        _vals[1]=static_cast<double>(1.0*random()/RAND_MAX);
        _vals[2]=static_cast<double>(1.0*random()/RAND_MAX);
    }
    //***********************************************
    //*** for some basic math operators 
    //***********************************************
    //*** for =
    inline Vector3d& operator=(const double &val){
        _vals[0]=val;_vals[1]=val;_vals[2]=val;
        return *this;
    }
    //*** for +
    inline Vector3d operator+(const double &val){
        Vector3d temp(0.0);
        temp._vals[0]=_vals[0]+val;
        temp._vals[1]=_vals[1]+val;
        temp._vals[2]=_vals[2]+val;
        return temp;
    }
    inline Vector3d operator+(const Vector3d &a){
        Vector3d temp(0.0);
        temp._vals[0]=_vals[0]+a._vals[0];
        temp._vals[1]=_vals[1]+a._vals[1];
        temp._vals[2]=_vals[2]+a._vals[2];
        return temp;
    }
    //*** for +=
    inline Vector3d& operator+=(const double &val){
        _vals[0]=_vals[0]+val;
        _vals[1]=_vals[1]+val;
        _vals[2]=_vals[2]+val;
        return (*this);
    }
    inline Vector3d& operator+=(const Vector3d &a){
        _vals[0]=_vals[0]+a._vals[0];
        _vals[1]=_vals[1]+a._vals[1];
        _vals[2]=_vals[2]+a._vals[2];
        return (*this);
    }
    //***********************************
    //*** for -
    //***********************************
    inline Vector3d operator-(const double &val){
        Vector3d temp(0.0);
        temp._vals[0]=_vals[0]-val;
        temp._vals[1]=_vals[1]-val;
        temp._vals[2]=_vals[2]-val;
        return temp;
    }
    inline Vector3d operator-(const Vector3d &a){
        Vector3d temp(0.0);
        temp._vals[0]=_vals[0]-a._vals[0];
        temp._vals[1]=_vals[1]-a._vals[1];
        temp._vals[2]=_vals[2]-a._vals[2];
        return temp;
    }
    //*** for -=
    inline Vector3d& operator-=(const double &val){
        _vals[0]=_vals[0]-val;
        _vals[1]=_vals[1]-val;
        _vals[2]=_vals[2]-val;
        return (*this);
    }
    inline Vector3d& operator-=(const Vector3d &a){
        _vals[0]=_vals[0]-a._vals[0];
        _vals[1]=_vals[1]-a._vals[1];
        _vals[2]=_vals[2]-a._vals[2];
        return (*this);
    }
    //*** for *
    inline Vector3d operator*(const double &val){
        Vector3d temp(0.0);
        temp._vals[0]=_vals[0]*val;
        temp._vals[1]=_vals[1]*val;
        temp._vals[2]=_vals[2]*val;
        return temp;
    }
    // Please put all the friend funs to the cpp file !!!
    friend Vector3d operator*(const double &val,const Vector3d &a);

    inline double operator*(const Vector3d &a)const{
        return _vals[0]*a._vals[0]
              +_vals[1]*a._vals[1]
              +_vals[2]*a._vals[2];
    }
    inline double ODot(const Vector3d &a)const{
        return _vals[0]*a._vals[0]
            +_vals[1]*a._vals[1]
            +_vals[2]*a._vals[2];
    }
    //*** for *=
    inline Vector3d& operator*=(const double &val){
        _vals[0]=_vals[0]*val;
        _vals[1]=_vals[1]*val;
        _vals[2]=_vals[2]*val;
        return (*this);
    }
    //*** for /
    inline Vector3d operator/(const double &val){
        Vector3d temp(0.0);
        temp._vals[0]=_vals[0]/val;
        temp._vals[1]=_vals[1]/val;
        temp._vals[2]=_vals[2]/val;
        return temp;
    }
    //*** for /=
    inline Vector3d& operator/=(const double &val){
        _vals[0]=_vals[0]/val;
        _vals[1]=_vals[1]/val;
        _vals[2]=_vals[2]/val;
        return (*this);
    }
    //*** for different norm calculation
    inline double normsq()const{
        return _vals[0]*_vals[0]+_vals[1]*_vals[1]+_vals[2]*_vals[2];
    }
    inline double norm()const{
        return sqrt(_vals[0]*_vals[0]+_vals[1]*_vals[1]+_vals[2]*_vals[2]);
    }
    //*** for output operator
    friend ostream& operator<<(ostream& os,const Vector3d &vec);
    inline void Print()const{
        PetscPrintf(PETSC_COMM_WORLD,"%14.6e ,%14.6e ,%14.6e\n",_vals[0],_vals[1],_vals[2]);
    }
private:
    double _vals[3];
};