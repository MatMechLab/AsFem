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
//+++ Date   : 2022.05.13
//+++ Purpose: defines the vector with only 3-components, this will
//+++          be frequently used in shape function calculation.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once


#include <cmath>

#include "Utils/MessagePrinter.h"

using std::sqrt;
using std::abs;
using std::fill;


/**
 * This class defines the vector with only 3-components
 */
class Vector3d{
public:
    /**
     * constructor
     */
    Vector3d();
    /**
     * @param val the scalar value
    */
    Vector3d(const double &val);
    /**
     * @param a the right hand vector
    */
    Vector3d(const Vector3d &a);
    //****************************************************
    //*** for operators
    //****************************************************
    /**
     * () operator
     * @param i index
     */
    inline double& operator()(const int &i){
        if(i<1||i>3){
            MessagePrinter::printErrorTxt(to_string(i)+" is out of range for Vector3d");
            MessagePrinter::exitAsFem();
        }
        return m_Vals[i-1];
    }
    /**
     * const () operator
     * @param i index
     */
    inline double operator()(const int &i)const{
        if(i<1||i>3){
            MessagePrinter::printErrorTxt(to_string(i)+" is out of range for Vector3d");
            MessagePrinter::exitAsFem();
        }
        return m_Vals[i-1];
    }
    /**
     * = operator
     * @param val right hand side double value
     */
    inline Vector3d& operator=(const double &val){
        m_Vals[0]=val;m_Vals[1]=val;m_Vals[2]=val;
        return *this;
    }
    /**
     * = operator
     * @param val right hand side Vector3 value
     */
    inline Vector3d& operator=(const Vector3d &a){
        m_Vals[0]=a.m_Vals[0];m_Vals[1]=a.m_Vals[1];m_Vals[2]=a.m_Vals[2];
        return *this;
    }
    //***********************************************
    /**
     * + operator
     * @param val right hand side double value
     */
    inline Vector3d operator+(const double &val)const{
        Vector3d temp(0.0);
        temp.m_Vals[0]=m_Vals[0]+val;temp.m_Vals[1]=m_Vals[1]+val;temp.m_Vals[2]=m_Vals[2]+val;
        return temp;
    }
    /**
     * + operator
     * @param val right hand side double value
     */
    inline Vector3d operator+(const Vector3d &a)const{
        Vector3d temp(0.0);
        temp.m_Vals[0]=m_Vals[0]+a.m_Vals[0];temp.m_Vals[1]=m_Vals[1]+a.m_Vals[1];temp.m_Vals[2]=m_Vals[2]+a.m_Vals[2];
        return temp;
    }
    //*************************************************
    /**
     * += operator
     * @param val right hand side double value
     */
    inline Vector3d& operator+=(const double &val){
        m_Vals[0]+=val;m_Vals[1]+=val;m_Vals[2]+=val;
        return *this;
    }
    /**
     * += operator
     * @param val right hand side Vector3 value
     */
    inline Vector3d& operator+=(const Vector3d &a){
        m_Vals[0]+=a.m_Vals[0];m_Vals[1]+=a.m_Vals[1];m_Vals[2]+=a.m_Vals[2];
        return *this;
    }
    //***********************************************
    /**
     * - operator
     * @param val right hand side double value
     */
    inline Vector3d operator-(const double &val)const{
        Vector3d temp(0.0);
        temp.m_Vals[0]=m_Vals[0]-val;temp.m_Vals[1]=m_Vals[1]-val;temp.m_Vals[2]=m_Vals[2]-val;
        return temp;
    }
    /**
     * - operator
     * @param val right hand side double value
     */
    inline Vector3d operator-(const Vector3d &a)const{
        Vector3d temp(0.0);
        temp.m_Vals[0]=m_Vals[0]-a.m_Vals[0];temp.m_Vals[1]=m_Vals[1]-a.m_Vals[1];temp.m_Vals[2]=m_Vals[2]-a.m_Vals[2];
        return temp;
    }
    //*************************************************
    /**
     * -= operator
     * @param val right hand side double value
     */
    inline Vector3d& operator-=(const double &val){
        m_Vals[0]-=val;m_Vals[1]-=val;m_Vals[2]-=val;
        return *this;
    }
    /**
     * -= operator
     * @param val right hand side Vector3 value
     */
    inline Vector3d& operator-=(const Vector3d &a){
        m_Vals[0]-=a.m_Vals[0];m_Vals[1]-=a.m_Vals[1];m_Vals[2]-=a.m_Vals[2];
        return *this;
    }
    //***********************************************
    /**
     * * operator
     * @param val right hand side double value
     */
    inline Vector3d operator*(const double &val)const{
        Vector3d temp(0.0);
        temp.m_Vals[0]=m_Vals[0]*val;temp.m_Vals[1]=m_Vals[1]*val;temp.m_Vals[2]=m_Vals[2]*val;
        return temp;
    }
    /**
     * * operator
     * @param val right hand side Vector3 value
     */
    inline double operator*(const Vector3d &a)const{
        double sum=static_cast<double>(m_Vals[0]*a.m_Vals[0]
                                      +m_Vals[1]*a.m_Vals[1]
                                      +m_Vals[2]*a.m_Vals[2]);
        return sum;
    }
    /**
     * left hand side * operator with scalar
    */
    friend Vector3d operator*(const double &val,const Vector3d &a);
    /**
     * same function as '*' operator between two vector3d
     * @param a the input vector3d
     */
    inline double odot(const Vector3d &a)const{
        // Thanks Qingchen&Jie for pointing out the Floating-point underflow issue
        // Generally using double is enough.
        double sum = 0.0;
        if(m_Vals[0] >= std::numeric_limits<double>::epsilon()/a.m_Vals[0]||
           a.m_Vals[0] >= std::numeric_limits<double>::epsilon()/m_Vals[0]) {
            sum += m_Vals[0] * a.m_Vals[0];
        }
        if(m_Vals[1] >= std::numeric_limits<double>::epsilon()/a.m_Vals[1] ||
           a.m_Vals[1] >= std::numeric_limits<double>::epsilon()/m_Vals[1]){
            sum += m_Vals[1] * a.m_Vals[1];
        } 
        if(m_Vals[2] >= std::numeric_limits<double>::epsilon()/a.m_Vals[2]||
           a.m_Vals[2] >= std::numeric_limits<double>::epsilon()/m_Vals[2]){
            sum += m_Vals[2] * a.m_Vals[2];
        } 
        return sum;
    }
    //*************************************************
    /**
     * *= operator
     * @param val right hand side double value
     */
    inline Vector3d& operator*=(const double &val){
        m_Vals[0]*=val;m_Vals[1]*=val;m_Vals[2]*=val;
        return *this;
    }
    //***********************************************
    /**
     * / operator
     * @param val right hand side double value
     */
    inline Vector3d operator/(const double &val)const{
        Vector3d temp(0.0);
        if(abs(val)<1.0e-15){
            MessagePrinter::printErrorTxt("val= "+to_string(val)+" is singular for / operator in Vector3");
            MessagePrinter::exitAsFem();
        }
        temp.m_Vals[0]=m_Vals[0]/val;temp.m_Vals[1]=m_Vals[1]/val;temp.m_Vals[2]=m_Vals[2]/val;
        return temp;
    }
    //*****************************************************
    //*** for other math funs
    //*****************************************************
    /**
     * return the L2 norm of current vector3 array
     */
    inline double norm()const{
        double sum=static_cast<double>(m_Vals[0]*m_Vals[0]+m_Vals[1]*m_Vals[1]+m_Vals[2]*m_Vals[2]);
        return sqrt(sum);
    }
    /**
     * return the squared norm of current vector3 array
     */
    inline double normsq()const{
        double sum=static_cast<double>(m_Vals[0]*m_Vals[0]+m_Vals[1]*m_Vals[1]+m_Vals[2]*m_Vals[2]);
        return sum;
    }


private:
    double m_Vals[3];/**< components vector, its size is always 3! */
};