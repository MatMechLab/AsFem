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
//+++ Date   : 2022.06.04
//+++ Purpose: implement the API for PETSc and other package's vector
//+++          for dense vector or small size vector, please use
//+++          VectorXd, this is mainly used for the system solution
//+++          and system residual
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "petsc.h"
#include "Utils/MessagePrinter.h"

/**
 * This class implements the general vector class with API to PETSc and other packages. This vector
 * should be used for global right hand side vector(with large size)
 */
class Vector{
public:
    /**
     * constructor
     */
    Vector();
    /**
     * constructor
     * @param n integer for the size of current vector
     */
    Vector(const int &n);
    /**
     * constructor
     * @param n integer for the size of current vector
     * @param val double value for the initial value
     */
    Vector(const int &n,const double &val);
    /**
     * constructor
     * @param a vector class
     */
    Vector(const Vector &a);

    /**
     * allocate memory for current vector, this will re-allocate memory and remove previous values
     */
    void setup();
    /**
     * resize the vector, this will re-allocate memory and destroy the old structure
     * @param n integer for the size of current vector
     */
    void resize(const int &n);
    /**
     * resize the vector with given value
     * @param n integer for the size of current vector
     * @param val the given initial value
     */
    void resize(const int &n,const double &val);
    //********************************************************
    //*** operators
    //********************************************************
    /**
     * get the i-th component of current vecor
     * @param i i-th index, start from 1
     */
    inline double operator()(const int &i)const{
        double val;int index[1];
        index[0]=i-1;
        VecGetValues(m_vector,1,index,&val);
        return val;
    }
    /**
     * = operator for vector with scalar right hand side value
     * @param val the right hand side scalar value
     */
    inline Vector& operator=(const double &val){
        if(!m_allocated){
            MessagePrinter::printErrorTxt("can\'t apply = to an empty vector, you must initialize it before use");
            MessagePrinter::exitAsFem();
        }
        VecSet(m_vector,val);
        assemble();
        return *this;
    }
    /**
     * = operator for vectors with vector type right hand side
     * @param a the right hand side vector class
     */
    inline Vector& operator=(const Vector &a){
        if(!m_allocated){
            MessagePrinter::printErrorTxt("can\'t apply = to an empty vector, you must initialize it before use");
            MessagePrinter::exitAsFem();
        }
        if(m_size!=a.getSize()){
            MessagePrinter::printErrorTxt("can\'t apply = to two different vectors, they must have the same size");
            MessagePrinter::exitAsFem();
        }
        VecCopy(a.m_vector,m_vector);
        assemble();
        return *this;
    }
    /**
     * = operator for vectors with vec type right hand side
     * @param a the right hand side vector class
     */
    inline Vector& operator=(const Vec &a){
        if(!m_allocated){
            MessagePrinter::printErrorTxt("can\'t apply = to an empty vector, you must initialize it before use");
            MessagePrinter::exitAsFem();
        }
        int size;
        VecGetSize(a,&size);
        if(m_size!=size){
            MessagePrinter::printErrorTxt("can\'t apply = to two different vectors, they must have the same size");
            MessagePrinter::exitAsFem();
        }
        VecCopy(a,m_vector);
        assemble();
        return *this;
    }
    /**
     * copy entities from another Vector class
     * @param a the vector class
    */
    inline void copyFrom(const Vector &a){
        if(!m_allocated){
            MessagePrinter::printErrorTxt("can\'t apply = to an empty vector, you must initialize it before use");
            MessagePrinter::exitAsFem();
        }
        if(m_size!=a.getSize()){
            MessagePrinter::printErrorTxt("can\'t apply = to two different vectors, they must have the same size");
            MessagePrinter::exitAsFem();
        }
        VecCopy(a.m_vector,m_vector);
        assemble();
    }
    /**
     * copy entities from PETSc vec
     * @param a the PETSc vec object
    */
    inline void copyFrom(const Vec &a){
        if(!m_allocated){
            MessagePrinter::printErrorTxt("can\'t apply copyFrom to an empty vector, you must initialize it before use");
            MessagePrinter::exitAsFem();
        }
        int size;
        VecGetSize(a,&size);
        if(m_size!=size){
            MessagePrinter::printErrorTxt("can\'t apply copyFrom to two different vectors, they must have the same size");
            MessagePrinter::exitAsFem();
        }
        VecCopy(a,m_vector);
        assemble();
    }
    //**********************************************
    /**
     * + operator for vector
     * @param val the right hand side scalar
     */
    inline Vector operator+(const double &val)const{
        Vector anew;
        VecDuplicate(m_vector,&anew.m_vector);
        VecCopy(m_vector,anew.m_vector);
        anew.m_size=m_size;
        anew.m_allocated=true;
        VecShift(anew.m_vector,val);
        VecAssemblyBegin(anew.m_vector);
        VecAssemblyEnd(anew.m_vector);
        return anew;
    }
    /**
     * + operator for vector
     * @param a right hand side vector
     */
    inline Vector operator+(const Vector &a)const{
        if(a.getSize()!=m_size){
            MessagePrinter::printErrorTxt("can\'t appy + to two different vectors, they must have the same size");
            MessagePrinter::exitAsFem(); 
        }
        Vector anew;
        VecDuplicate(m_vector,&anew.m_vector);
        VecCopy(m_vector,anew.m_vector);
        anew.m_size=m_size;
        anew.m_allocated=true;
        //VecAXPY(Vec y,PetscScalar a,Vec x);//y = y + a*x
        VecAXPY(anew.m_vector,1.0,a.m_vector);
        VecAssemblyBegin(anew.m_vector);
        VecAssemblyEnd(anew.m_vector);
        return anew;
    }
    /**
     * + operator for vec
     * @param a right hand side vector
     */
    inline Vector operator+(const Vec &a)const{
        int size;
        VecGetSize(a,&size);
        if(m_size!=size){
            MessagePrinter::printErrorTxt("can\'t appy + to two different vectors, they must have the same size");
            MessagePrinter::exitAsFem(); 
        }
        Vector anew;
        VecDuplicate(m_vector,&anew.m_vector);
        VecCopy(m_vector,anew.m_vector);
        anew.m_size=m_size;
        anew.m_allocated=true;
        //VecAXPY(Vec y,PetscScalar a,Vec x);//y = y + a*x
        VecAXPY(anew.m_vector,1.0,a);
        VecAssemblyBegin(anew.m_vector);
        VecAssemblyEnd(anew.m_vector);
        return anew;
    }
    //***********************************
    /**
     * += operator for vector
     * @param val the right hand side scalar
     */
    inline Vector& operator+=(const double &val){
        VecShift(m_vector,val);
        assemble();
        return *this;
    }
    /**
     * += operator for vector
     * @param a the right hand side vector class
     */
    inline Vector& operator+=(const Vector &a){
        if(m_size!=a.getSize()){
            MessagePrinter::printErrorTxt("can\'t apply += to two different vectors, they must have the same size");
            MessagePrinter::exitAsFem();
        }
        VecAXPY(m_vector,1.0,a.m_vector);
        assemble();
        return *this;
    }
    /**
     * += operator for Vec
     * @param a the right hand side vector class
     */
    inline Vector& operator+=(const Vec &a){
        int size;
        VecGetSize(a,&size);
        if(m_size!=size){
            MessagePrinter::printErrorTxt("can\'t apply += to two different vectors, they must have the same size");
            MessagePrinter::exitAsFem();
        }
        VecAXPY(m_vector,1.0,a);
        assemble();
        return *this;
    }
    //**********************************************
    /**
     * - operator for vector
     * @param val the right hand side scalar
     */
    inline Vector operator-(const double &val)const{
        Vector anew;
        VecDuplicate(m_vector,&anew.m_vector);
        VecCopy(m_vector,anew.m_vector);
        anew.m_size=m_size;
        anew.m_allocated=true;
        VecShift(anew.m_vector,-1.0*val);
        VecAssemblyBegin(anew.m_vector);
        VecAssemblyEnd(anew.m_vector);
        return anew;
    }
    /**
     * - operator for vector
     * @param a right hand side vector class 
     */
    inline Vector operator-(const Vector &a)const{
        if(a.getSize()!=m_size){
            MessagePrinter::printErrorTxt("can\'t appy - to two different vectors, they must have the same size");
            MessagePrinter::exitAsFem(); 
        }
        Vector anew;
        VecDuplicate(m_vector,&anew.m_vector);
        VecCopy(m_vector,anew.m_vector);
        anew.m_size=m_size;
        anew.m_allocated=true;
        //VecAXPY(Vec y,PetscScalar a,Vec x);//y = y + a*x
        VecAXPY(anew.m_vector,-1.0,a.m_vector);
        VecAssemblyBegin(anew.m_vector);
        VecAssemblyEnd(anew.m_vector);
        return anew;
    }
    /**
     * - operator for Vec
     * @param a right hand side vector class 
     */
    inline Vector operator-(const Vec &a)const{
        int size;
        VecGetSize(a,&size);
        if(size!=m_size){
            MessagePrinter::printErrorTxt("can\'t appy - to two different vectors, they must have the same size");
            MessagePrinter::exitAsFem(); 
        }
        Vector anew;
        VecDuplicate(m_vector,&anew.m_vector);
        VecCopy(m_vector,anew.m_vector);
        anew.m_size=m_size;
        anew.m_allocated=true;
        //VecAXPY(Vec y,PetscScalar a,Vec x);//y = y + a*x
        VecAXPY(anew.m_vector,-1.0,a);
        VecAssemblyBegin(anew.m_vector);
        VecAssemblyEnd(anew.m_vector);
        return anew;
    }
    //***********************************
    /**
     * -= operator for vector
     * @param val the right hand side scalar
     */
    inline Vector& operator-=(const double &val){
        VecShift(m_vector,-1.0*val);
        assemble();
        return *this;
    }
    /**
     * -= operator for vector
     * @param a the right hand side vector class
     */
    inline Vector& operator-=(const Vector &a){
        if(m_size!=a.getSize()){
            MessagePrinter::printErrorTxt("can\'t apply -= to two different vectors, they must have the same size");
            MessagePrinter::exitAsFem();
        }
        VecAXPY(m_vector,-1.0,a.m_vector);
        assemble();
        return *this;
    }
    /**
     * -= operator for Vec
     * @param a the right hand side vector class
     */
    inline Vector& operator-=(const Vec &a){
        int size;
        VecGetSize(a,&size);
        if(m_size!=size){
            MessagePrinter::printErrorTxt("can\'t apply -= to two different vectors, they must have the same size");
            MessagePrinter::exitAsFem();
        }
        VecAXPY(m_vector,-1.0,a);
        assemble();
        return *this;
    }
    //**********************************************
    /**
     * * operator for vector
     * @param val the right hand side scalar
     */
    inline Vector operator*(const double &val)const{
        Vector anew;
        VecDuplicate(m_vector,&anew.m_vector);
        VecCopy(m_vector,anew.m_vector);
        anew.m_size=m_size;
        anew.m_allocated=true;
        VecScale(anew.m_vector,val);
        VecAssemblyBegin(anew.m_vector);
        VecAssemblyEnd(anew.m_vector);
        return anew;
    }
    /**
     * dot for double contruction between two vectors 
     * @param a right hand side vector class 
     */
    inline double dot(const Vector &a)const{
        if(a.getSize()!=m_size){
            MessagePrinter::printErrorTxt("can\'t appy dot to two different vectors, they must have the same size");
            MessagePrinter::exitAsFem(); 
        }
        double val;
        VecDot(m_vector,a.m_vector,&val);
        return val;
    }
    /**
     * dot for double contruction between two vectors 
     * @param a right hand side vector class 
     */
    inline double dot(const Vec &a)const{
        int size;
        VecGetSize(a,&size);
        if(size!=m_size){
            MessagePrinter::printErrorTxt("can\'t appy dot to two different vectors, they must have the same size");
            MessagePrinter::exitAsFem(); 
        }
        double val;
        VecDot(m_vector,a,&val);
        return val;
    }
    //***********************************
    /**
     * *= operator for vector
     * @param val the right hand side scalar
     */
    inline Vector& operator*=(const double &val){
        VecScale(m_vector,val);
        assemble();
        return *this;
    }
    //**********************************************
    /**
     * / operator for vector
     * @param val the right hand side scalar
     */
    inline Vector operator/(const double &val)const{
        if(abs(val)<1.0e-15){
            MessagePrinter::printErrorTxt("can\'t appy / to a singular value, val must be nonzero");
            MessagePrinter::exitAsFem(); 
        }
        Vector anew;
        VecDuplicate(m_vector,&anew.m_vector);
        VecCopy(m_vector,anew.m_vector);
        anew.m_size=m_size;
        anew.m_allocated=true;
        VecScale(anew.m_vector,1.0/val);
        VecAssemblyBegin(anew.m_vector);
        VecAssemblyEnd(anew.m_vector);
        return anew;
    }
    //***********************************
    /**
     * /= operator for vector
     * @param val the right hand side scalar
     */
    inline Vector& operator/=(const double &val){
        if(abs(val)<1.0e-15){
            MessagePrinter::printErrorTxt("can\'t appy /= to a singular value, val must be nonzero");
            MessagePrinter::exitAsFem(); 
        }
        VecScale(m_vector,1.0/val);
        assemble();
        return *this;
    }

    //********************************************************
    //*** general settings
    //********************************************************
    /**
     * set the size of current vector
     */
    inline void setSize(const int &n){m_size=n;}
    /**
     * set the current vector's elements to zero, but keep the storage
     */
    inline void setToZero(){
        if(m_size) {
            VecSet(m_vector,0.0);
            assemble();
        }
    }
    /**
     * set the current vector's elements to random value, but keep the storage
     */
    inline void setToRandom(){
        PetscRandom    rand;
        PetscRandomCreate(PETSC_COMM_WORLD, &rand);
        PetscRandomSetFromOptions(rand);
        VecSetRandom(m_vector,rand);
        PetscRandomDestroy(&rand);
        assemble();
    }
    /**
     * add single value to current vector
     * @param i i-th index, start form 1
     * @param val scalar value
     */
    inline void addValue(const int &i,const double &val){
        VecSetValue(m_vector,i-1,val,ADD_VALUES);
    }
    /**
     * add multiple values to current vector
     * @param n the number of elements to be added
     * @param index the integer vector(it should start from 0)
     * @param vals the values vector
     */
    inline void addValues(const int &n,const int *index,const double *vals){
        VecSetValues(m_vector,n,index,vals,ADD_VALUES);
    }
    /**
     * add multiple values to current vector
     * @param index the integer index vector(it should start from 0)
     * @param vals the values vector(it must have the same size as index)
     */
    inline void addValues(const vector<int> &index,const vector<double> &vals){
        VecSetValues(m_vector,index.size(),index.data(),vals.data(),ADD_VALUES);
    }
    //****************************************************************
    /**
     * insert single value to current vector
     * @param i i-th index, start form 1
     * @param val scalar value
     */
    inline void insertValue(const int &i,const double &val){
        VecSetValue(m_vector,i-1,val,INSERT_VALUES);
    }
    /**
     * insert multiple values to current vector
     * @param n the number of elements to be added
     * @param index the integer vector(it should start from 0)
     * @param vals the values vector
     */
    inline void insertValues(const int &n,const int *index,const double *vals){
        VecSetValues(m_vector,n,index,vals,INSERT_VALUES);
    }
    /**
     * insert multiple values to current vector
     * @param index the integer index vector(it should start from 0)
     * @param vals the values vector(it must have the same size as index)
     */
    inline void insertValues(const vector<int> &index,const vector<double> &vals){
        VecSetValues(m_vector,index.size(),index.data(),vals.data(),INSERT_VALUES);
    }
    /**
     * assemble function, you can add/insert values as you like, but do not forget to do the
     * 'assemble' call, otherwise, no value can be added into your vector
     */
    inline void assemble(){
        if(m_size<1) return;
        VecAssemblyBegin(m_vector);
        VecAssemblyEnd(m_vector);
    }

    //********************************************************
    //*** general gettings
    //********************************************************
    /**
     * get the size of current vector
     */
    inline int getSize()const{return m_size;}
    /**
     * get the L2 norm of current vector
     */
    inline double getNorm()const{
        double val;
        VecNorm(m_vector,NORM_2,&val);
        return val;
    }
    /**
     * get the sum of current vector
     */
    inline double getSum()const{
        double val;
        VecSum(m_vector,&val);
        return val;
    }
    /**
     * get the reference of vector 
     */
    inline Vec& getVectorRef(){return m_vector;}
    /**
     * get the reference of vector 
     */
    inline Vec  getVectorCopy()const{return m_vector;}
    /**
     * copy current vector to PETSc vec
     * @param a the PETSc Vec class
     */
    inline void copy2Vec(Vec &a){
        int size;
        VecGetSize(a,&size);
        if(size!=m_size){
            MessagePrinter::printErrorTxt("can\'t copy current vector to PETSc vec, the size does not match");
            MessagePrinter::exitAsFem();
        }
        VecCopy(m_vector,a);
        VecAssemblyBegin(a);
        VecAssemblyEnd(a);
    }
    /**
     * de-allocate the vector and release the memory
     */
    void releaseMemory(){
        if(m_allocated) VecDestroy(&m_vector);
        m_size=0;
        m_allocated=false;
        if(m_ghostallocated){
            VecScatterDestroy(&m_scatter);
            VecDestroy(&m_vector_ghost);
            m_ghostallocated=false;
        }
    }
    /**
     * print the vector's elements
     * @param txt preset content for terminal output
     */
    void printVec(const string &txt="")const;

    //********************************************************
    //*** for ghost access(for multi-core case)
    //********************************************************
    /**
     * make a ghost copy of current vector for the cross-core access, after use, one must
     * call destroyGhostCopy !!!
     */
    void makeGhostCopy();
    /**
     * destroy the ghost copy of current vector and release the allocated memory
     */
    void destroyGhostCopy();
    /**
     * get the i-th element value from a ghost copied vector
     * @param i integer for element index
     */
    inline double getIthValueFromGhost(const int &i)const{
        if(i<1||i>getSize()){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of range("+to_string(getSize())+") for your vector");
            MessagePrinter::exitAsFem();
        }
        double val;int index[1];
        index[0]=i-1;
        VecGetValues(m_vector_ghost,1,index,&val);
        return val;
    }


private:
    bool m_allocated=false;/**< boolean flag for the status of allocation */
    Vec m_vector;/** petsc vector for current vector class */
    int m_size;/** the size of current vector */

    bool m_ghostallocated=false;/**< boolean flag for ghost processor allocation */
    Vec m_vector_ghost;/**< the ghost copy of m_vector, which should be destroy after use */
    VecScatter m_scatter;/**< scatter used for ghost copy */
};
