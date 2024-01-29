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
//+++ Date   : 2022.06.13
//+++ Purpose: implement the API for PETSc and other package's sparse
//+++          matrix, for small size matrix please use MatrixXd
//+++          this is mainly used for the system K matrix
//+++ TODO   : in the future, we should keep the same function name
//+++          but offer APIs to different packages, i.e., trilinos.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "petsc.h"
#include "Utils/MessagePrinter.h"
#include "DofHandler/DofHandler.h"


/**
 * This class implements the general operations for sparse matrix
 */
class SparseMatrix{
public:
    /**
     * constructor for different purpose
     */
    SparseMatrix();
    /**
     * constructor with sparse matrix
     */
    SparseMatrix(const SparseMatrix &a);
    /**
     * init from dof handler
     * @param t_dofHandler the dofHandler class
     */
    void init(const DofHandler &t_dofHandler);
    /**
     * resize the matrix, this will re-allocate memory and destroy the previously defined structure
     * @param m integer for the 1st dimension
     * @param n integer for the 2nd dimension
     * @param maxrownnz integer for the maximum non-zero elements of each row
     */
    inline void resize(const int &m,const int &n,const int &maxrownnz){
        if(m<0||n<0||m!=n){
            MessagePrinter::printErrorTxt("either you m<0 or n<0 or m!=n detected in resize function");
            MessagePrinter::exitAsFem();
        }
        if(m_allocated){
            MatDestroy(&m_matrix);
            MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,m,n,maxrownnz,NULL,maxrownnz,NULL,&m_matrix);
            m_m=m;m_n=n;
            m_allocated=true;
        }
        else{
            MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,m,n,maxrownnz,NULL,maxrownnz,NULL,&m_matrix);
            m_m=m;m_n=n;
            m_allocated=true;
        }
        MatSetOption(m_matrix,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);//allow new element insertion
    }
    //****************************************************************
    //*** operators
    //****************************************************************
    /**
     * = operator for current sparse matrix, if you want to set matrix to zero, please use setToZero,
     * which is much faster than =0 operator!
     * @param a right hand side scalar
     */
    SparseMatrix& operator=(const double &a);
    /**
     * = operator for current sparse matrix
     * @param a right hand side sparse matrix
     */
    SparseMatrix& operator=(const SparseMatrix &a);
    //******************************************
    /**
     * + operator for current sparse matrix
     * @param a right hand side double scalar
     */
    SparseMatrix operator+(const double &a)const;
    /**
     * + operator for current sparse matrix
     * @param a right hand side sparse matrix
     */
    SparseMatrix operator+(const SparseMatrix &a)const;
    /**
     * += operator for current sparse matrix
     * @param a right hand side double scalar
     */
    SparseMatrix& operator+=(const double &a);
    /**
     * += operator for current sparse matrix
     * @param a right hand side sparse matrix
     */
    SparseMatrix& operator+=(const SparseMatrix &a);
    //******************************************
    /**
     * - operator for current sparse matrix
     * @param a right hand side double scalar
     */
    SparseMatrix operator-(const double &a)const;
    /**
     * - operator for current sparse matrix
     * @param a right hand side sparse matrix
     */
    SparseMatrix operator-(const SparseMatrix &a)const;
    /**
     * -= operator for current sparse matrix
     * @param a right hand side double scalar
     */
    SparseMatrix& operator-=(const double &a);
    /**
     * -= operator for current sparse matrix
     * @param a right hand side sparse matrix
     */
    SparseMatrix& operator-=(const SparseMatrix &a);
    //***********************************************
    /**
     * * operator for current sparse matrix
     * @param a right hand side double scalar
     */
    SparseMatrix operator*(const double &a)const;
    /**
     * *= operator for current sparse matrix
     * @param a right hand side scalar
     */
    SparseMatrix& operator*=(const double &a);
    //***********************************************
    /**
     * / operator for current sparse matrix
     * @param a right hand side double scalar
     */
    SparseMatrix operator/(const double &a)const;
    /**
     * /= operator for current sparse matrix
     * @param a right hand side double scalar
     */
    SparseMatrix& operator/=(const double &a);


    /**
     * set all the elements to zero but keep the sparsity pattern
     */
    inline void setToZero(){
        MatZeroEntries(m_matrix);
    }
    /**
     * set all elements to random value but keep the sparsity pattern
     */
    inline void setToRandom(){
        PetscRandom    rand;
        PetscRandomCreate(PETSC_COMM_WORLD, &rand);
        PetscRandomSetFromOptions(rand);
        MatSetRandom(m_matrix,rand);
        PetscRandomDestroy(&rand);
    }
    inline void disableReallocation(){
        MatSetOption(m_matrix,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);//disable new element insertion
    }
    //****************************************************************
    //*** add/insert value operations
    //****************************************************************
    /**
     * add single value to matrix
     * @param i integer for 1st dim, start from 1
     * @param j integer for 2nd dim, start from 1
     * @param val the double type scalar
     */
    inline void addValue(const int &i,const int &j,const double &val){
        MatSetValue(m_matrix,i-1,j-1,val,ADD_VALUES);
    }
    /**
     * add values to matrix
     * @param iInd integer index vector for 1st dim, start from 0
     * @param jInd integer index vector for 2nd dim, start from 0
     * @param val the double type scalar vector
     */
    inline void addValues(const vector<int> &iInd,const vector<int> &jInd,const vector<double> &vals){
        MatSetValues(m_matrix,iInd.size(),iInd.data(),jInd.size(),jInd.data(),vals.data(),ADD_VALUES);
    }
    /**
     * add values to matrix with predefined size
     * @param isize the size of 1st vector to be assembled
     * @param iInd integer index vector for 1st dim, start from 0
     * @param jsize the size of 2nd vector to be assembled
     * @param jInd integer index vector for 2nd dim, start from 0
     * @param val the double type scalar vector
     */
    inline void addValues(const int &isize,const vector<int> &iInd,const int &jsize,const vector<int> &jInd,const vector<double> &vals){
        MatSetValues(m_matrix,isize,iInd.data(),jsize,jInd.data(),vals.data(),ADD_VALUES);
    }
    /**
     * add values to matrix with predefined size
     * @param isize the size of 1st vector to be assembled
     * @param iInd integer index vector for 1st dim, start from 0
     * @param jsize the size of 2nd vector to be assembled
     * @param jInd integer index vector for 2nd dim, start from 0
     * @param val the double type scalar vector
     */
    inline void addValues(const int &isize,const int (&iInd)[],const int &jsize,const int (&jInd)[],const vector<double> &vals){
        MatSetValues(m_matrix,isize,iInd,jsize,jInd,vals.data(),ADD_VALUES);
    }
    /**
     * add values to matrix with predefined size
     * @param isize the size of 1st vector to be assembled
     * @param iInd integer index vector for 1st dim, start from 0
     * @param jsize the size of 2nd vector to be assembled
     * @param jInd integer index vector for 2nd dim, start from 0
     * @param val the double type scalar vector
     */
    inline void addValues(const int &isize,const int *iInd,const int &jsize,const int *jInd,const double *vals){
        MatSetValues(m_matrix,isize,iInd,jsize,jInd,vals,ADD_VALUES);
    }
    //************************************************
    /**
     * insert single value to matrix
     * @param i integer for 1st dim, start from 1
     * @param j integer for 2nd dim, start from 1
     * @param val the double type scalar
     */
    inline void insertValue(const int &i,const int &j,const double &val){
        MatSetValue(m_matrix,i-1,j-1,val,INSERT_VALUES);
    }
    /**
     * insert values to matrix
     * @param iInd integer index vector for 1st dim, start from 0
     * @param jInd integer index vector for 2nd dim, start from 0
     * @param val the double type scalar vector
     */
    inline void insertValues(const vector<int> &iInd,const vector<int> &jInd,const vector<double> &vals){
        MatSetValues(m_matrix,iInd.size(),iInd.data(),jInd.size(),jInd.data(),vals.data(),INSERT_VALUES);
    }
    /**
     * insert values to matrix with predefined size
     * @param isize the size of 1st vector to be assembled
     * @param iInd integer index vector for 1st dim, start from 0
     * @param jsize the size of 2nd vector to be assembled
     * @param jInd integer index vector for 2nd dim, start from 0
     * @param val the double type scalar vector
     */
    inline void insertValues(const int &isize,const vector<int> &iInd,const int &jsize,const vector<int> &jInd,const vector<double> &vals){
        MatSetValues(m_matrix,isize,iInd.data(),jsize,jInd.data(),vals.data(),INSERT_VALUES);
    }
    /**
     * insert values to matrix with predefined size
     * @param isize the size of 1st vector to be assembled
     * @param iInd integer index vector for 1st dim, start from 0
     * @param jsize the size of 2nd vector to be assembled
     * @param jInd integer index vector for 2nd dim, start from 0
     * @param val the double type scalar vector
     */
    inline void insertValues(const int &isize,const int (&iInd)[],const int &jsize,const int (&jInd)[],const vector<double> &vals){
        MatSetValues(m_matrix,isize,iInd,jsize,jInd,vals.data(),INSERT_VALUES);
    }
    //********************************************************
    /**
     * assemble the sparse matrix, this should be called after all the intermediate assemble
     */
    inline void assemble(){
        MatAssemblyBegin(m_matrix,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(m_matrix,MAT_FINAL_ASSEMBLY);
    }
    //****************************************************************
    //*** general gettings
    //****************************************************************
    /**
     * get the size of current sparse matrix(return 1st dim, 1st dim=2nd dim)
     */
    inline int getSize()const{return m_m;}
    /**
     * get the L2 norm of current sparse matrix
     */
    inline double getNorm()const{
        double norm=0.0;
        MatNorm(m_matrix,NORM_FROBENIUS,&norm);
        return norm;
    }
    /**
     * get the reference of current matrix
     */
    inline Mat& getReference(){return m_matrix;}
    /**
     * get the copy of current matrix
     */
    inline Mat getCopy()const{return m_matrix;}
    /**
     * copy current SparseMatrix to PETSc Mat
     * @param a the PETSc Mat class
     */
    inline void copy2Mat(Mat &a){
        int m,n;
        MatGetSize(a,&m,&n);
        if(m!=m_m || n!=m_n){
            MessagePrinter::printErrorTxt("can\'t copy current sparse matrix to PETSc Mat, their sizes do not match");
            MessagePrinter::exitAsFem();
        }
        MatCopy(m_matrix,a,SAME_NONZERO_PATTERN);
        MatSetOption(a,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);// disable reallocation
    }
    /**
     * release the memory and destroy mat object
     */
    inline void releaseMemory(){
        if(m_allocated){
            MatDestroy(&m_matrix);
            m_allocated=false;
            m_m=m_n=0;
        }
    }

    /**
     * print out the sparse matrix
     * @param txt preset text for terminal print
     */
    void printMatrix(const string &txt="")const;

private:
    bool m_allocated=false;/**< boolean for allocation status */
    Mat m_matrix;/**< sparse matrix class */
    int m_m;/**< the size of 1st dim */
    int m_n;/**< the size of 2nd dim */
};