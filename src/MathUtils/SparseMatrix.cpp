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
//+++ Date   : 2022.06.13
//+++ Purpose: implement the API for PETSc and other package's sparse
//+++          matrix, for small size matrix please use MatrixXd
//+++          this is mainly used for the system K matrix
//+++ TODO   : in the future, we should keep the same function name
//+++          but offer APIs to different packages, i.e., trilinos.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MathUtils/SparseMatrix.h"

SparseMatrix::SparseMatrix(){
    m_m=0;m_n=0;
    m_allocated=false;
    m_matrix=NULL;
}
SparseMatrix::SparseMatrix(const SparseMatrix &a){
    m_m=a.m_m;
    m_n=a.m_n;
    MatDuplicate(a.m_matrix,MAT_SHARE_NONZERO_PATTERN,&m_matrix);
    MatCopy(a.m_matrix,m_matrix,SAME_NONZERO_PATTERN);
    m_allocated=true;
}
//************************************************************
//*** operator overload for different purpose
//************************************************************
// for = operator
SparseMatrix& SparseMatrix::operator=(const double &a){
    // if you want to set your matrix to zero(but keep the structure), please use setToZero, 
    // which is faster than this one!!!
    if(!m_allocated){
        MessagePrinter::printErrorTxt("can\'t apply = to your sparse matrix, it is not allocated");
        MessagePrinter::exitAsFem();
    }
    int start,end,r,j;
    int ncols;
    const int *cols;
    const double *vals;
    Mat temp;
    MatDuplicate(m_matrix,MAT_SHARE_NONZERO_PATTERN,&temp);
    MatZeroEntries(temp);
    MatGetOwnershipRange(m_matrix,&start,&end);
    for(r=start;r<end;r++){
        MatGetRow(m_matrix,r,&ncols,&cols,&vals);
        for(j=0;j<ncols;j++){
            MatSetValues(temp,1,&r,1,&cols[j],&a,INSERT_VALUES);
        }
        MatRestoreRow(m_matrix,r,&ncols,&cols,&vals);
    }
    MatAssemblyBegin(temp,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(temp,MAT_FINAL_ASSEMBLY);

    MatCopy(temp,m_matrix,SAME_NONZERO_PATTERN);
    MatDestroy(&temp);

    return *this;
}
SparseMatrix& SparseMatrix::operator=(const SparseMatrix &a){
    if(!m_allocated){
        m_m=a.m_m;
        m_n=a.m_n;
        MatDuplicate(a.m_matrix,MAT_SHARE_NONZERO_PATTERN,&m_matrix);
        MatCopy(a.m_matrix,m_matrix,SAME_NONZERO_PATTERN);
        m_allocated=true;
    }
    else{
        if(m_m!=a.m_m || m_n!=a.m_n){
            MessagePrinter::printErrorTxt("can\'t apply = to sparse matrix with different size");
            MessagePrinter::exitAsFem();
        }
        MatCopy(a.m_matrix,m_matrix,SAME_NONZERO_PATTERN);
    }
    return *this;
}
//*****************************
//*** + and += operator
//*****************************
SparseMatrix SparseMatrix::operator+(const double &a)const{
    SparseMatrix temp;
    temp=*this;
    int start,end,r,j;
    int ncols;
    const int *cols;
    MatGetOwnershipRange(m_matrix,&start,&end);
    for(r=start;r<end;r++){
        MatGetRow(m_matrix,r,&ncols,&cols,NULL);
        for(j=0;j<ncols;j++){
            MatSetValues(temp.m_matrix,1,&r,1,&cols[j],&a,ADD_VALUES);
        }
        MatRestoreRow(m_matrix,r,&ncols,&cols,NULL);
    }
    MatAssemblyBegin(temp.m_matrix,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(temp.m_matrix,MAT_FINAL_ASSEMBLY);
    return temp;
}
SparseMatrix SparseMatrix::operator+(const SparseMatrix &a)const{
    if(a.m_m!=m_m || a.m_n!=m_n){
        MessagePrinter::printErrorTxt("can\'t apply + to sparse matrix with different size");
        MessagePrinter::exitAsFem();
    }
    SparseMatrix temp;
    temp=*this;
    //Computes Y = a*X + Y.
    MatAXPY(temp.m_matrix,1.0,a.m_matrix,SAME_NONZERO_PATTERN);
    return temp;
}
SparseMatrix& SparseMatrix::operator+=(const double &a){
    int start,end,r,j;
    int ncols;
    const int *cols;
    Mat temp=NULL;
    MatDuplicate(m_matrix,MAT_SHARE_NONZERO_PATTERN,&temp);
    MatCopy(m_matrix,temp,SAME_NONZERO_PATTERN);
    MatGetOwnershipRange(m_matrix,&start,&end);
    for(r=start;r<end;r++){
        MatGetRow(m_matrix,r,&ncols,&cols,NULL);
        for(j=0;j<ncols;j++){
            // MatSetValues(m_matrix,1,&r,1,&cols[j],&a,ADD_VALUES);
            MatSetValue(temp,r,cols[j],a,ADD_VALUES);
        }
        MatRestoreRow(m_matrix,r,&ncols,&cols,NULL);
    }
    MatAssemblyBegin(temp,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(temp,MAT_FINAL_ASSEMBLY);
    MatCopy(temp,m_matrix,SAME_NONZERO_PATTERN);
    MatDestroy(&temp);
    return *this;
}
SparseMatrix& SparseMatrix::operator+=(const SparseMatrix &a){
    if(a.m_m!=m_m || a.m_n!=m_n){
        MessagePrinter::printErrorTxt("can\'t apply += to sparse matrix with different size");
        MessagePrinter::exitAsFem();
    }
    //Computes Y = a*X + Y.
    MatAXPY(m_matrix,1.0,a.m_matrix,SAME_NONZERO_PATTERN);
    return *this;
}
//*****************************
//*** - and -= operator
//*****************************
SparseMatrix SparseMatrix::operator-(const double &a)const{
    SparseMatrix temp;
    double aa;
    int start,end,r,j;
    int ncols;
    const int *cols;
    const double *vals;
    aa=-1.0*a;
    temp.m_m=m_m;
    temp.m_n=m_n;
    MatDuplicate(m_matrix,MAT_SHARE_NONZERO_PATTERN,&temp.m_matrix);
    MatCopy(m_matrix,temp.m_matrix,SAME_NONZERO_PATTERN);
    temp.m_allocated=true;
    MatGetOwnershipRange(m_matrix,&start,&end);
    for(r=start;r<end;r++){
        MatGetRow(m_matrix,r,&ncols,&cols,&vals);
        for(j=0;j<ncols;j++){
            MatSetValues(temp.m_matrix,1,&r,1,&cols[j],&aa,ADD_VALUES);
        }
        MatRestoreRow(m_matrix,r,&ncols,&cols,&vals);
    }
    MatAssemblyBegin(temp.m_matrix,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(temp.m_matrix,MAT_FINAL_ASSEMBLY);
    return temp;
}
SparseMatrix SparseMatrix::operator-(const SparseMatrix &a)const{
    if(a.m_m!=m_m || a.m_n!=m_n){
        MessagePrinter::printErrorTxt("can\'t apply - to sparse matrix with different size");
        MessagePrinter::exitAsFem();
    }
    SparseMatrix temp;
    temp=*this;
    //Computes Y = a*X + Y.
    MatAXPY(temp.m_matrix,-1.0,a.m_matrix,SAME_NONZERO_PATTERN);
    return temp;
}
SparseMatrix& SparseMatrix::operator-=(const double &a){
    int start,end,r,j;
    int ncols;
    const int *cols;
    const double *vals;
    double aa;
    aa=-1.0*a;
    Mat temp;
    MatDuplicate(m_matrix,MAT_SHARE_NONZERO_PATTERN,&temp);
    MatCopy(m_matrix,temp,SAME_NONZERO_PATTERN);
    MatGetOwnershipRange(m_matrix,&start,&end);
    for(r=start;r<end;r++){
        MatGetRow(m_matrix,r,&ncols,&cols,&vals);
        for(j=0;j<ncols;j++){
            MatSetValues(temp,1,&r,1,&cols[j],&aa,ADD_VALUES);// getrow and setvalues can not be used at the same time
                                                              // for the same Mat !!!
        }
        MatRestoreRow(m_matrix,r,&ncols,&cols,&vals);
    }
    MatAssemblyBegin(temp,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(temp,MAT_FINAL_ASSEMBLY);
    MatCopy(temp,m_matrix,SAME_NONZERO_PATTERN);
    MatDestroy(&temp);
    return *this;
}
SparseMatrix& SparseMatrix::operator-=(const SparseMatrix &a){
    if(a.m_m!=m_m || a.m_n!=m_n){
        MessagePrinter::printErrorTxt("can\'t apply -= to sparse matrix with different size");
        MessagePrinter::exitAsFem();
    }
    //Computes Y = a*X + Y.
    MatAXPY(m_matrix,-1.0,a.m_matrix,SAME_NONZERO_PATTERN);
    return *this;
}
//*****************************
//*** * and *= operator
//*****************************
SparseMatrix SparseMatrix::operator*(const double &a)const{
    SparseMatrix anew;
    anew.m_m=m_m;
    anew.m_n=m_n;
    MatDuplicate(m_matrix,MAT_SHARE_NONZERO_PATTERN,&anew.m_matrix);
    MatAXPY(anew.m_matrix,a,m_matrix,SAME_NONZERO_PATTERN);//Y = Y + a â X
    anew.m_allocated=true;
    return anew;
}
SparseMatrix& SparseMatrix::operator*=(const double &a){
    MatScale(m_matrix,a);
    return *this;
}
//*****************************
//*** / and /= operator
//*****************************
SparseMatrix SparseMatrix::operator/(const double &a)const{
    if(abs(a)<1.0e-15){
        MessagePrinter::printErrorTxt("can\'t apply A/0 to sparse matrix, the scalar is singular");
        MessagePrinter::exitAsFem();
    }
    SparseMatrix anew;
    anew.m_m=m_m;
    anew.m_n=m_n;
    MatDuplicate(m_matrix,MAT_SHARE_NONZERO_PATTERN,&anew.m_matrix);
    MatAXPY(anew.m_matrix,1.0/a,m_matrix,SAME_NONZERO_PATTERN);//Y = Y + a â X
    anew.m_allocated=true;
    return anew;
}
SparseMatrix& SparseMatrix::operator/=(const double &a){
    if(abs(a)<1.0e-15){
        MessagePrinter::printErrorTxt("can\'t apply A/0 to sparse matrix, the scalar is singular");
        MessagePrinter::exitAsFem();
    }
    MatScale(m_matrix,1.0/a);
    return *this;
}
//**************************************************
void SparseMatrix::printMatrix(const string &txt)const{
    MessagePrinter::printStars();
    if(txt.size()>0){
        MessagePrinter::printNormalTxt(txt);
    }
    MatView(m_matrix,PETSC_VIEWER_STDOUT_WORLD);
    MessagePrinter::printStars();
}