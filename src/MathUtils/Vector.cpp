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
//+++ Date   : 2022.06.04
//+++ Purpose: implement the API for PETSc and other package's vector
//+++          for dense vector or small size vector, please use
//+++          VectorXd, this is mainly used for the system solution
//+++          and system residual
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MathUtils/Vector.h"

Vector::Vector(){
    m_allocated=false;
    m_size=0;
    m_ghostallocated=false;
}
Vector::Vector(const int &n){
    m_size=n;
    VecCreate(PETSC_COMM_WORLD,&m_vector);
    VecSetSizes(m_vector,PETSC_DECIDE,m_size);
    VecSetFromOptions(m_vector);
    VecSet(m_vector,0.0);
    assemble();
    m_allocated=true;
    m_ghostallocated=false;
}
Vector::Vector(const int &n,const double &val){
    m_size=n;
    VecCreate(PETSC_COMM_WORLD,&m_vector);
    VecSetSizes(m_vector,PETSC_DECIDE,m_size);
    VecSetFromOptions(m_vector);
    VecSet(m_vector,val);
    assemble();
    m_allocated=true;
    m_ghostallocated=false;
}
Vector::Vector(const Vector &a){
    VecDuplicate(a.getVectorCopy(),&m_vector);
    VecCopy(a.m_vector,m_vector);
    assemble();
    m_size=a.getSize();
    m_allocated=true;
    m_ghostallocated=false;
}
//**************************************************
void Vector::setup(){
    if(m_allocated){
        VecDestroy(&m_vector);
    }
    VecCreate(PETSC_COMM_WORLD,&m_vector);
    VecSetSizes(m_vector,PETSC_DECIDE,m_size);
    VecSetFromOptions(m_vector);
    VecSet(m_vector,0.0);
    assemble();
    m_allocated=true;
    m_ghostallocated=false;
}
void Vector::resize(const int &n){
    m_size=n;
    if(m_allocated) VecDestroy(&m_vector);
    
    VecCreate(PETSC_COMM_WORLD,&m_vector);
    VecSetSizes(m_vector,PETSC_DECIDE,m_size);
    VecSetFromOptions(m_vector);
    VecSet(m_vector,0.0);
    assemble();
    m_allocated=true;
    m_ghostallocated=false;
}
void Vector::resize(const int &n,const double &val){
    m_size=n;
    if(m_allocated) VecDestroy(&m_vector);
    
    VecCreate(PETSC_COMM_WORLD,&m_vector);
    VecSetSizes(m_vector,PETSC_DECIDE,m_size);
    VecSetFromOptions(m_vector);
    VecSet(m_vector,val);
    assemble();
    m_allocated=true;
    
    m_ghostallocated=false;
}
//**************************************************
void Vector::makeGhostCopy(){
    if(m_ghostallocated){
        VecDestroy(&m_vector_ghost);
        VecScatterDestroy(&m_scatter);
        m_ghostallocated=false;
    }
    if(m_size){
        VecScatterCreateToAll(m_vector,&m_scatter,&m_vector_ghost);
        VecScatterBegin(m_scatter,m_vector,m_vector_ghost,INSERT_VALUES,SCATTER_FORWARD);
        VecScatterEnd(m_scatter,m_vector,m_vector_ghost,INSERT_VALUES,SCATTER_FORWARD);
        m_ghostallocated=true;
    }
    else{
        m_ghostallocated=false;
    }
}
void Vector::destroyGhostCopy(){
    if(m_ghostallocated){
        VecDestroy(&m_vector_ghost);
        VecScatterDestroy(&m_scatter);
        m_ghostallocated=false;
    }
}
//**************************************************
void Vector::printVec(const string &txt)const{
    if(txt.size()>0){
        MessagePrinter::printNormalTxt(txt);
    }
    Vec seqvec;
    VecScatter scatter;
    VecScatterCreateToAll(m_vector,&scatter,&seqvec);
    VecScatterBegin(scatter,m_vector,seqvec,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(scatter,m_vector,seqvec,INSERT_VALUES,SCATTER_FORWARD);
    double val;
    for(int i=0;i<m_size;i++){
        VecGetValues(seqvec,1,&i,&val);
        PetscPrintf(PETSC_COMM_WORLD,"*** i=%8d, value=%14.6e\n",i+1,val);
    }
    MessagePrinter::printStars();
    VecDestroy(&seqvec);
    VecScatterDestroy(&scatter);
}
