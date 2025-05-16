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

#include "MathUtils/Vector.h"

Vector::Vector(){
    m_Allocated=false;
    m_Size=0;
    m_GhostAllocated=false;
}
Vector::Vector(const int &n){
    m_Size=n;
    VecCreate(PETSC_COMM_WORLD,&m_Vector);
    VecSetSizes(m_Vector,PETSC_DECIDE,m_Size);
    VecSetFromOptions(m_Vector);
    VecSet(m_Vector,0.0);
    assemble();
    m_Allocated=true;
    m_GhostAllocated=false;
}
Vector::Vector(const int &n,const double &val){
    m_Size=n;
    VecCreate(PETSC_COMM_WORLD,&m_Vector);
    VecSetSizes(m_Vector,PETSC_DECIDE,m_Size);
    VecSetFromOptions(m_Vector);
    VecSet(m_Vector,val);
    assemble();
    m_Allocated=true;
    m_GhostAllocated=false;
}
Vector::Vector(const Vector &a){
    VecDuplicate(a.getVectorCopy(),&m_Vector);
    VecCopy(a.m_Vector,m_Vector);
    assemble();
    m_Size=a.getSize();
    m_Allocated=true;
    m_GhostAllocated=false;
}
//**************************************************
void Vector::setup(){
    if(m_Allocated){
        VecDestroy(&m_Vector);
    }
    VecCreate(PETSC_COMM_WORLD,&m_Vector);
    VecSetSizes(m_Vector,PETSC_DECIDE,m_Size);
    VecSetFromOptions(m_Vector);
    VecSet(m_Vector,0.0);
    assemble();
    m_Allocated=true;
    m_GhostAllocated=false;
}
void Vector::resize(const int &n){
    m_Size=n;
    if(m_Allocated) VecDestroy(&m_Vector);
    
    VecCreate(PETSC_COMM_WORLD,&m_Vector);
    VecSetSizes(m_Vector,PETSC_DECIDE,m_Size);
    VecSetFromOptions(m_Vector);
    VecSet(m_Vector,0.0);
    assemble();
    m_Allocated=true;
    m_GhostAllocated=false;
}
void Vector::resize(const int &n,const double &val){
    m_Size=n;
    if(m_Allocated) VecDestroy(&m_Vector);
    
    VecCreate(PETSC_COMM_WORLD,&m_Vector);
    VecSetSizes(m_Vector,PETSC_DECIDE,m_Size);
    VecSetFromOptions(m_Vector);
    VecSet(m_Vector,val);
    assemble();
    m_Allocated=true;
    
    m_GhostAllocated=false;
}
//**************************************************
void Vector::makeGhostCopy(){
    if(m_GhostAllocated){
        VecDestroy(&m_VectorGhost);
        VecScatterDestroy(&m_Scatter);
        m_GhostAllocated=false;
    }
    if(m_Size){
        VecScatterCreateToAll(m_Vector,&m_Scatter,&m_VectorGhost);
        VecScatterBegin(m_Scatter,m_Vector,m_VectorGhost,INSERT_VALUES,SCATTER_FORWARD);
        VecScatterEnd(m_Scatter,m_Vector,m_VectorGhost,INSERT_VALUES,SCATTER_FORWARD);
        m_GhostAllocated=true;
    }
    else{
        m_GhostAllocated=false;
    }
}
void Vector::destroyGhostCopy(){
    if(m_GhostAllocated){
        VecDestroy(&m_VectorGhost);
        VecScatterDestroy(&m_Scatter);
        m_GhostAllocated=false;
    }
}
//**************************************************
void Vector::printVec(const string &txt)const{
    if(txt.size()>0){
        MessagePrinter::printNormalTxt(txt);
    }
    Vec seqvec;
    VecScatter scatter;
    VecScatterCreateToAll(m_Vector,&scatter,&seqvec);
    VecScatterBegin(scatter,m_Vector,seqvec,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(scatter,m_Vector,seqvec,INSERT_VALUES,SCATTER_FORWARD);
    double val;
    for(int i=0;i<m_Size;i++){
        VecGetValues(seqvec,1,&i,&val);
        PetscPrintf(PETSC_COMM_WORLD,"*** i=%8d, value=%14.6e\n",i+1,val);
    }
    MessagePrinter::printStars();
    VecDestroy(&seqvec);
    VecScatterDestroy(&scatter);
}
