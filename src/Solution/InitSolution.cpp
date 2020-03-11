//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "Solution/Solution.h"

void Solution::InitSolution(const PetscInt &ndofs,const PetscInt &nelmts,const PetscInt &nnodes,const PetscInt &ngp){
    
    _nDofs=ndofs;
    _nElmts=nelmts;
    _nNodes=nnodes;
    _nGPointsPerBulkElmt=ngp;

    VecCreate(PETSC_COMM_WORLD,&_Uold);
    VecSetSizes(_Uold,PETSC_DECIDE,ndofs);
    VecSetUp(_Uold);// must call this, otherwise PETSc will have memory segmentation error!!!
   

    VecCreate(PETSC_COMM_WORLD,&_Uolder);
    VecSetSizes(_Uolder,PETSC_DECIDE,ndofs);
    VecSetUp(_Uolder);

    VecCreate(PETSC_COMM_WORLD,&_Unew);
    VecSetSizes(_Unew,PETSC_DECIDE,ndofs);
    VecSetUp(_Unew);

    VecCreate(PETSC_COMM_WORLD,&_Utemp);
    VecSetSizes(_Utemp,PETSC_DECIDE,ndofs);
    VecSetUp(_Utemp);

    VecCreate(PETSC_COMM_WORLD,&_dU);
    VecSetSizes(_dU,PETSC_DECIDE,ndofs);
    VecSetUp(_dU);

    VecSet(_Uold,0.0);
    VecSet(_Uolder,0.0);
    VecSet(_Unew,0.0);
    VecSet(_Utemp,0.0);
    VecSet(_dU,0.0);
    //*******************************
    VecCreate(PETSC_COMM_WORLD,&_V);
    VecSetSizes(_V,PETSC_DECIDE,ndofs);
    VecSetUp(_V);

    VecCreate(PETSC_COMM_WORLD,&_Vold);
    VecSetSizes(_Vold,PETSC_DECIDE,ndofs);
    VecSetUp(_Vold);
    
    VecCreate(PETSC_COMM_WORLD,&_Volder);
    VecSetSizes(_Volder,PETSC_DECIDE,ndofs);
    VecSetUp(_Volder);


    VecSet(_V,0.0);
    VecSet(_Vold,0.0);
    VecSet(_Volder,0.0);


    //*****************************************************
    //*** For projection array and history variable 
    //*****************************************************
    if(!_HasProjNameList){
        _nProjPerNode=9;
        _ProjectionNameList.clear();
        _ProjectionNameList.reserve(_nProjPerNode);
        string name;
        for(int i=1;i<=_nProjPerNode;++i){
            name="Proj-"+to_string(i);
            _ProjectionNameList.push_back(name);
        }
    }
    else{
        _nProjPerNode=(PetscInt)_ProjectionNameList.size();
    }
    VecCreate(PETSC_COMM_WORLD,&_Proj);
    VecSetSizes(_Proj,PETSC_DECIDE,_nNodes*(1+_nProjPerNode));
    VecSetUp(_Proj);
    VecSet(_Proj,0.0);


    _nHistPerGPoint=6;

    VecCreate(PETSC_COMM_WORLD,&_Hist);
    VecSetSizes(_Hist,PETSC_DECIDE,_nElmts*(_nGPointsPerBulkElmt*_nHistPerGPoint));
    VecSetUp(_Hist);
    VecSet(_Hist,0.0);

    VecDuplicate(_Hist,&_HistOld);

    VecSet(_HistOld,0.0);

    // VecCreate(PETSC_COMM_WORLD,&_HistOld);
    // VecSetSizes(_HistOld,PETSC_DECIDE,_nElmts*(_nGPointsPerBulkElmt*_nHistPerGPoint));
    // VecSetUp(_HistOld);
    // VecSet(_HistOld,0.0);

    
    _IsInit=true;
}