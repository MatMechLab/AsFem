//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.07.12
//+++ Purpose: initialize our solution system, allocate the necessary
//+++          memory for the vectors.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "SolutionSystem/SolutionSystem.h"

void SolutionSystem::InitSolution(const int &ndofs,const int &nelmts,const int &nnodes,const int &ngp){
    _nDofs=ndofs;
    _nElmts=nelmts;
    _nNodes=nnodes;
    _nGPointsPerBulkElmt=ngp;

    VecCreate(PETSC_COMM_WORLD,&_U);
    VecSetSizes(_U,PETSC_DECIDE,ndofs);
    VecSetUp(_U);// must call this, otherwise PETSc will have memory segmentation error!!!

    VecCreate(PETSC_COMM_WORLD,&_Unew);
    VecSetSizes(_Unew,PETSC_DECIDE,ndofs);
    VecSetUp(_Unew);


    VecCreate(PETSC_COMM_WORLD,&_dU);
    VecSetSizes(_dU,PETSC_DECIDE,ndofs);
    VecSetUp(_dU);

    VecSet(_Unew,0.0);
    VecSet(_U,0.0);
    VecSet(_dU,0.0);
    //*******************************
    VecCreate(PETSC_COMM_WORLD,&_V);
    VecSetSizes(_V,PETSC_DECIDE,ndofs);
    VecSetUp(_V);

    VecSet(_V,0.0);
    //*****************************************************
    //*** For projection array and history variable 
    //*****************************************************
    if(!_HasProjNameList){
        _nProjPerNode=0;
        _ProjectionNameList.clear();
        _IsProjection=false;
    }
    else{
        _nProjPerNode=static_cast<int>(_ProjectionNameList.size());
        _IsProjection=true;
    }
    VecCreate(PETSC_COMM_WORLD,&_Proj);
    VecSetSizes(_Proj,PETSC_DECIDE,_nNodes*(1+_nProjPerNode));
    VecSetUp(_Proj);
    VecSet(_Proj,0.0);
    //*****************************************************
    //*** for scalar type projection quantities
    //*****************************************************
    if(!_HasScalarMateProjName){
        _nScalarProjPerNode=0;
        _ScalarMateProjectionNameList.clear();
    }
    else{
        _nScalarProjPerNode=static_cast<int>(_ScalarMateProjectionNameList.size());
    }
    VecCreate(PETSC_COMM_WORLD,&_ProjScalarMate);
    VecSetSizes(_ProjScalarMate,PETSC_DECIDE,_nNodes*(1+_nScalarProjPerNode));
    VecSetUp(_ProjScalarMate);
    VecSet(_ProjScalarMate,0.0);
    //*****************************************************
    //*** for vector type projection quantities
    //*****************************************************
    if(!_HasVectorMateProjName){
        _nVectorProjPerNode=0;
        _VectorMateProjctionNameList.clear();
    }
    else{
        _nVectorProjPerNode=static_cast<int>(_VectorMateProjctionNameList.size());
    }
    VecCreate(PETSC_COMM_WORLD,&_ProjVectorMate);
    VecSetSizes(_ProjVectorMate,PETSC_DECIDE,_nNodes*(1+_nVectorProjPerNode*3));
    VecSetUp(_ProjVectorMate);
    VecSet(_ProjVectorMate,0.0);
    //*****************************************************
    //*** for rank-2 tensor type projection quantities
    //*****************************************************
    if(!_HasRank2MateProjName){
        _nRank2ProjPerNode=0;
        _Rank2MateProjectionNameList.clear();
    }
    else{
        _nRank2ProjPerNode=static_cast<int>(_Rank2MateProjectionNameList.size());
    }
    VecCreate(PETSC_COMM_WORLD,&_ProjRank2Mate);
    VecSetSizes(_ProjRank2Mate,PETSC_DECIDE,_nNodes*(1+_nRank2ProjPerNode*9));
    VecSetUp(_ProjRank2Mate);
    VecSet(_ProjRank2Mate,0.0);
    //*****************************************************
    //*** for rank-4 tensor type projection quantities
    //*****************************************************
    if(!_HasRank4MateProjName){
        _nRank4ProjPerNode=0;
        _Rank4MateProjectionNameList.clear();
    }
    else{
        _nRank4ProjPerNode=static_cast<int>(_Rank4MateProjectionNameList.size());
    }
    VecCreate(PETSC_COMM_WORLD,&_ProjRank4Mate);
    VecSetSizes(_ProjRank4Mate,PETSC_DECIDE,_nNodes*(1+_nRank4ProjPerNode*36));
    VecSetUp(_ProjRank4Mate);
    VecSet(_ProjRank4Mate,0.0);


    VecCreate(PETSC_COMM_WORLD,&_Hist);
    VecSetSizes(_Hist,PETSC_DECIDE,_nElmts*(_nGPointsPerBulkElmt*_nHistPerGPoint));
    VecSetUp(_Hist);
    VecSet(_Hist,0.0);

    VecDuplicate(_Hist,&_HistOld);

    VecSet(_HistOld,0.0);

    
    _IsInit=true;
}