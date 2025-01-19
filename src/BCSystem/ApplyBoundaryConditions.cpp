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
//+++ Date   : 2020.12.26
//+++ Purpose: here we apply the boundary conditions we defined in
//+++          each [bcs] sub blocks
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/BCSystem.h"

void BCSystem::applyBoundaryConditions(const FECalcType &CalcType,
                                       const double &t,const double (&Ctan)[3],
                                       const FECell &t_FECell,
                                       const DofHandler &t_DofHandler,
                                       FE &t_FE,
                                       Vector &U,
                                       Vector &Ucopy,
                                       Vector &Uold,
                                       Vector &Uolder,
                                       Vector &V,
                                       SparseMatrix &AMATRIX,
                                       Vector &RHS){

    double bcvalue;
    for(const auto &it:m_BCBlockList){
        bcvalue=it.m_BCValue;
        m_LocalElmtInfo.m_T=t;
        if(it.m_IsTimeDependent) bcvalue=t*it.m_BCValue;
        if(it.m_BCType==BCType::DIRICHLETBC||
           it.m_BCType==BCType::ROTATEDDIRICHLETBC||
           it.m_BCType==BCType::CYCLICDIRICHLETBC||
           it.m_BCType==BCType::USER1DIRICHLETBC||
           it.m_BCType==BCType::USER2DIRICHLETBC||
           it.m_BCType==BCType::USER3DIRICHLETBC||
           it.m_BCType==BCType::USER4DIRICHLETBC||
           it.m_BCType==BCType::USER5DIRICHLETBC||
           it.m_BCType==BCType::POISSON2DBENCHMARKBC){
            applyDirichletBC(CalcType,bcvalue,it.m_BCType,it.m_JsonParams,it.m_DofIDs,it.m_BoundaryNameList,t_FECell,t_DofHandler,U,Ucopy,Uold,Uolder,V,AMATRIX,RHS);
        }
        else if(it.m_BCType==BCType::NODALDIRICHLETBC){
            continue;
        }
        else if(it.m_BCType==BCType::NODALNEUMANNBC){
            continue;
        }
        else if(it.m_BCType==BCType::NODALFLUXBC){
            continue;
        }
        else if(it.m_BCType==BCType::NODALFORCEBC){
            continue;
        }
        else if(it.m_BCType==BCType::NULLBC){
            continue;
        }
        else{
            applyIntegratedBC(CalcType,bcvalue,it.m_BCType,Ctan,
                              it.m_JsonParams,
                              it.m_DofIDs,
                              it.m_BoundaryNameList,
                              t_FECell,t_DofHandler,t_FE,
                              U,Uold,Uolder,V,
                              AMATRIX,RHS);
        }// end-of-boundary-type-choose
        
    }// end-of-boundary-block-loop

}
//****************************************************************
//*** for preset dirichelt type boundary condition
//****************************************************************
void BCSystem::applyPresetBoundaryConditions(const FECalcType &CalcType,
                                             const double &t,
                                             const FECell &t_FECell,
                                             const DofHandler &t_DofHandler,
                                             Vector &U,
                                             Vector &Ucopy,
                                             Vector &Uold,
                                             Vector &Uolder,
                                             Vector &V,
                                             SparseMatrix &AMATRIX,
                                             Vector &RHS){
    double bcvalue;
    for(const auto &it:m_BCBlockList){
        bcvalue=it.m_BCValue;
        m_LocalElmtInfo.m_T=t;
        if(it.m_IsTimeDependent) bcvalue=t*it.m_BCValue;
        if(it.m_BCType==BCType::DIRICHLETBC||
           it.m_BCType==BCType::ROTATEDDIRICHLETBC||
           it.m_BCType==BCType::CYCLICDIRICHLETBC||
           it.m_BCType==BCType::USER1DIRICHLETBC||
           it.m_BCType==BCType::USER2DIRICHLETBC||
           it.m_BCType==BCType::USER3DIRICHLETBC||
           it.m_BCType==BCType::USER4DIRICHLETBC||
           it.m_BCType==BCType::USER5DIRICHLETBC||
           it.m_BCType==BCType::POISSON2DBENCHMARKBC){
            applyDirichletBC(CalcType,bcvalue,it.m_BCType,it.m_JsonParams,it.m_DofIDs,it.m_BoundaryNameList,t_FECell,t_DofHandler,U,Ucopy,Uold,Uolder,V,AMATRIX,RHS);
        }
        else if(it.m_BCType==BCType::NODALDIRICHLETBC){
            continue;
        }
        else if(it.m_BCType==BCType::NODALNEUMANNBC){
            continue;
        }
        else if(it.m_BCType==BCType::NODALFLUXBC){
            continue;
        }
        else if(it.m_BCType==BCType::NODALFORCEBC){
            continue;
        }
        else if(it.m_BCType==BCType::NULLBC){
            continue;
        }
        else{
            continue;
        }
    }

}