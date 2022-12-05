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
//+++ Date   : 2020.12.26
//+++ Purpose: here we apply the boundary conditions we defined in
//+++          each [bcs] sub blocks
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/BCSystem.h"

void BCSystem::applyBoundaryConditions(const FECalcType &t_calctype,
                                       const double &t,const double (&ctan)[3],
                                       const Mesh &t_mesh,
                                       const DofHandler &t_dofhandler,
                                       FE &t_fe,
                                       Vector &U,Vector &Ucopy,Vector &Uold,Vector &Uolder,
                                       Vector &V,
                                       SparseMatrix &AMATRIX,
                                       Vector &RHS){

    double bcvalue;
    for(const auto &it:m_bclock_list){
        bcvalue=it.m_bcValue;
        m_local_elmtinfo.m_t=t;
        if(it.m_isTimeDependent) bcvalue=t*it.m_bcValue;
        if(it.m_bcType==BCType::DIRICHLETBC||
           it.m_bcType==BCType::ROTATEDDIRICHLETBC||
           it.m_bcType==BCType::CYCLICDIRICHLETBC||
           it.m_bcType==BCType::USER1DIRICHLETBC||
           it.m_bcType==BCType::USER2DIRICHLETBC||
           it.m_bcType==BCType::USER3DIRICHLETBC||
           it.m_bcType==BCType::USER4DIRICHLETBC||
           it.m_bcType==BCType::USER5DIRICHLETBC||
           it.m_bcType==BCType::POISSON2DBENCHMARKBC){
            applyDirichletBC(t_calctype,bcvalue,it.m_bcType,it.m_json_params,it.m_dofIDs,it.m_boundaryNameList,t_mesh,t_dofhandler,U,Ucopy,Uold,Uolder,V,AMATRIX,RHS);
        }
        else if(it.m_bcType==BCType::NODALDIRICHLETBC){
            continue;
        }
        else if(it.m_bcType==BCType::NODALNEUMANNBC){
            continue;
        }
        else if(it.m_bcType==BCType::NODALFLUXBC){
            continue;
        }
        else if(it.m_bcType==BCType::NODALFORCEBC){
            continue;
        }
        else if(it.m_bcType==BCType::NULLBC){
            continue;
        }
        else{
            applyIntegratedBC(t_calctype,bcvalue,it.m_bcType,ctan,
                              it.m_json_params,
                              it.m_dofIDs,
                              it.m_boundaryNameList,
                              t_mesh,t_dofhandler,t_fe,
                              U,Uold,Uolder,V,
                              AMATRIX,RHS);
        }// end-of-boundary-type-choose
        
    }// end-of-boundary-block-loop

}
//****************************************************************
//*** for preset dirichelt type boundary condition
//****************************************************************
void BCSystem::applyPresetBoundaryConditions(const FECalcType &t_calctype,
                                 const double &t,
                                 const Mesh &t_mesh,
                                 const DofHandler &t_dofhandler,
                                 Vector &U,Vector &Ucopy,Vector &Uold,Vector &Uolder,
                                 Vector &V,
                                 SparseMatrix &AMATRIX,
                                 Vector &RHS){
    double bcvalue;
    for(const auto &it:m_bclock_list){
        bcvalue=it.m_bcValue;
        m_local_elmtinfo.m_t=t;
        if(it.m_isTimeDependent) bcvalue=t*it.m_bcValue;
        if(it.m_bcType==BCType::DIRICHLETBC||
           it.m_bcType==BCType::ROTATEDDIRICHLETBC||
           it.m_bcType==BCType::CYCLICDIRICHLETBC||
           it.m_bcType==BCType::USER1DIRICHLETBC||
           it.m_bcType==BCType::USER2DIRICHLETBC||
           it.m_bcType==BCType::USER3DIRICHLETBC||
           it.m_bcType==BCType::USER4DIRICHLETBC||
           it.m_bcType==BCType::USER5DIRICHLETBC||
           it.m_bcType==BCType::POISSON2DBENCHMARKBC){
            applyDirichletBC(t_calctype,bcvalue,it.m_bcType,it.m_json_params,it.m_dofIDs,it.m_boundaryNameList,t_mesh,t_dofhandler,U,Ucopy,Uold,Uolder,V,AMATRIX,RHS);
        }
        else if(it.m_bcType==BCType::NODALDIRICHLETBC){
            continue;
        }
        else if(it.m_bcType==BCType::NODALNEUMANNBC){
            continue;
        }
        else if(it.m_bcType==BCType::NODALFLUXBC){
            continue;
        }
        else if(it.m_bcType==BCType::NODALFORCEBC){
            continue;
        }
        else if(it.m_bcType==BCType::NULLBC){
            continue;
        }
        else{
            continue;
        }
    }

}