//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "FE/FE.h"

void FE::PrintQPointInfo() const{
    PetscPrintf(PETSC_COMM_WORLD,"*** Gauss integration point information:                              ***\n");
    if(_qp_bulk.GetQpType()=="gauss"){
        PetscPrintf(PETSC_COMM_WORLD,"***   bulk qpoint order=%2d, qpoint number=%4d, type=gauss            ***\n",_qp_bulk.GetQpOrder(),_qp_bulk.GetQpPointsNum());
        if(_qp_bulk.GetDim()==2){
            PetscPrintf(PETSC_COMM_WORLD,"***   line qpoint order=%2d, qpoint number=%4d, type=gauss            ***\n",_qp_line.GetQpOrder(),_qp_line.GetQpPointsNum());
        }
        if(_qp_bulk.GetDim()==3){
            PetscPrintf(PETSC_COMM_WORLD,"***   surf qpoint order=%2d, qpoint number=%4d, type=gauss            ***\n",_qp_surface.GetQpOrder(),_qp_surface.GetQpPointsNum());
        }
    }
    else if(_qp_bulk.GetQpType()=="gausslobatto"){
        PetscPrintf(PETSC_COMM_WORLD,"***   bulk qpoint order=%2d, qpoint number=%4d, type=gausslobatto    ***\n",_qp_bulk.GetQpOrder(),_qp_bulk.GetQpPointsNum());
        if(_qp_bulk.GetDim()==2){
            PetscPrintf(PETSC_COMM_WORLD,"***   line qpoint order=%2d, qpoint number=%4d, type=gausslobatto     ***\n",_qp_line.GetQpOrder(),_qp_line.GetQpPointsNum());
        }
        if(_qp_bulk.GetDim()==3){
            PetscPrintf(PETSC_COMM_WORLD,"***   surf qpoint order=%2d, qpoint number=%4d, type=gausslobatto     ***\n",_qp_surface.GetQpOrder(),_qp_surface.GetQpPointsNum());
        }
    }
    PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n");
       
}