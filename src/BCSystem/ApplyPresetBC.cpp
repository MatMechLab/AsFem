//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2022.01.16
//+++ Purpose: here we apply the dirichlet type boundary conditions
//+++          to 'disp' solution, then, one can use the correct 'disp'
//+++          in the Residual and Jacobian calculation.
//+++          In short, this function should be called before one 
//+++          goes into the FormFE function!
//+++          The old one 'ApplyInitialBC' is deprecated !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/BCSystem.h"
#include "DofHandler/DofHandler.h"

void BCSystem::ApplyPresetBC(const Mesh &mesh,const DofHandler &dofHandler,const FECalcType &calctype,const double &t,const double (&ctan)[3],Vec &U,Mat &AMATRIX,Vec &RHS){

    double bcvalue;
    vector<string> bcnamelist;
    vector<int> DofsIndex;
    if(ctan[0]){}

    _elmtinfo.t=t;
    _elmtinfo.dt=0.0;
    for(auto it:_BCBlockList){
        bcvalue=it._BCValue;
        if(it._IsTimeDependent) bcvalue=it._BCValue*t;
        DofsIndex=it._DofIDs;
        bcnamelist=it._BoundaryNameList;
        if(it._BCType==BCType::DIRICHLETBC||
           it._BCType==BCType::CYCLICDIRICHLETBC||
           it._BCType==BCType::USER1DIRICHLETBC||
           it._BCType==BCType::USER2DIRICHLETBC||
           it._BCType==BCType::USER3DIRICHLETBC||
           it._BCType==BCType::USER4DIRICHLETBC||
           it._BCType==BCType::USER5DIRICHLETBC){
            ApplyDirichletBC(calctype,it._BCType,bcnamelist,DofsIndex,bcvalue,it._Parameters,mesh,dofHandler,U,AMATRIX,RHS);
        }
        else if(it._BCType==BCType::NODALDIRICHLETBC){
            ApplyNodalDirichletBC(calctype,it._BCType,bcnamelist,DofsIndex,bcvalue,it._Parameters,mesh,dofHandler,U,AMATRIX,RHS);
        }
        else{
            continue;
        }//===> end-of-boundary-type-if-else-condition
    }//===> end-of-bcblock-loop

    VecAssemblyBegin(U);
    VecAssemblyEnd(U);
    VecAssemblyBegin(RHS);
    VecAssemblyEnd(RHS);

    MatAssemblyBegin(AMATRIX,MAT_FINAL_ASSEMBLY);

}
