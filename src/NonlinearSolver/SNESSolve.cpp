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
//+++ Purpose: here we call the SNES solver from PETSc to solve
//+++          our nonlinear system equation R(x)->0
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "NonlinearSolver/SNESSolver.h"
#include "TimeStepping/TimeSteppingTool.h"

//***************************************************************
//*** here we define a monitor to print out the iteration info
//***************************************************************
PetscErrorCode myMonitor(SNES snes,PetscInt iters,PetscReal rnorm,void* ctx){
    char buff[68];
    string str;
    MonitorCtx *user=(MonitorCtx*)ctx;
    user->iters=iters;
    user->rnorm=rnorm;
    SNESGetUpdateNorm(snes,&user->dunorm);
    if(iters==0) user->dunorm=1.0;

    user->enorm=rnorm*user->dunorm;
    if(iters==0){
        user->rnorm0=rnorm;
        user->dunorm0=user->dunorm;
        user->enorm0=user->enorm;
    }
    if(user->IsDepDebug){
        snprintf(buff,68,"  SNES solver: iters=%4d, |R|=%12.5e, |dU|=%12.5e",iters,rnorm,user->dunorm);
        str=buff;
        MessagePrinter::printNormalTxt(str);
    }
    return 0;
}
//***************************************************************
//*** here we setup the subroutine for residual 
//***************************************************************
PetscErrorCode computeResidual(SNES snes,Vec U,Vec RHS,void *ctx){
    AppCtx *user=(AppCtx*)ctx;
    int i;
    SNESGetMaxNonlinearStepFailures(snes,&i);// just to get rid of unused snes warning

    computeTimeDerivatives(*user->_fectrlinfo,U,*user->_solutionSystem);// copy U vec to Vector of AsFem

    user->_feSystem->formBulkFE(FECalcType::COMPUTERESIDUAL,
                                user->_fectrlinfo->T+user->_fectrlinfo->Dt,
                                user->_fectrlinfo->Dt,
                                user->_fectrlinfo->Ctan,
                                *user->_fecell,
                                *user->_dofHandler,
                                *user->_fe,
                                *user->_elmtSystem,
                                *user->_mateSystem,
                                *user->_solutionSystem,
                                user->_equationSystem->m_AMATRIX,
                                user->_equationSystem->m_RHS);
    
    user->_bcSystem->setDirichletPenalty(user->_feSystem->getMaxCoefOfKMatrix()*1.0e10);

    user->_solutionSystem->m_Ucopy.copyFrom(U);
    user->_bcSystem->applyBoundaryConditions(FECalcType::COMPUTERESIDUAL,
                                             user->_fectrlinfo->T+user->_fectrlinfo->Dt,
                                             user->_fectrlinfo->Ctan,
                                             *user->_fecell,
                                             *user->_dofHandler,
                                             *user->_fe,
                                             user->_solutionSystem->m_Utemp,
                                             user->_solutionSystem->m_Ucopy,
                                             user->_solutionSystem->m_Uold,
                                             user->_solutionSystem->m_Uolder,
                                             user->_solutionSystem->m_V,
                                             user->_equationSystem->m_AMATRIX,
                                             user->_equationSystem->m_RHS);

    user->_equationSystem->m_RHS.copy2Vec(RHS);// please overwrite RHS to make sure it is the correct one !

    return 0;
}

//***************************************************************
//*** here we setup the subroutine for jacobian 
//***************************************************************
PetscErrorCode computeJacobian(SNES snes,Vec U,Mat Jac,Mat B,void *ctx){
    AppCtx *user=(AppCtx*)ctx;
    
    user->_feSystem->resetMaxKMatrixCoeff();
    
    computeTimeDerivatives(*user->_fectrlinfo,U,*user->_solutionSystem);// compute the Utemp vector, and other derivatives
    
    user->_feSystem->formBulkFE(FECalcType::COMPUTEJACOBIAN,
                                user->_fectrlinfo->T+user->_fectrlinfo->Dt,
                                user->_fectrlinfo->Dt,
                                user->_fectrlinfo->Ctan,
                                *user->_fecell,
                                *user->_dofHandler,
                                *user->_fe,
                                *user->_elmtSystem,
                                *user->_mateSystem,
                                *user->_solutionSystem,
                                user->_equationSystem->m_AMATRIX,
                                user->_equationSystem->m_RHS);
    
    user->_bcSystem->setDirichletPenalty(user->_feSystem->getMaxCoefOfKMatrix()*1.0e10);

    user->_solutionSystem->m_Ucopy.copyFrom(U);
    user->_bcSystem->applyBoundaryConditions(FECalcType::COMPUTEJACOBIAN,
                                             user->_fectrlinfo->T+user->_fectrlinfo->Dt,
                                             user->_fectrlinfo->Ctan,
                                             *user->_fecell,
                                             *user->_dofHandler,
                                             *user->_fe,
                                             user->_solutionSystem->m_Utemp,
                                             user->_solutionSystem->m_Ucopy,
                                             user->_solutionSystem->m_Uold,
                                             user->_solutionSystem->m_Uolder,
                                             user->_solutionSystem->m_V,
                                             user->_equationSystem->m_AMATRIX,
                                             user->_equationSystem->m_RHS);

    int i;
    SNESGetMaxNonlinearStepFailures(snes,&i);

    if(Jac!=B){
        MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY);
    }

    return 0;
}
//***************************************************************
//*** here we solve our nonlinear equations for R(x)->0 
//***************************************************************
bool SNESSolver::solve(FECell &fecell,DofHandler &dofhandler,FE &fe,
                       ElmtSystem &elmtsystem,MateSystem &matesystem,
                       FESystem &fesystem,
                       BCSystem &bcsystem,
                       SolutionSystem &solutionsystem,
                       EquationSystem &equationsystem,
                       FEControlInfo &fectrlinfo){

    solutionsystem.m_Ucopy.copyFrom(solutionsystem.m_Ucurrent);
    
    m_appctx=AppCtx{&fecell,&dofhandler,
                   &bcsystem,
                   &elmtsystem,&matesystem,
                   &solutionsystem,&equationsystem,
                   &fe,&fesystem,
                   &fectrlinfo
                   };
    m_monctx=MonitorCtx{0.0,1.0,
                0.0,1.0,
                0.0,1.0,
                0,
                fectrlinfo.IsDepDebug};


    m_appctx._bcSystem->applyPresetBoundaryConditions(FECalcType::UPDATEU,
                                             m_appctx._fectrlinfo->T+m_appctx._fectrlinfo->Dt,
                                             *m_appctx._fecell,
                                             *m_appctx._dofHandler,
                                             m_appctx._solutionSystem->m_Ucurrent,
                                             m_appctx._solutionSystem->m_Ucopy,
                                             m_appctx._solutionSystem->m_Uold,
                                             m_appctx._solutionSystem->m_Uolder,
                                             m_appctx._solutionSystem->m_V,
                                             m_appctx._equationSystem->m_AMATRIX,
                                             m_appctx._equationSystem->m_RHS);

    SNESSetFunction(m_snes,m_appctx._equationSystem->m_RHS.getVectorCopy(),computeResidual,&m_appctx);
    SNESSetJacobian(m_snes,m_appctx._equationSystem->m_AMATRIX.getCopy(),m_appctx._equationSystem->m_AMATRIX.getCopy(),computeJacobian,&m_appctx);
    if (m_nlsolvertype==NonlinearSolverType::NEWTONAL) {
        SNESNewtonALSetFunction(m_snes,computeResidual,NULL);
    }
    SNESMonitorSet(m_snes,myMonitor,&m_monctx,0);
    SNESSetForceIteration(m_snes,PETSC_TRUE);
    SNESSetFromOptions(m_snes);
    SNESSolve(m_snes,NULL,m_appctx._solutionSystem->m_Ucurrent.getVectorRef());
    SNESGetConvergedReason(m_snes,&m_snesconvergereason);

    m_iterations=m_monctx.iters;
    m_rnorm=m_monctx.rnorm;
    m_abstol_du=m_monctx.dunorm;
    m_abstol_e=m_monctx.enorm;

    char buff[68];//77-12=65
    string str;

    if(m_snesconvergereason==SNES_CONVERGED_FNORM_ABS){
        if(fectrlinfo.IsDepDebug){
            snprintf(buff,68,"  Converged for |R|<atol, final iters=%3d",m_monctx.iters);
            str=buff;
            MessagePrinter::printNormalTxt(str);
        }
        else{
            snprintf(buff,68,"  SNES solver: iters=%3d,|R0|=%12.5e,|R|=%12.5e",m_iterations,m_monctx.rnorm0,m_rnorm);
            str=buff;
            MessagePrinter::printNormalTxt(str);
        }
        return true;
    }
    else if(m_snesconvergereason==SNES_CONVERGED_FNORM_RELATIVE){
        if(fectrlinfo.IsDepDebug){
            snprintf(buff,68,"  Converged for |R|<rtol*|R0|, final iters=%3d",m_monctx.iters);
            str=buff;
            MessagePrinter::printNormalTxt(str);
        }
        else{
            snprintf(buff,68,"  SNES solver: iters=%3d,|R0|=%12.5e,|R|=%12.5e",m_iterations,m_monctx.rnorm0,m_rnorm);
            str=buff;
            MessagePrinter::printNormalTxt(str);
        }
        return true;
    }
    else if(m_snesconvergereason==SNES_CONVERGED_SNORM_RELATIVE){
        if(fectrlinfo.IsDepDebug){
            snprintf(buff,68,"  Converged for |delta x|<stol|x|, final iters=%3d",m_monctx.iters);
            str=buff;
            MessagePrinter::printNormalTxt(str);
        }
        else{
            snprintf(buff,68,"  SNES solver: iters=%3d,|R0|=%12.5e,|R|=%12.5e",m_iterations,m_monctx.rnorm0,m_rnorm);
            str=buff;
            MessagePrinter::printNormalTxt(str);
        }
        return true;
    }
    else{
        snprintf(buff,68,"  Divergent, SNES nonlinear solver failed, iters=%3d",m_monctx.iters);
        str=buff;
        MessagePrinter::printNormalTxt(str);
        return false;
    }
    return false;
}
