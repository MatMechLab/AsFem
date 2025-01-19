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
bool SNESSolver::solve(FECell &t_FECell,
                       DofHandler &t_DofHandler,
                       FE &t_FE,
                       ElmtSystem &t_ElmtSystem,
                       MateSystem &t_MateSystem,
                       FESystem &t_FESystem,
                       BCSystem &t_BCSystem,
                       SolutionSystem &t_SolnSystem,
                       EquationSystem &t_EqSystem,
                       LinearSolver &t_LinearSolver,
                       FEControlInfo &t_FECtrlInfo){

    if (t_LinearSolver.getIterationNumber()){}
    t_SolnSystem.m_Ucopy.copyFrom(t_SolnSystem.m_Ucurrent);
    
    m_AppCtx=AppCtx{&t_FECell,&t_DofHandler,
                   &t_BCSystem,
                   &t_ElmtSystem,&t_MateSystem,
                   &t_SolnSystem,&t_EqSystem,
                   &t_FE,&t_FESystem,
                   &t_FECtrlInfo
                   };
    m_MonCtx=MonitorCtx{0.0,1.0,
                0.0,1.0,
                0.0,1.0,
                0,
                t_FECtrlInfo.IsDepDebug};


    m_AppCtx._bcSystem->applyPresetBoundaryConditions(FECalcType::UPDATEU,
                                             m_AppCtx._fectrlinfo->T+m_AppCtx._fectrlinfo->Dt,
                                             *m_AppCtx._fecell,
                                             *m_AppCtx._dofHandler,
                                             m_AppCtx._solutionSystem->m_Ucurrent,
                                             m_AppCtx._solutionSystem->m_Ucopy,
                                             m_AppCtx._solutionSystem->m_Uold,
                                             m_AppCtx._solutionSystem->m_Uolder,
                                             m_AppCtx._solutionSystem->m_V,
                                             m_AppCtx._equationSystem->m_AMATRIX,
                                             m_AppCtx._equationSystem->m_RHS);

    SNESSetFunction(m_SNES,m_AppCtx._equationSystem->m_RHS.getVectorCopy(),computeResidual,&m_AppCtx);
    SNESSetJacobian(m_SNES,m_AppCtx._equationSystem->m_AMATRIX.getCopy(),m_AppCtx._equationSystem->m_AMATRIX.getCopy(),computeJacobian,&m_AppCtx);
    if (m_NLSolverType==NonlinearSolverType::NEWTONAL) {
        SNESNewtonALSetFunction(m_SNES,computeResidual,NULL);
    }
    SNESMonitorSet(m_SNES,myMonitor,&m_MonCtx,0);
    SNESSetForceIteration(m_SNES,PETSC_TRUE);
    SNESSetFromOptions(m_SNES);
    SNESSolve(m_SNES,NULL,m_AppCtx._solutionSystem->m_Ucurrent.getVectorRef());
    SNESGetConvergedReason(m_SNES,&m_SNESConvergeReason);

    m_Iterations=m_MonCtx.iters;
    m_RNorm=m_MonCtx.rnorm;
    m_AbsTolDu=m_MonCtx.dunorm;
    m_AbsTolE=m_MonCtx.enorm;

    char buff[68];//77-12=65
    string str;

    if(m_SNESConvergeReason==SNES_CONVERGED_FNORM_ABS){
        if(t_FECtrlInfo.IsDepDebug){
            snprintf(buff,68,"  Converged for |R|<atol, final iters=%3d",m_MonCtx.iters);
            str=buff;
            MessagePrinter::printNormalTxt(str);
        }
        else{
            snprintf(buff,68,"  SNES solver: iters=%3d,|R0|=%12.5e,|R|=%12.5e",m_Iterations,m_MonCtx.rnorm0,m_RNorm);
            str=buff;
            MessagePrinter::printNormalTxt(str);
        }
        return true;
    }
    else if(m_SNESConvergeReason==SNES_CONVERGED_FNORM_RELATIVE){
        if(t_FECtrlInfo.IsDepDebug){
            snprintf(buff,68,"  Converged for |R|<rtol*|R0|, final iters=%3d",m_MonCtx.iters);
            str=buff;
            MessagePrinter::printNormalTxt(str);
        }
        else{
            snprintf(buff,68,"  SNES solver: iters=%3d,|R0|=%12.5e,|R|=%12.5e",m_Iterations,m_MonCtx.rnorm0,m_RNorm);
            str=buff;
            MessagePrinter::printNormalTxt(str);
        }
        return true;
    }
    else if(m_SNESConvergeReason==SNES_CONVERGED_SNORM_RELATIVE){
        if(t_FECtrlInfo.IsDepDebug){
            snprintf(buff,68,"  Converged for |delta x|<stol|x|, final iters=%3d",m_MonCtx.iters);
            str=buff;
            MessagePrinter::printNormalTxt(str);
        }
        else{
            snprintf(buff,68,"  SNES solver: iters=%3d,|R0|=%12.5e,|R|=%12.5e",m_Iterations,m_MonCtx.rnorm0,m_RNorm);
            str=buff;
            MessagePrinter::printNormalTxt(str);
        }
        return true;
    }
    else{
        snprintf(buff,68,"  Divergent, SNES nonlinear solver failed, iters=%3d",m_MonCtx.iters);
        str=buff;
        MessagePrinter::printNormalTxt(str);
        return false;
    }
    return false;
}
