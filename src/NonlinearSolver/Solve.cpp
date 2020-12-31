//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.12.26
//+++ Purpose: here we call the SNES solver from PETSc to solve
//+++          our nonlinear system equation R(x)->0
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "NonlinearSolver/NonlinearSolver.h"


//***************************************************************
//*** here we define a monitor to print out the iteration info
//***************************************************************
PetscErrorCode Monitor(SNES snes,PetscInt iters,PetscReal rnorm,void* ctx){
    char buff[68];
    string str;
    MonitorCtx *user=(MonitorCtx*)ctx;
    user->iters=iters;
    user->rnorm=rnorm;
    if(iters==0){
        SNESGetSolutionNorm(snes,&user->dunorm);
    }
    else{
        SNESGetUpdateNorm(snes,&user->dunorm);
    }
    user->enorm=rnorm*user->dunorm;
    if(iters==0){
        user->rnorm0=rnorm;
        user->dunorm0=user->dunorm;
        user->enorm0=user->enorm;
    }
    if(user->IsDepDebug){
        snprintf(buff,68,"  SNES solver: iters=%3d,|R|=%12.5e,|dU|=%12.5e",iters,rnorm,user->dunorm);
        str=buff;
        MessagePrinter::PrintNormalTxt(str);
    }
    return 0;
}

//***************************************************************
//*** here we setup the subroutine for residual 
//***************************************************************
PetscErrorCode ComputeResidual(SNES snes,Vec U,Vec RHS,void *ctx){
    AppCtx *user=(AppCtx*)ctx;
    int i;
    SNESGetMaxNonlinearStepFailures(snes,&i);// just to get rid of unused snes warning
    user->_bcSystem.ApplyInitialBC(user->_mesh,user->_dofHandler,user->_fectrlinfo.t,U);

    // calculate the current velocity
    VecWAXPY(user->_solutionSystem._V,-1.0,user->_solutionSystem._Uold,U);//V=-Uold+Unew
    VecScale(user->_solutionSystem._V,user->_fectrlinfo.ctan[1]);//V=V*1.0/dt
    
    user->_feSystem.FormBulkFE(FECalcType::ComputeResidual,
                        user->_fectrlinfo.t,user->_fectrlinfo.dt,user->_fectrlinfo.ctan,
                        user->_mesh,user->_dofHandler,user->_fe,
                        user->_elmtSystem,user->_mateSystem,
                        U,user->_solutionSystem._V,
                        user->_solutionSystem._Hist,user->_solutionSystem._HistOld,
                        user->_solutionSystem._Proj,
                        user->_equationSystem._AMATRIX,RHS);
    
    user->_bcSystem.SetBCPenaltyFactor(user->_feSystem.GetMaxAMatrixValue()*1.0e8);

    user->_bcSystem.ApplyBC(user->_mesh,user->_dofHandler,user->_fe,
                    FECalcType::ComputeResidual,user->_fectrlinfo.t,user->_fectrlinfo.ctan,U,
                    user->_equationSystem._AMATRIX,RHS);
    
    return 0;
}

//***************************************************************
//*** here we setup the subroutine for jacobian 
//***************************************************************
PetscErrorCode ComputeJacobian(SNES snes,Vec U,Mat A,Mat B,void *ctx){
    AppCtx *user=(AppCtx*)ctx;
    int i;

    user->_feSystem.ResetMaxAMatrixValue();
    user->_bcSystem.ApplyInitialBC(user->_mesh,user->_dofHandler,user->_fectrlinfo.t,U);

    //*** calculate the current velocity
    VecWAXPY(user->_solutionSystem._V,-1.0,user->_solutionSystem._Uold,U);//V=-Uold+Unew
    VecScale(user->_solutionSystem._V,user->_fectrlinfo.ctan[1]);//V=V*1.0/dt

    user->_feSystem.FormBulkFE(FECalcType::ComputeJacobian,
                        user->_fectrlinfo.t,user->_fectrlinfo.dt,user->_fectrlinfo.ctan,
                        user->_mesh,user->_dofHandler,user->_fe,
                        user->_elmtSystem,user->_mateSystem,
                        U,user->_solutionSystem._V,
                        user->_solutionSystem._Hist,user->_solutionSystem._HistOld,
                        user->_solutionSystem._Proj,
                        A,user->_equationSystem._RHS);
    
    if(user->_feSystem.GetMaxAMatrixValue()>1.0e12){
        user->_bcSystem.SetBCPenaltyFactor(1.0e20);
    }
    else if(user->_feSystem.GetMaxAMatrixValue()>1.0e6&&user->_feSystem.GetMaxAMatrixValue()<=1.0e12){
        user->_bcSystem.SetBCPenaltyFactor(user->_feSystem.GetMaxAMatrixValue()*1.0e8);
    }
    else if(user->_feSystem.GetMaxAMatrixValue()>1.0e3&&user->_feSystem.GetMaxAMatrixValue()<=1.0e6){
        user->_bcSystem.SetBCPenaltyFactor(user->_feSystem.GetMaxAMatrixValue()*1.0e12);
    }
    else{
        user->_bcSystem.SetBCPenaltyFactor(user->_feSystem.GetMaxAMatrixValue()*1.0e16);
    }

    user->_bcSystem.ApplyBC(user->_mesh,user->_dofHandler,user->_fe,
                    FECalcType::ComputeJacobian,user->_fectrlinfo.t,user->_fectrlinfo.ctan,U,
                    A,user->_equationSystem._RHS);

    MatGetSize(B,&i,&i);
    SNESGetMaxNonlinearStepFailures(snes,&i);

    return 0;
}

//***************************************************************
//*** here we solve our nonlinear equations for R(x)->0 
//***************************************************************
bool NonlinearSolver::Solve(Mesh &mesh,DofHandler &dofHandler,
                        ElmtSystem &elmtSystem,MateSystem &mateSystem,
                        BCSystem &bcSystem,ICSystem &icSystem,
                        SolutionSystem &solutionSystem,EquationSystem &equationSystem,
                        FE &fe,FESystem &feSystem,
                        FEControlInfo &fectrlinfo){
    
    _appctx=AppCtx{mesh,dofHandler,
                   bcSystem,icSystem,
                   elmtSystem,mateSystem,
                   solutionSystem,equationSystem,
                   fe,feSystem,
                   fectrlinfo
                   };
    
    _monctx=MonitorCtx{0.0,1.0,
            0.0,1.0,
            0.0,1.0,
            0,
            fectrlinfo.IsDepDebug};


    _appctx._bcSystem.ApplyInitialBC(_appctx._mesh,_appctx._dofHandler,1.0,_appctx._solutionSystem._Unew);
    
    SNESSetFunction(_snes,_appctx._equationSystem._RHS,ComputeResidual,&_appctx);

    SNESSetJacobian(_snes,_appctx._equationSystem._AMATRIX,_appctx._equationSystem._AMATRIX,ComputeJacobian,&_appctx);

    SNESMonitorSet(_snes,Monitor,&_monctx,0);

    SNESSetForceIteration(_snes,PETSC_TRUE);

    SNESSetFromOptions(_snes);

    
    SNESSolve(_snes,NULL,_appctx._solutionSystem._Unew);
    

    SNESGetConvergedReason(_snes,&_snesreason);
    
    _Iters=_monctx.iters;
    _Rnorm=_monctx.rnorm;
    _dUnorm=_monctx.dunorm;
    _Enorm=_monctx.enorm;

    char buff[65];//77-12=65
    char buffnew[68];
    string str;

    if(_snesreason==SNES_CONVERGED_FNORM_ABS){
        if(fectrlinfo.IsDepDebug){
            snprintf(buff,65,"  Converged for |R|<atol, final iters=%3d",_monctx.iters+1);
            str=buff;
            MessagePrinter::PrintShortTxt(str);
        }
        else{
            snprintf(buffnew,68,"  SNES solver: iters=%3d,|R0|=%12.5e,|R|=%12.5e",_Iters+1,_monctx.rnorm0,_Rnorm);
            str=buffnew;
            MessagePrinter::PrintNormalTxt(str);
        }
        return true;
    }
    else if(_snesreason==SNES_CONVERGED_FNORM_RELATIVE){
        if(fectrlinfo.IsDepDebug){
            snprintf(buff,65,"  Converged for |R|<rtol*|R0|, final iters=%3d",_monctx.iters+1);
            str=buff;
            MessagePrinter::PrintShortTxt(str);
        }
        else{
            snprintf(buffnew,68,"  SNES solver: iters=%3d,|R0|=%12.5e,|R|=%12.5e",_Iters+1,_monctx.rnorm0,_Rnorm);
            str=buffnew;;
            MessagePrinter::PrintNormalTxt(str);
        }
        return true;
    }
    else if(_snesreason==SNES_CONVERGED_SNORM_RELATIVE){
        if(fectrlinfo.IsDepDebug){
            snprintf(buff,65,"  Converged for |delta x|<stol|x|, final iters=%3d",_monctx.iters+1);
            str=buff;
            MessagePrinter::PrintShortTxt(str);
        }
        else{
            snprintf(buffnew,68,"  SNES solver: iters=%3d,|R0|=%12.5e,|R|=%12.5e",_Iters+1,_monctx.rnorm0,_Rnorm);
            str=buffnew;
            MessagePrinter::PrintNormalTxt(str);
        }
        return true;
    }
    else{
        snprintf(buff,65,"  Divergent, SNES nonlinear solver failed, iters=%3d",_monctx.iters+1);
        str=buff;
        MessagePrinter::PrintShortTxt(str);
        return false;
    }
    return false;
}