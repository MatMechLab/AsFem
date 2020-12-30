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
//+++ Date   : 2020.12.30
//+++ Purpose: here we call the TS+SNES solver from PETSc to solve
//+++          our nonlinear system equation R(x,xdot,t)->0
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "TimeStepping/TimeStepping.h"


//***************************************************************
//*** here we define a monitor to print out the iteration info
//***************************************************************
PetscErrorCode MySNESMonitor(SNES snes,PetscInt iters,PetscReal rnorm,void* ctx){
    char buff[70];
    string str;
    TSAppCtx *user=(TSAppCtx*)ctx;
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
        snprintf(buff,70,"  SNES solver:iters=%3d,|R|=%11.4e,|dU|=%11.4e,dt=%7.2e",iters,rnorm,user->dunorm,user->dt);
        str=buff;
        MessagePrinter::PrintNormalTxt(str);
    }
    return 0;
}
//**********************************
PetscErrorCode MyTSMonitor(TS ts,PetscInt step,PetscReal time,Vec U,void *ctx){
    TSAppCtx *user=(TSAppCtx*)ctx;
    char buff[68];
    string str;
    double dt;

    TSGetTimeStep(ts,&dt);

    user->time=time;
    user->step=step;
    user->dt=dt;
    
    snprintf(buff,68,"Time step=%8d, time=%13.5e, dt=%13.5e",step,time,dt);
    str=buff;
    MessagePrinter::PrintNormalTxt(str);
    if(!user->IsDebug){
        snprintf(buff,68," SNES solver: iters=%3d,|R0|=%12.5e,|R|=%12.5e",user->iters+1,user->rnorm0,user->rnorm);
        str=buff;
        MessagePrinter::PrintNormalTxt(str);
    }

    if(step%user->_outputSystem.GetIntervalNum()==0){
        if(user->_fectrlinfo.IsProjection){
            user->_feSystem.FormBulkFE(FECalcType::Projection,time,dt,user->_fectrlinfo.ctan,
            user->_mesh,user->_dofHandler,user->_fe,user->_elmtSystem,user->_mateSystem,
            U,user->_solutionSystem._V,user->_solutionSystem._Hist,user->_solutionSystem._HistOld,
            user->_solutionSystem._Proj,user->_equationSystem._AMATRIX,user->_equationSystem._RHS);

            user->_outputSystem.WriteResultToFile(step,user->_mesh,user->_dofHandler,U,
            user->_solutionSystem.GetProjNumPerNode(),user->_solutionSystem.GetProjNameVec(),
            user->_solutionSystem._Proj);
        }
        else{
            user->_outputSystem.WriteResultToFile(step,user->_mesh,user->_dofHandler,U);
        }
        MessagePrinter::PrintNormalTxt("Write result to "+user->_outputSystem.GetOutputFileName());
        MessagePrinter::PrintDashLine();
    }

    if(user->IsAdaptive){
        if(user->iters+1<=user->OptiIters){
            dt=user->dt*user->GrowthFactor;
            if(dt>user->DtMax) dt=user->DtMax;
        }
        else{
            dt=user->dt*user->CutbackFactor;
            if(dt<user->DtMin) dt=user->DtMin;
        }
        TSSetTimeStep(ts,dt);
    }

    return 0;
}
//***************************************************************
//*** for our Residual and Jacobian calculation
//***************************************************************
PetscErrorCode ComputeIResidual(TS ts,PetscReal t,Vec U,Vec V,Vec RHS,void *ctx){
    TSAppCtx *user=(TSAppCtx*)ctx;

    TSGetTimeStep(ts,&user->dt);
    // apply the initial dirichlet boundary condition
    user->_bcSystem.ApplyInitialBC(user->_mesh,user->_dofHandler,t,U);

    user->_feSystem.FormBulkFE(FECalcType::ComputeResidual,
                        t,user->_fectrlinfo.dt,user->_fectrlinfo.ctan,
                        user->_mesh,user->_dofHandler,user->_fe,
                        user->_elmtSystem,user->_mateSystem,
                        U,V,
                        user->_solutionSystem._Hist,user->_solutionSystem._HistOld,
                        user->_solutionSystem._Proj,
                        user->_equationSystem._AMATRIX,RHS);
    
    user->_bcSystem.SetBCPenaltyFactor(user->_feSystem.GetMaxAMatrixValue()*1.0e8);

    user->_bcSystem.ApplyBC(user->_mesh,user->_dofHandler,user->_fe,
                    FECalcType::ComputeResidual,t,user->_fectrlinfo.ctan,U,
                    user->_equationSystem._AMATRIX,RHS);
    
    return 0;
}
//******************************************************
PetscErrorCode ComputeIJacobian(TS ts,PetscReal t,Vec U,Vec V,PetscReal s,Mat A,Mat B,void *ctx){
    TSAppCtx *user=(TSAppCtx*)ctx;
    int i;

    TSGetTimeStep(ts,&user->_fectrlinfo.dt);
    TSGetTimeStep(ts,&user->dt);
    user->_feSystem.ResetMaxAMatrixValue();// we reset the penalty factor
    user->_bcSystem.ApplyInitialBC(user->_mesh,user->_dofHandler,t,U);

    user->_fectrlinfo.ctan[0]=1.0;
    user->_fectrlinfo.ctan[1]=s;// dUdot/dU

    user->_feSystem.FormBulkFE(FECalcType::ComputeJacobian,
                        t,user->_fectrlinfo.dt,user->_fectrlinfo.ctan,
                        user->_mesh,user->_dofHandler,user->_fe,
                        user->_elmtSystem,user->_mateSystem,
                        U,V,
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
                    FECalcType::ComputeJacobian,t,user->_fectrlinfo.ctan,U,
                    A,user->_equationSystem._RHS);

    MatGetSize(B,&i,&i);

    return 0;
}

//*************************************************************************
bool TimeStepping::Solve(Mesh &mesh,DofHandler &dofHandler,
            ElmtSystem &elmtSystem,MateSystem &mateSystem,
            BCSystem &bcSystem,ICSystem &icSystem,
            SolutionSystem &solutionSystem,EquationSystem &equationSystem,
            FE &fe,FESystem &feSystem,
            OutputSystem &outputSystem,
            FEControlInfo &fectrlinfo){

    _appctx=TSAppCtx{mesh,dofHandler,
                   bcSystem,icSystem,
                   elmtSystem,mateSystem,
                   solutionSystem,equationSystem,
                   fe,feSystem,
                   outputSystem,
                   fectrlinfo,
                   //*****************
                   0.0,0.0,
                   0.0,0.0,
                   0.0,0.0,
                   0,
                   fectrlinfo.IsDebug,fectrlinfo.IsDepDebug,
                   0.0,0.0,0,
                    //********************************
                    _Adaptive,
                    _OptIters,
                    _GrowthFactor,_CutBackFactor,
                    _DtMin,_DtMax
                   };
    


    _appctx._icSystem.ApplyIC(_appctx._mesh,_appctx._dofHandler,_appctx._solutionSystem._Unew);
    _appctx._bcSystem.ApplyInitialBC(_appctx._mesh,_appctx._dofHandler,1.0,_appctx._solutionSystem._Unew);

    TSSetIFunction(_ts,_appctx._equationSystem._RHS,ComputeIResidual,&_appctx);
    TSSetIJacobian(_ts,_appctx._equationSystem._AMATRIX,_appctx._equationSystem._AMATRIX,ComputeIJacobian,&_appctx);
    
    TSMonitorSet(_ts,MyTSMonitor,&_appctx,NULL);
    SNESMonitorSet(_snes,MySNESMonitor,&_appctx,0);

    SNESSetForceIteration(_snes,PETSC_TRUE);

    SNESSetFromOptions(_snes);

    TSSetFromOptions(_ts);

    TSSolve(_ts,_appctx._solutionSystem._Unew);

    return true;
}