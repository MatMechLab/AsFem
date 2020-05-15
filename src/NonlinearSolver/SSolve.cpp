//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "NonlinearSolver/NonlinearSolver.h"

PetscErrorCode Monitor(SNES snes,PetscInt iters,PetscReal rnorm,void* ctx){
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
        // PetscPrintf(PETSC_COMM_WORLD,"***    SNES solver: iters=%3d , |R|=%14.6e                    ***\n",iters,rnorm);
        PetscPrintf(PETSC_COMM_WORLD,"***    SNES solver: iters=%3d ,|R|=%14.6e,|dU|=%14.6e ***\n",iters,rnorm,user->dunorm);
    }
    // PetscPrintf(PETSC_COMM_WORLD,"***    SNES solver: iters=%3d ,|R|=%14.6e,|E|=%14.6e,|dU|=%14.6e***\n",iters,rnorm,user->enorm,user->dunorm);
    return 0;
}


PetscErrorCode FormResidual(SNES snes,Vec U,Vec RHS,void *ctx){
    AppCtx *user=(AppCtx*)ctx;
    int i;

    user->_bcSystem.ApplyInitialBC(user->_mesh,user->_dofHandler,user->_fectrlinfo.t,U);

    // calculate the current velocity
    VecWAXPY(user->_solution._V,-1.0,user->_solution._Uold,U);//V=-Uold+Unew
    VecScale(user->_solution._V,user->_fectrlinfo.ctan[1]);//V=V*1.0/dt

    if(user->_fectrlinfo.timesteppingtype==TimeSteppingType::BackWardEuler){
        // VecCopy(Vec x, Vec y); y = x
        VecCopy(U,user->_solution._Utemp);
    }
    else if(user->_fectrlinfo.timesteppingtype==TimeSteppingType::CrankNicolson){
        //VecWAXPY(Vec w,PetscScalar a,Vec x,Vec y); w = a ∗ x + y
        VecWAXPY(user->_solution._Utemp,1.0,user->_solution._Uold,U);
        VecScale(user->_solution._Utemp,user->_fectrlinfo.ctan[0]);
    }
    user->_feSystem.FormFE(3,user->_fectrlinfo.t,user->_fectrlinfo.dt,user->_fectrlinfo.ctan,
                           user->_mesh,user->_dofHandler,user->_fe,user->_elmtSystem,
                           user->_mateSystem,
                           user->_solution._Utemp,user->_solution._V,
                           user->_solution._Hist,user->_solution._HistOld,user->_solution._Proj,
                           user->_equationSystem._AMATRIX,RHS);
    
    user->_bcSystem.SetBCPenaltyFactor(user->_feSystem.GetMaxAMatrixValue()*1.0e8);

    user->_bcSystem.ApplyBC(user->_mesh,user->_dofHandler,user->_fe,
                        user->_fectrlinfo.t,user->_fectrlinfo.ctan,
                        user->_equationSystem._AMATRIX,RHS,user->_solution._Utemp);

    SNESGetMaxNonlinearStepFailures(snes,&i);

    return 0;
}

PetscErrorCode FormJacobian(SNES snes,Vec U,Mat A,Mat B,void *ctx){
    AppCtx *user=(AppCtx*)ctx;
    int i;

    user->_feSystem.ResetMaxAMatrixValue();

    user->_bcSystem.ApplyInitialBC(user->_mesh,user->_dofHandler,user->_fectrlinfo.t,U);

    //*** calculate the current velocity
    VecWAXPY(user->_solution._V,-1.0,user->_solution._Uold,U);//V=-Uold+Unew
    VecScale(user->_solution._V,user->_fectrlinfo.ctan[1]);//V=V*1.0/dt

    if(user->_fectrlinfo.timesteppingtype==TimeSteppingType::BackWardEuler){
        // VecCopy(Vec x, Vec y); y = x
        VecCopy(U,user->_solution._Utemp);
    }
    else if(user->_fectrlinfo.timesteppingtype==TimeSteppingType::CrankNicolson){
        //VecWAXPY(Vec w,PetscScalar a,Vec x,Vec y); w = a ∗ x + y
        VecWAXPY(user->_solution._Utemp,1.0,user->_solution._Uold,U);
        VecScale(user->_solution._Utemp,user->_fectrlinfo.ctan[0]);
    }

    // int rank;
    // chrono::high_resolution_clock::time_point mystart,myend;
    // MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    // if(rank==0){
    //     mystart=chrono::high_resolution_clock::now();
    // }

    user->_feSystem.FormFE(6,user->_fectrlinfo.t,user->_fectrlinfo.dt,user->_fectrlinfo.ctan,
                           user->_mesh,user->_dofHandler,user->_fe,user->_elmtSystem,
                           user->_mateSystem,
                           user->_solution._Utemp,user->_solution._V,
                           user->_solution._Hist,user->_solution._HistOld,user->_solution._Proj,
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
                        user->_fectrlinfo.t,user->_fectrlinfo.ctan,
                        A,user->_equationSystem._RHS,user->_solution._Utemp);

   
    // if(rank==0){
    //     myend=chrono::high_resolution_clock::now();
    // }
    // PetscPrintf(PETSC_COMM_WORLD,"*** System assemble using time=%14.6e [s] \n",
    // chrono::duration_cast<std::chrono::microseconds>(myend-mystart).count()/1.0e6);

    // MatScale(A,-1.0);
    MatGetSize(B,&i,&i);
    SNESGetMaxNonlinearStepFailures(snes,&i);

    return 0;
}


bool NonlinearSolver::SSolve(Mesh &mesh,DofHandler &dofHandler,
               ElmtSystem &elmtSystem,MateSystem &mateSystem,
               BCSystem &bcSystem,ICSystem &icSystem,
               Solution &solution,EquationSystem &equationSystem,
               FE &fe,FESystem &feSystem,
               FeCtrlInfo &fectrlinfo){

    _appctx=AppCtx{mesh,dofHandler,
                   bcSystem,icSystem,
                   elmtSystem,mateSystem,
                   solution,equationSystem,
                   fe,feSystem,
                   fectrlinfo
                   };

    _monctx=MonitorCtx{0.0,1.0,
            0.0,1.0,
            0.0,1.0,
            0,
            fectrlinfo.IsDepDebug};


    _appctx._bcSystem.ApplyInitialBC(_appctx._mesh,_appctx._dofHandler,1.0,_appctx._solution._Unew);
    
    SNESSetFunction(_snes,_appctx._equationSystem._RHS,FormResidual,&_appctx);

    SNESSetJacobian(_snes,_appctx._equationSystem._AMATRIX,_appctx._equationSystem._AMATRIX,FormJacobian,&_appctx);

    SNESMonitorSet(_snes,Monitor,&_monctx,0);

    SNESSetForceIteration(_snes,PETSC_TRUE);

    SNESSetFromOptions(_snes);

    SNESSolve(_snes,NULL,_appctx._solution._Unew);

    SNESGetConvergedReason(_snes,&_snesreason);
    
    _Iters=_monctx.iters;
    _Rnorm=_monctx.rnorm;
    _dUnorm=_monctx.dunorm;
    _Enorm=_monctx.enorm;

    if(_snesreason==SNES_CONVERGED_FNORM_ABS){
        if(fectrlinfo.IsDepDebug){
            PetscPrintf(PETSC_COMM_WORLD,"*** Convergent for |R|<atol, final iters=%3d                    !!!   ***\n",_monctx.iters+1);
        }
        return true;
    }
    else if(_snesreason==SNES_CONVERGED_FNORM_RELATIVE){
        if(fectrlinfo.IsDepDebug){
            PetscPrintf(PETSC_COMM_WORLD,"*** Convergent for |R|<rtol*|R0|, final iters=%3d               !!!   ***\n",_monctx.iters+1);
        }
        return true;
    }
    else if(_snesreason==SNES_CONVERGED_SNORM_RELATIVE){
        if(fectrlinfo.IsDepDebug){
            PetscPrintf(PETSC_COMM_WORLD,"*** Convergent for |delta x|<stol|x|, final iters=%3d           !!!   ***\n",_monctx.iters+1);
        }
        return true;
    }
    else{
        PetscPrintf(PETSC_COMM_WORLD,"*** Divergent, SNES nonlinear solver failed, iters=%3d          !!!   ***\n",_monctx.iters+1);
        // PetscPrintf(PETSC_COMM_WORLD,"*************************************************************************\n");
        return false;
    }

}

