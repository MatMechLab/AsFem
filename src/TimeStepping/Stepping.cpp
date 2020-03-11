//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "TimeStepping/TimeStepping.h"

// PetscErrorCode Residual(TS ts,PetscReal t,Vec U,Vec Udot,Vec RHS,void *ctx){
//     TSAppCtx *user=(TSAppCtx *)ctx;
//     double dt;
//     double ctan[]={1.0,1.0};
//     TSGetTimeStep(ts,&dt);

//     user->_feSystem.ResetMaxAMatrixValue();
//     user->_bcSystem.ApplyInitialBC(user->_mesh,user->_dofHandler,t,U);
//     user->_feSystem.FormFE(3,t,dt,ctan,
//                         user->_mesh,user->_dofHandler,
//                         user->_fe,user->_elmtSystem,user->_mateSystem,U,Udot,
//                         user->_solution._Hist,user->_solution._HistOld,user->_solution._Proj,
//                         user->_equationSystem._AMATRIX,RHS);
    
//     user->_bcSystem.SetBCPenaltyFactor(user->_feSystem.GetMaxAMatrixValue());
//     user->_bcSystem.ApplyBC(user->_mesh,user->_dofHandler,user->_fe,
//                         t,ctan,user->_equationSystem._AMATRIX,RHS,U);
    
//     return 0;
// }

// //**************************************
// PetscErrorCode Jacobian(TS ts,PetscReal t,Vec U,Vec Udot,PetscReal s,Mat A,Mat B,void *ctx){
//     TSAppCtx *user=(TSAppCtx *)ctx;
//     double dt;
//     double ctan[]={1.0,1.0};
//     int i;
//     TSGetTimeStep(ts,&dt);
//     ctan[1]=s;


//     user->_bcSystem.ApplyInitialBC(user->_mesh,user->_dofHandler,t,U);

//     user->_feSystem.ResetMaxAMatrixValue();
//     user->_feSystem.FormFE(6,t,dt,ctan,
//             user->_mesh,user->_dofHandler,user->_fe,user->_elmtSystem,
//             user->_mateSystem,U,Udot,
//             user->_solution._Hist,user->_solution._HistOld,
//             user->_solution._Proj,A,user->_equationSystem._RHS);
    
    

//     user->_bcSystem.SetBCPenaltyFactor(user->_feSystem.GetMaxAMatrixValue());
//     user->_bcSystem.ApplyBC(user->_mesh,user->_dofHandler,user->_fe,
//                     t,ctan,A,user->_equationSystem._RHS,U);
    
//     MatScale(A,-1.0);
//     MatGetSize(B,&i,&i);
//     return 0;
    
// }
// //*******************************************
// //*** for related monitor
// //*******************************************
// PetscErrorCode MySNESMonitor(SNES snes,PetscInt iters,PetscReal rnorm,void* ctx){
//     TSAppCtx *user=(TSAppCtx*)ctx;
//     user->iters=iters;
//     user->rnorm=rnorm;
//     SNESGetSolutionNorm(snes,&user->dunorm);
//     user->enorm=rnorm*user->dunorm;
//     if(iters==0){
//         user->rnorm0=rnorm;
//         user->dunorm0=user->dunorm;
//         user->enorm0=user->enorm;
//     }
//     if(user->IsDebug){
//         if(user->IsDepDebug){
//             PetscPrintf(PETSC_COMM_WORLD,"***    SNES solver: iters=%3d , |R|=%14.6e                    ***\n",iters,rnorm);
//         }
//     }

//     return 0;
// }
// PetscErrorCode MyTSMonitor(TS ts,PetscInt step,PetscReal time,Vec U,void *ctx){
//     TSAppCtx *user=(TSAppCtx*)ctx;
//     user->step=step;
//     user->time=time;
//     double dt;

//     TSGetTimeStep(ts,&dt);
//     double ctan[]={1.0,1.0};
//     user->dt=dt;
//     TS2GetSolution(ts,&user->_solution._Uold,&user->_solution._V);
    

//     if(user->IsDebug){
//         if(user->IsDepDebug){
//             PetscPrintf(PETSC_COMM_WORLD,"***-Step=%6d, time=%13.5e,  dt=%13.5e                ***\n",step,time,user->dt);
//             PetscPrintf(PETSC_COMM_WORLD,"***              |R0|=%13.5e, |R|=%13.5e                ***\n",user->rnorm0,user->rnorm);
//         }
//         else{
//             PetscPrintf(PETSC_COMM_WORLD,"***-Step=%6d, time=%13.5e,  dt=%13.5e                ***\n",step,time,user->dt);
//             PetscPrintf(PETSC_COMM_WORLD,"***              |R0|=%13.5e, |R|=%13.5e, iters=%4d    ***\n",user->rnorm0,user->rnorm,user->iters);
//         }
//     }
    
//     if(step%user->interval==0&&step>-1){
//         if(user->IsProjection){
//             user->_feSystem.FormFE(9,time,user->dt,ctan,
//                         user->_mesh,user->_dofHandler,user->_fe,
//                         user->_elmtSystem,user->_mateSystem,
//                         U,user->_solution._V,
//                         user->_solution._Hist,user->_solution._HistOld,
//                         user->_solution._Proj,
//                         user->_equationSystem._AMATRIX,user->_equationSystem._RHS);
            
//             //const int &step,Mesh &mesh,DofHandler &dofHandler,const Vec &U,const int &nproj,const vector<string> &namelist,const Vec &Proj
//             user->_outputSystem.WriteResultToVTU(step,user->_mesh,user->_dofHandler,U,
//                                             user->_solution.GetProjNumPerNode(),
//                                             user->_solution.GetProjNameVec(),
//                                             user->_solution._Proj);
//         }
//         else{
//             user->_outputSystem.WriteResultToVTU(step,user->_mesh,user->_dofHandler,U);
//         }

//         PetscPrintf(PETSC_COMM_WORLD,"*** Write result to [%41s] !!!   ***\n",user->_outputSystem.GetVTUFileName().c_str());
        
//         // PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n");
//         // PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n");
//     }

//     //update history variables
//     user->_feSystem.FormFE(8,time,user->dt,ctan,
//         user->_mesh,user->_dofHandler,user->_fe,
//         user->_elmtSystem,user->_mateSystem,
//         U,user->_solution._V,
//         user->_solution._Hist,user->_solution._HistOld,
//         user->_solution._Proj,
//         user->_equationSystem._AMATRIX,user->_equationSystem._RHS);
//     VecCopy(user->_solution._Hist,user->_solution._HistOld);

//     // if(user->iters>5){
//     //     dt=dt*0.8;
//     //     TSSetTimeStep(ts,dt);       
//     // }
//     // else{
//     //     dt=dt*1.1;
//     //     if(dt>0.1){
//     //         dt=0.1;
//     //     }
//     //     TSSetTimeStep(ts,dt);
//     // }

//     return 0;
// }

void TimeStepping::Stepping(/*Mesh &mesh,DofHandler &dofHandler,
                BCSystem &bcSystem,ICSystem &icSystem,
                ElmtSystem &elmtSystem,MateSystem &mateSystem,
                EquationSystem &equationSystem,Solution &solution,
                FE &fe,FESystem &feSystem,
                OutputSystem &outputSystem,
                FeCtrlInfo &fectrl*/){
    // _appctx=TSAppCtx{mesh,dofHandler,
    //             bcSystem,icSystem,
    //             elmtSystem,mateSystem,
    //             solution,equationSystem,
    //             fe,feSystem,outputSystem,
    //             0.0,0.0,
    //             0.0,0.0,
    //             0.0,0.0,
    //             0.1,1.0e-5,
    //             0,
    //             1,
    //             1,
    //             true,
    //             false,
    //             false};
    // _appctx.IsDebug=fectrl.IsDebug;
    // _appctx.IsDepDebug=fectrl.IsDepDebug;
    // _appctx.interval=_interval;
    // _appctx.IsProjection=fectrl.IsProjection;

    // _appctx._icSystem.ApplyIC(_appctx._mesh,_appctx._dofHandler,_appctx._solution._Unew);
    // _appctx._bcSystem.ApplyInitialBC(_appctx._mesh,_appctx._dofHandler,_dt0,_appctx._solution._Unew);

    // TSSetIFunction(_ts,_appctx._solution._Unew,Residual,&_appctx);
    // TSSetIJacobian(_ts,_appctx._equationSystem._AMATRIX,_appctx._equationSystem._AMATRIX,
    //             Jacobian,&_appctx);
    
    // TSGetSNES(_ts,&_snes);
    // SNESMonitorSet(_snes,MySNESMonitor,&_appctx,0);

    // TSMonitorSet(_ts,MyTSMonitor,&_appctx,0);

    // TSSetSolution(_ts,_appctx._solution._Unew);

    // TS2SetSolution(_ts,_appctx._solution._Unew,_appctx._solution._V);

    // TSSetUp(_ts);

    // TSSolve(_ts,_appctx._solution._Unew);
}