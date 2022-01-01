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
        snprintf(buff,68,"  SNES solver: iters=%4d, |R|=%12.5e, |dU|=%12.5e",iters,rnorm,user->dunorm);
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
    user->_bcSystem->ApplyInitialBC(*user->_mesh,*user->_dofHandler,
                                   user->_fectrlinfo->t+user->_fectrlinfo->dt,// t in fectrlinfo is still the previous one
                                   U);

    // set the ctan array according to different time stepping method
    if(user->_fectrlinfo->_timesteppingtype==TimeSteppingType::STATIC){
        // for static analysis
        user->_fectrlinfo->ctan[0]=1.0;
        user->_fectrlinfo->ctan[1]=0.0;
        user->_fectrlinfo->ctan[2]=0.0;
        VecCopy(U,user->_solutionSystem->_Utemp);// in FormBulkFE, we always use Utemp for the calculation
    }
    else if(user->_fectrlinfo->_timesteppingtype==TimeSteppingType::BACKWARDEULER){
        user->_fectrlinfo->ctan[0]=1.0;
        user->_fectrlinfo->ctan[1]=1.0/user->_fectrlinfo->dt;
        user->_fectrlinfo->ctan[2]=0.0;
        // calculate the current velocity
        VecWAXPY(user->_solutionSystem->_V,-1.0,user->_solutionSystem->_U,U);//V=-Uold+Unew
        VecScale(user->_solutionSystem->_V,user->_fectrlinfo->ctan[1]);//V=V*1.0/dt
        VecCopy(U,user->_solutionSystem->_Utemp);// Used in FormBulkFE
    }
    else if(user->_fectrlinfo->_timesteppingtype==TimeSteppingType::CRANCKNICLSON){
        user->_fectrlinfo->ctan[0]=0.5;
        user->_fectrlinfo->ctan[1]=1.0/user->_fectrlinfo->dt;
        user->_fectrlinfo->ctan[2]=0.0;
        // calculate the current velocity
        VecWAXPY(user->_solutionSystem->_V,-1.0,user->_solutionSystem->_U,U);//V=-Uold+Unew
        VecScale(user->_solutionSystem->_V,user->_fectrlinfo->ctan[1]);//V=V*1.0/dt
        
        // for Cranck-Nicolson Utemp=0.5*(Unew+U)
        VecWAXPY(user->_solutionSystem->_Utemp,1.0,U,user->_solutionSystem->_U);// Utemp=Unew+Uold
        VecScale(user->_solutionSystem->_Utemp,0.5);// Utemp=0.5*(Unew+U)
    }
    else if(user->_fectrlinfo->_timesteppingtype==TimeSteppingType::BDF2){
        // the expression for BDF2 is:
        // [Un+2-(4/3)Un+1+(1/3)Un]/dt=(2/3)f(Un+1)

        if(user->_fectrlinfo->CurrentStep<2){
            user->_fectrlinfo->ctan[0]=1.0;
            user->_fectrlinfo->ctan[1]=1.0/user->_fectrlinfo->dt;
            user->_fectrlinfo->ctan[2]=0.0;
            // calculate the current velocity
            VecWAXPY(user->_solutionSystem->_V,-1.0,user->_solutionSystem->_U,U);//V=-Uold+Unew
            VecScale(user->_solutionSystem->_V,user->_fectrlinfo->ctan[1]);//V=V*1.0/dt
        }
        else{
            user->_fectrlinfo->ctan[0]=2.0/3.0;
            user->_fectrlinfo->ctan[1]=1.0/user->_fectrlinfo->dt;
            user->_fectrlinfo->ctan[2]=0.0;
            // calculate the current velocity
            VecWAXPY(user->_solutionSystem->_V,-4.0/3.0,user->_solutionSystem->_U,U);//V=Unew-(4/3)*U
            VecAXPY(user->_solutionSystem->_V,1.0/3.0,user->_solutionSystem->_Uold);// V=V+(1/3)Uold
            VecScale(user->_solutionSystem->_V,user->_fectrlinfo->ctan[1]);//V=V*1.0/dt
        }
        // for Utemp used in FormBulkFE
        VecCopy(U,user->_solutionSystem->_Utemp);
    }
    else{
        MessagePrinter::PrintErrorTxt("Unsupported time stepping method in nonlinear solver");
        MessagePrinter::AsFem_Exit();
    }

    user->_feSystem->FormBulkFE(FECalcType::ComputeResidual,
                        user->_fectrlinfo->t+user->_fectrlinfo->dt,
                        user->_fectrlinfo->dt,
                        user->_fectrlinfo->ctan,
                        *user->_mesh,*user->_dofHandler,*user->_fe,
                        *user->_elmtSystem,*user->_mateSystem,
                        *user->_solutionSystem,
                        user->_equationSystem->_AMATRIX,RHS);
    
    user->_bcSystem->SetBCPenaltyFactor(user->_feSystem->GetMaxAMatrixValue()*1.0e8);


    user->_bcSystem->ApplyBC(*user->_mesh,*user->_dofHandler,*user->_fe,
            FECalcType::ComputeResidual,
            user->_fectrlinfo->t+user->_fectrlinfo->dt,
            user->_fectrlinfo->ctan,
            U,user->_solutionSystem->_V,
            user->_equationSystem->_AMATRIX,RHS);

    return 0;
}

//***************************************************************
//*** here we setup the subroutine for jacobian 
//***************************************************************
PetscErrorCode ComputeJacobian(SNES snes,Vec U,Mat Jac,Mat B,void *ctx){
    AppCtx *user=(AppCtx*)ctx;
    
    user->_feSystem->ResetMaxAMatrixValue();
    user->_bcSystem->ApplyInitialBC(*user->_mesh,*user->_dofHandler,
                                   user->_fectrlinfo->t+user->_fectrlinfo->dt,
                                   U);

    // set the ctan array according to different time stepping method
    if(user->_fectrlinfo->_timesteppingtype==TimeSteppingType::STATIC){
        // for static analysis
        user->_fectrlinfo->ctan[0]=1.0;
        user->_fectrlinfo->ctan[1]=0.0;
        user->_fectrlinfo->ctan[2]=0.0;
    }
    else if(user->_fectrlinfo->_timesteppingtype==TimeSteppingType::BACKWARDEULER){
        user->_fectrlinfo->ctan[0]=1.0;
        user->_fectrlinfo->ctan[1]=1.0/user->_fectrlinfo->dt;
        user->_fectrlinfo->ctan[2]=0.0;
        // calculate the current velocity
        VecWAXPY(user->_solutionSystem->_V,-1.0,user->_solutionSystem->_U,U);//V=-Uold+Unew
        VecScale(user->_solutionSystem->_V,user->_fectrlinfo->ctan[1]);//V=V*1.0/dt
    }   
    else if(user->_fectrlinfo->_timesteppingtype==TimeSteppingType::CRANCKNICLSON){
        user->_fectrlinfo->ctan[0]=0.5;
        user->_fectrlinfo->ctan[1]=1.0/user->_fectrlinfo->dt;
        user->_fectrlinfo->ctan[2]=0.0;
        // calculate the current velocity
        VecWAXPY(user->_solutionSystem->_V,-1.0,user->_solutionSystem->_U,U);//V=-Uold+Unew
        VecScale(user->_solutionSystem->_V,user->_fectrlinfo->ctan[1]);//V=V*1.0/dt
           
        // for Cranck-Nicolson Utemp=0.5*(Unew+U)
        VecWAXPY(user->_solutionSystem->_Utemp,1.0,U,user->_solutionSystem->_U);// Utemp=Unew+Uold
        VecScale(user->_solutionSystem->_Utemp,0.5);// Utemp=0.5*(Unew+U)
    }   
    else if(user->_fectrlinfo->_timesteppingtype==TimeSteppingType::BDF2){
        // the expression for BDF2 is:
        // [Un+2-(4/3)Un+1+(1/3)Un]/dt=(2/3)f(Un+1)
   
        if(user->_fectrlinfo->CurrentStep<2){
            user->_fectrlinfo->ctan[0]=1.0;
            user->_fectrlinfo->ctan[1]=1.0/user->_fectrlinfo->dt;
            user->_fectrlinfo->ctan[2]=0.0;
            // calculate the current velocity
            VecWAXPY(user->_solutionSystem->_V,-1.0,user->_solutionSystem->_U,U);//V=-Uold+Unew
            VecScale(user->_solutionSystem->_V,user->_fectrlinfo->ctan[1]);//V=V*1.0/dt
        }   
        else{   
            user->_fectrlinfo->ctan[0]=2.0/3.0;
            user->_fectrlinfo->ctan[1]=1.0/user->_fectrlinfo->dt;
            user->_fectrlinfo->ctan[2]=0.0;
            // calculate the current velocity
            VecWAXPY(user->_solutionSystem->_V,-4.0/3.0,user->_solutionSystem->_U,U);//V=Unew-(4/3)*U
            VecAXPY(user->_solutionSystem->_V,1.0/3.0,user->_solutionSystem->_Uold);// V=V+(1/3)Uold
            VecScale(user->_solutionSystem->_V,user->_fectrlinfo->ctan[1]);//V=V*1.0/dt
        }   
        // for Utemp used in FormBulkFE
        VecCopy(U,user->_solutionSystem->_Utemp);
    }   
    else{   
        MessagePrinter::PrintErrorTxt("Unsupported time stepping method in nonlinear solver");
        MessagePrinter::AsFem_Exit();
    }   
   
   
    user->_feSystem->FormBulkFE(FECalcType::ComputeJacobian,
                        user->_fectrlinfo->t+user->_fectrlinfo->dt,
                        user->_fectrlinfo->dt,
                        user->_fectrlinfo->ctan,
                        *user->_mesh,*user->_dofHandler,*user->_fe,
                        *user->_elmtSystem,*user->_mateSystem,
                        *user->_solutionSystem,
                        B,user->_equationSystem->_RHS);
    
    user->_bcSystem->ApplyBC(*user->_mesh,*user->_dofHandler,*user->_fe,
                            FECalcType::ComputeJacobian,
                            user->_fectrlinfo->t+user->_fectrlinfo->dt,
                            user->_fectrlinfo->ctan,
                            U,user->_solutionSystem->_V,
                            B,user->_equationSystem->_RHS);

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
bool NonlinearSolver::Solve(Mesh &mesh,DofHandler &dofHandler,
                        ElmtSystem &elmtSystem,MateSystem &mateSystem,
                        BCSystem &bcSystem,
                        SolutionSystem &solutionSystem,EquationSystem &equationSystem,
                        FE &fe,FESystem &feSystem,
                        FEControlInfo &fectrlinfo){
    
    _appctx=AppCtx{&mesh,&dofHandler,
                   &bcSystem,
                   &elmtSystem,&mateSystem,
                   &solutionSystem,&equationSystem,
                   &fe,&feSystem,
                   &fectrlinfo
                   };
    
    _monctx=MonitorCtx{0.0,1.0,
            0.0,1.0,
            0.0,1.0,
            0,
            fectrlinfo.IsDepDebug};


    _appctx._bcSystem->ApplyInitialBC(*_appctx._mesh,*_appctx._dofHandler,1.0,_appctx._solutionSystem->_Unew);
    
    SNESSetFunction(_snes,_appctx._equationSystem->_RHS,ComputeResidual,&_appctx);

    SNESSetJacobian(_snes,_appctx._equationSystem->_AMATRIX,_appctx._equationSystem->_AMATRIX,ComputeJacobian,&_appctx);

    SNESMonitorSet(_snes,Monitor,&_monctx,0);

    SNESSetForceIteration(_snes,PETSC_TRUE);

    SNESSetFromOptions(_snes);
    
    SNESSolve(_snes,NULL,_appctx._solutionSystem->_Unew);
   
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
